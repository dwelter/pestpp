#include <random>
#include <iomanip>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "Ensemble.h"
#include "EnsembleSmoother.h"
#include "ObjectiveFunc.h"
#include "covariance.h"
#include "RedSVD-h.h"
#include "SVDPackage.h"


PhiStats::PhiStats(const map<string, double> &_phi_map)
{
	update(_phi_map);
}
void PhiStats::update(const map<string, double> &_phi_map)
{
	phi_map = _phi_map;
	map<string, double> stats;
	Eigen::VectorXd phi_vec;
	size = phi_map.size();
	phi_vec.resize(size);
	PhiComponets pc;
	int i = 0;
	for (auto &pi : phi_map)
	{
		phi_vec[i] = pi.second;
		i++;
	}
	mean = phi_vec.mean();
	Eigen::VectorXd centered = (phi_vec - (Eigen::VectorXd::Ones(phi_vec.size()) * mean));
	double var = (centered.cwiseProduct(centered)).sum();
	std = sqrt(var / double(phi_vec.size()));
	max = phi_vec.maxCoeff();
	min = phi_vec.minCoeff();
}

void PhiStats::rec_report(ofstream &f_rec)
{
	f_rec << "   ensemble size:                   " << setw(20) << left << size << endl;
	f_rec << "   ensemble phi mean:               " << setw(20) << left << mean << endl;
	f_rec << "   ensemble phi standard deviation: " << setw(20) << left << std << endl;
	f_rec << "   ensemble phi min:                " << setw(20) << left << min << endl;
	f_rec << "   ensemble phi max:                " << setw(20) << left << max << endl;
	cout << "   ensemble size:                   " << setw(20) << left << size << endl;
	cout << "   ensemble phi mean:               " << setw(20) << left << mean << endl;
	cout << "   ensemble phi standard deviation: " << setw(20) << left << std << endl;
	cout << "   ensemble phi min:                " << setw(20) << left << min << endl;
	cout << "   ensemble phi max:                " << setw(20) << left << max << endl;

}

void PhiStats::initialize_csv(ofstream &csv, const vector<string> &names)
{
	csv << "iteration,ensemble_size,phi_mean,phi_standard_deviation,phi_min,phi_max";
	for (auto &name : names)
		csv << ',' << name;
	csv << endl;
}

void PhiStats::csv_report(ofstream &csv, const int iter, const vector<string> &names)
{
	csv << iter << ',' << size << ',' << mean << ',' << std << ',' << min << ',' << max << ',';

	//vector<string>::iterator start = oe_real_names.begin(), end = oe_real_names.end();
	for (auto &rname : names)
	{
		if (phi_map.find(rname) != phi_map.end())
			csv << phi_map[rname];
		csv << ',';
	}
}


IterEnsembleSmoother::IterEnsembleSmoother(Pest &_pest_scenario, FileManager &_file_manager,
	OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log,
	RunManagerAbstract* _run_mgr_ptr) : pest_scenario(_pest_scenario), file_manager(_file_manager),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr)
{
	pe.set_pest_scenario(&pest_scenario);
	oe.set_pest_scenario(&pest_scenario);
	fphi_name = file_manager.get_base_filename() + ".phi.csv";
	fphi.open(fphi_name);
	if (!fphi.good())
		throw_ies_error("initialize: unable to open phi csv file " + fphi_name);
}

void IterEnsembleSmoother::throw_ies_error(string &message)
{
	file_manager.rec_ofstream() << "IterEnsembleSmoother error: " << message << endl;
	performance_log->log_event("IterEnsembleSmoother error: " + message);
	cout << "***" << endl << "IterEnsembleSmoother error: " << message << endl;
	file_manager.close_file("rec");
	performance_log->~PerformanceLog();
	throw runtime_error("IterEnsembleSmoother error: " + message);
}

void IterEnsembleSmoother::initialize()
{
	iter = 0;
	ofstream &frec = file_manager.rec_ofstream();
	last_best_mean = 1.0E+30;
	last_best_std = 1.0e+30;
	lambda_max = 1.0E+10;
	lambda_min = 1.0E-10;
	lam_mults = pest_scenario.get_pestpp_options_ptr()->get_ies_lam_mults();
	if (lam_mults.size() == 0)
		lam_mults.push_back(1.0);

	subset_size = pest_scenario.get_pestpp_options().get_ies_subset_size();
	

	performance_log->log_event("processing par csv");
	stringstream ss;
	try
	{
		pe.from_csv(pest_scenario.get_pestpp_options().get_ies_par_csv());
	}
	catch (const exception &e)
	{
		ss << "error processing par csv: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error processing par csv"));
	}

	performance_log->log_event("processing obs csv");
	try
	{
		oe.from_csv(pest_scenario.get_pestpp_options().get_ies_obs_csv());
	}
	catch (const exception &e)
	{
		ss << "error processing obs csv: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error processing obs csv"));
	}
	if (pe.shape().first != oe.shape().first)
	{
		ss.str("");
		ss << "parameter ensemble rows (" << pe.shape().first << ") not equal to observation ensemble rows (" << oe.shape().first << ")";
		throw_ies_error(ss.str());
	}

	
	if (subset_size > pe.shape().first)
	{
		//ss.str("");
		//ss << "subset_size " << subset_size << "greater than initial ensemble size " << pe.shape().first;
		//throw_ies_error(ss.str());
		use_subset = false;
	}
	else
	{
		frec << "  ---  using first " << subset_size << " realizations in ensemble lambda testing" << endl;
		cout << "  ---  using first " << subset_size << " realizations in ensemble lambda testing" << endl;
		use_subset = true;
	}
		
	oe_org_real_names = oe.get_real_names();
	pe_org_real_names = pe.get_real_names();

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	/*fphi << "iteration,num_active,";
	for (auto &s : phi_stat_names)
		fphi << s << ',';
	for (auto &rname : oe.get_real_names())
		fphi << rname << ',';
	fphi << endl;*/



	PhiStats::initialize_csv(fphi, oe_org_real_names);

	oe_base = oe; //copy
	//reorder this for later...
	oe_base.reorder(vector<string>(), act_obs_names);

	pe_base = pe;
	//reorder this for later
	pe_base.reorder(vector<string>(), act_par_names);

	performance_log->log_event("load obscov");
	Covariance obscov;
    obscov.from_observation_weights(pest_scenario);
	obscov = obscov.get(act_obs_names);
	obscov.inv_ip();
	//NEED SQRT operator to matrix
	obscov_inv_sqrt = obscov.get_matrix().diagonal().cwiseSqrt().asDiagonal();
	
	performance_log->log_event("load parcov");
	Covariance parcov;
	string parcov_filename = pest_scenario.get_pestpp_options().get_parcov_filename();
	if (parcov_filename.size() == 0)
		parcov.from_parameter_bounds(pest_scenario);
	else
	{
		string extension = parcov_filename.substr(parcov_filename.size() - 3, 3);
		pest_utils::upper_ip(extension);
		if (extension.compare("COV") == 0)
			parcov.from_ascii(parcov_filename);
		else if ((extension.compare("JCO") == 0) || (extension.compare("JCB") == 0))
			parcov.from_binary(parcov_filename);
		else if (extension.compare("UNC") == 0)
			parcov.from_uncertainty_file(parcov_filename);
		else
			throw_ies_error("unrecognized parcov_filename extension: " + extension);
	}
	parcov = parcov.get(act_par_names);
	parcov.inv_ip();
	Eigen::SparseMatrix<double> parcov_eigen = parcov.get_matrix();

	performance_log->log_event("calculate prior par diff");
	//no scaling...chen and oliver scale this...
	
	if (!pest_scenario.get_pestpp_options().get_ies_use_approx()) //eventually ies_use_approx
	{
		performance_log->log_event("calculating 'Am' matrix for full solution");
		double scale = (1.0 / (sqrt(double(pe.shape().first - 1))));
		Eigen::MatrixXd par_diff = scale * pe.get_eigen_mean_diff();
		RedSVD::RedSVD<Eigen::MatrixXd> rsvd(par_diff);
		Am = rsvd.matrixU() * rsvd.singularValues().inverse();
	}

	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();	
	//no restart
	if (obs_restart_csv.size() == 0)
	{
		performance_log->log_event("running initial ensemble");
		run_ensemble(pe, oe);
	}
	else
	{
		performance_log->log_event("restart IES with existing obs ensemble csv: " + obs_restart_csv);
		//ObservationEnsemble restart_oe(&pest_scenario);
		try
		{
			oe.from_csv(obs_restart_csv);
		}
		catch (const exception &e)
		{
			ss << "error loading restart obs ensemble csv:" << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error loading restart obs ensemble csv"));
		}
		//check that restart oe is in sync
		stringstream ss;
		vector<string> oe_real_names = oe.get_real_names(),oe_base_real_names = oe_base.get_real_names();
		vector<string>::const_iterator start, end;
		vector<string> missing;
		start = oe_base_real_names.begin();
		end = oe_base_real_names.end();
		for (auto &rname : oe_real_names)
			if (find(start, end, rname) == end)
				missing.push_back(rname);
		if (missing.size() > 0)
		{
			ss << "the following realization names were not found in the restart obs csv:";
			for (auto &m : missing)
				ss << m << ",";
			throw_ies_error(ss.str());

		}
		if (oe.shape().first < oe_base.shape().first) //maybe some runs failed...
		{		
			//find which realizations are missing and reorder oe_base, pe and pe_base
			vector<string> pe_real_names;
			start = oe_base_real_names.begin();
			end = oe_base_real_names.end();	
			vector<string>::const_iterator it;
			int iit;
			for (int i = 0; i < oe.shape().first; i++)
			{
				it = find(start, end, oe_real_names[i]);
				if (it != end)
				{
					iit = it - start;
					pe_real_names.push_back(pe_org_real_names[iit]);
				}
			}
			try
			{
				oe_base.reorder(oe_real_names, vector<string>());
			}
			catch (exception &e)
			{
				ss << "error reordering oe_base with restart oe:" << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error reordering oe_base with restart oe"));
			}
			try
			{
				pe.reorder(pe_real_names, vector<string>());
			}
			catch (exception &e)
			{
				ss << "error reordering pe with restart oe:" << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error reordering pe with restart oe"));
			}
			

		}
		else if (oe.shape().first > oe_base.shape().first) //something is wrong
		{
			ss << "restart oe has too many rows: " << oe.shape().first << " compared to oe_base: " << oe_base.shape().first;
			throw_ies_error(ss.str());
		}
	}

	map<string, double> phi_vec = get_phi_map(oe);
	PhiStats phistats = report_and_save();
	last_best_mean = phistats.get_mean();
	last_best_std = phistats.get_std();
	last_best_lam = pest_scenario.get_pestpp_options().get_ies_init_lam();
	if (last_best_lam <= 0.0)
	{
		double x = last_best_mean / (2.0 * double(oe.shape().second));
		last_best_lam = pow(10.0,(floor(log10(x))));
	}
	file_manager.rec_ofstream() << "  ---  current lambda: " << setw(20) << left << last_best_lam << endl;
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
}


void IterEnsembleSmoother::solve()
{
	iter++;
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();
	frec << endl << endl << endl << "  ---  starting solve for iteration: " << iter << endl;
	cout << endl << endl << endl << "  ---  starting solve for iteration: " << iter << endl;
	ss << "starting solve for iteration: " << iter;
	performance_log->log_event(ss.str());

	if (oe.shape().first < 2)
		throw_ies_error(string("less than 2 active realizations...cannot continue"));

	if ((use_subset) && (subset_size > pe.shape().first))
	{
		ss.str("");
		ss << "++ies_subset size (" << subset_size << ") greater than ensemble size (" << pe.shape().first << ")";
		frec << "  ---  " << ss.str() << endl;
		cout << "  ---  " << ss.str() << endl;
		frec << "  ...reducing ++ies_subset_size to " << pe.shape().first << endl;
		cout << "  ...reducing ++ies_subset_size to " << pe.shape().first << endl;
		subset_size = pe.shape().first;
	}

	double scale = (1.0 / (sqrt(double(oe.shape().first - 1))));
	
	performance_log->log_event("calculate residual matrix");
	//oe_base var_names should be ordered by act_obs_names, so only reorder real_names
	//oe should only include active realizations, so only reorder var_names
	Eigen::MatrixXd scaled_residual = obscov_inv_sqrt * (oe.get_eigen(vector<string>(), act_obs_names) - 
		oe_base.get_eigen(oe.get_real_names(), vector<string>())).transpose();
	
	performance_log->log_event("calculate scaled obs diff");
	Eigen::MatrixXd diff = oe.get_eigen_mean_diff(vector<string>(),act_obs_names).transpose();
	//cout << diff.rows() << ',' << diff.cols() << endl;
	cout << diff << endl;
	//cout << obscov_inv_sqrt.rows() << ',' << obscov_inv_sqrt.cols() << endl;
	Eigen::MatrixXd obs_diff = scale * (obscov_inv_sqrt * diff);

	performance_log->log_event("calculate scaled par diff");
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	diff = pe.get_eigen_mean_diff(vector<string>(), act_par_names).transpose();
	Eigen::MatrixXd par_diff = scale * diff;

//#ifdef _DEBUG
	cout << "scaled_residual" << endl << scaled_residual << endl << endl;
	cout << "par_diff" << endl << par_diff << endl << endl;
	cout << "obs_diff" << endl << obs_diff << endl << endl;
//#endif

	performance_log->log_event("SVD of obs diff");
	//RedSVD::RedSVD<Eigen::MatrixXd> rsvd1(obs_diff);	
	//Eigen::MatrixXd s1 = rsvd1.singularValues(), V1 = rsvd1.matrixV(), Ut1 = rsvd1.matrixU().transpose();
	Eigen::MatrixXd ivec, upgrade_1, s,V,Ut;

	SVD_REDSVD rsvd;
	rsvd.set_performance_log(performance_log);
	rsvd.solve_ip(obs_diff, s, Ut, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
	Ut.transposeInPlace();
	Eigen::MatrixXd s2 = s.cwiseProduct(s);
	vector<ParameterEnsemble> pe_lams;
	vector<double> lam_vals;
	for (auto &lam_mult : lam_mults)
	{
		ss.str("");
		double cur_lam = last_best_lam * lam_mult;
		ss << "starting calcs for lambda" << cur_lam;
		cout << "   ...starting calcs for lambda: " << cur_lam << endl;
		frec << "   ...starting calcs for lambda: " << cur_lam << endl;

		performance_log->log_event(ss.str());
		performance_log->log_event("form scaled identity matrix");
		ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();
		performance_log->log_event("calculate portion of upgrade_1 for localization");
		upgrade_1 = -1.0 * par_diff * V * s.asDiagonal() * ivec * Ut;
		
		//localization here...

		performance_log->log_event("apply residuals to upgrade_1");
		upgrade_1 = (upgrade_1 * scaled_residual).transpose();
		 
		ParameterEnsemble pe_lam = pe;//copy
		pe_lam.add_to_cols(upgrade_1, pe_base.get_var_names());

#ifdef _DEBUG
//		cout << "V:" << V.rows() << ',' << V.cols() << endl;
//		cout << V << endl << endl;
//		cout << "Ut" << Ut.rows() << ',' << Ut.cols() << endl;
//		cout << Ut << endl << endl;
//		cout << "s:" << s.size() << endl;
//		cout << s << endl << endl;
//		cout << "ivec:" << ivec.size() << endl;
//		cout << ivec << endl << endl;
//		cout << "par_diff:" << par_diff.rows() << ',' << par_diff.cols() << endl;
//		cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;
//		cout << upgrade_1 << endl << endl;
//		cout << "pe_lam:" << pe_lam.shape().first << "," << pe.shape().second << endl;
//		cout << *pe_lam.get_eigen_ptr() << endl << endl;
#endif

		if ((!pest_scenario.get_pestpp_options().get_ies_use_approx()) && (iter > 1))//eventually ies_use_approx
		{
			performance_log->log_event("calculating parameter correction (full solution)");

			performance_log->log_event("forming scaled par resid");
			Eigen::MatrixXd scaled_par_resid = pe.get_eigen(vector<string>(), act_par_names) - 
				pe_base.get_eigen(pe.get_real_names(), vector<string>());
			performance_log->log_event("forming x4");
			Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
			performance_log->log_event("forming x5");
			Eigen::MatrixXd x5 = Am * x4;
			performance_log->log_event("forming x6");
			Eigen::MatrixXd x6 = par_diff.transpose() * x5;
			performance_log->log_event("forming x7");
			Eigen::MatrixXd x7 = ivec * x6;
			performance_log->log_event("forming upgrade_2");
			Eigen::MatrixXd upgrade_2 = -1.0 * (par_diff * x7);
			pe_lam.set_eigen(*pe_lam.get_eigen_ptr() + upgrade_2.transpose());
			
		}

		pe_lam.enforce_bounds();
//#ifdef _DEBUG
//		cout << "pe_lam:" << pe_lam.shape().first << "," << pe.shape().second << endl;
//		cout << *pe_lam.get_eigen_ptr() << endl << endl;
//#endif
		pe_lams.push_back(pe_lam);
		lam_vals.push_back(cur_lam);
		frec << "   ...finished calcs for lambda: " << cur_lam << endl;
		cout << "   ...finished calcs for lambda: " << cur_lam << endl;

	}

	//queue up runs
	//ObservationEnsemble oe_lam = oe, oe_lam_best = oe;//two copies
	////drop the last rows of oe_lam and oe_lam_best if we are subsetting
	//if (subset_size < pe.shape().first) //++ies_subset
	//{
	//	vector<int> keep_idxs;
	//	for (int i = 0; i < subset_size; i++)
	//		keep_idxs.push_back(i);
	//	oe_lam.keep_rows(keep_idxs);
	//	oe_lam_best.keep_rows(keep_idxs);
	//}
	//ParameterEnsemble pe_lam_best = pe; //copy
	vector<map<int, int>> real_run_ids_lams;
	PhiStats phistats;
	int best_idx = -1;
	double best_mean = 1.0e+30, best_std = 1.0e+30;
	frec << "  ---  running lambda ensembles...  ---  " << endl;
	cout << "  ---  running lambda ensembles...  ---  " << endl;
	vector<ObservationEnsemble> oe_lams = run_lambda_ensembles(pe_lams);
	frec << "  ---  evaluting lambda ensemble results  --  " << endl;
	cout << "  ---  evaluting lambda ensemble results  --  " << endl;
	ObservationEnsemble oe_lam_best;
	for (int i=0;i<pe_lams.size();i++)
	{	
		//map<string, PhiComponets> temp = get_phi_info(oe_lams[i]);
		phistats.update(get_phi_map(oe_lams[i]));
		lam_test_report(lam_vals[i],phistats);
		if (phistats.get_mean() < best_mean)
		{
			oe_lam_best = oe_lams[i];
			best_mean = phistats.get_mean();
			best_idx = i;
		}
	}

	

	if ((use_subset) && (subset_size < pe.shape().first))//subset stuff here
	{

		ObservationEnsemble remaining_oe_lam = oe;//copy
		ParameterEnsemble remaining_pe_lam = pe_lams[best_idx];

		vector<int> real_idxs;
		vector<string> temp_idxs,org_pe_idxs,org_oe_idxs;
		int ireal = 0;
		for (int i = subset_size; i < pe.shape().first; i++)
		{
			real_idxs.push_back(i);
			ss.str("");
			ss << ireal;
			temp_idxs.push_back(ss.str());
			ireal++;
		}
		frec << "  ---  running remaining realizations from best lambda: " << lam_vals[best_idx] << endl;
		cout << "  ---  running remaining realizations from best lambda: " << lam_vals[best_idx] << endl;

		remaining_pe_lam.keep_rows(real_idxs);
		remaining_oe_lam.keep_rows(real_idxs);
		org_pe_idxs = remaining_pe_lam.get_real_names();
		org_oe_idxs = remaining_oe_lam.get_real_names();
		remaining_pe_lam.set_real_names(temp_idxs);
		remaining_oe_lam.set_real_names(temp_idxs);
		vector<int> fails = run_ensemble(remaining_pe_lam, remaining_oe_lam);
		//cout << *remaining_oe_lam.get_eigen_ptr() << endl;
		if (fails.size() > 0)
		{
			vector<string> new_pe_idxs, new_oe_idxs;
			vector<int>::iterator start = fails.begin(), end = fails.end();
			for (int i = 0; i < org_pe_idxs.size(); i++)
				if (find(start,end,i) == end)
				{
					new_pe_idxs.push_back(org_pe_idxs[i]);
					new_oe_idxs.push_back(org_oe_idxs[i]);
				}
			remaining_pe_lam.set_real_names(new_pe_idxs);
			remaining_oe_lam.set_real_names(new_oe_idxs);
		}
		else
		{
			remaining_pe_lam.set_real_names(org_pe_idxs);
			remaining_oe_lam.set_real_names(org_oe_idxs);
			
		}
		pe_lams[best_idx].drop_rows(real_idxs);
		pe_lams[best_idx].append_other_rows(remaining_pe_lam);
		oe_lam_best.append_other_rows(remaining_oe_lam);

	}
	//cout << oe_lam_best.get_var_names() << endl;
	//cout << *oe_lam_best.get_eigen_ptr() << endl;

	phistats.update(get_phi_map(oe_lam_best));
	frec << "  ---  last mean, current mean: " << last_best_mean << ',' << phistats.get_mean() << endl;
	frec << "  ---  last stdev, current stdev: " << last_best_std << ',' << phistats.get_std() << endl;
	cout << "  ---  last mean, current mean: " << last_best_mean << ',' << phistats.get_mean() << endl;
	cout << "  ---  last stdev, current stdev: " << last_best_std << ',' << phistats.get_std() << endl;

	if (phistats.get_mean() < last_best_mean * 1.1)
	{
		frec << "  ---  updating parameter ensemble  ---  " << endl;
		cout << "  ---  updating parameter ensemble  ---  " << endl;
		performance_log->log_event("updating parameter ensemble");
		last_best_mean = phistats.get_mean();
		pe = pe_lams[best_idx];
		oe = oe_lam_best;
		
		if (phistats.get_std() < last_best_std * 1.1)
		{
			double new_lam = lam_vals[best_idx] * 0.75;
			new_lam = (new_lam < lambda_min) ? lambda_min : new_lam;
			frec << "  ---  updating lambda from " << last_best_lam << " to " << new_lam << endl;
			cout << "  ---  updating lambda from " << last_best_lam << " to " << new_lam << endl;	
			last_best_lam = new_lam;
		}
		else
		{
			frec << "  ---  not updating lambda  --- " << endl;
		    cout << "  ---  not updating lambda  --- " << endl;
		}
	}

	else
	{
		frec << "  ---  not updating parameter ensemble  ---  " << endl;
		cout << "  ---  not updating parameter ensemble  ---  " << endl;

		double new_lam = last_best_mean * *max_element(lam_mults.begin(), lam_mults.end()) * 10.0;
		new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
		frec << " ---  increasing lambda from " << last_best_lam << " to " << new_lam << endl;
		cout << " ---  increasing lambda from " << last_best_lam << " to " << new_lam << endl;
		last_best_lam = new_lam;
	}
	report_and_save();
}

void IterEnsembleSmoother::lam_test_report(double lambda, PhiStats &phistats)
{
	ofstream &frec = file_manager.rec_ofstream();
	frec << "  ---  lambda: " << setw(10) << left << lambda;
	frec << ", mean: " << setw(10) << left << phistats.get_mean();
	frec << ", std: " << setw(10) << left << phistats.get_std() << endl;
	cout << "  ---  lambda: " << setw(10) << left << lambda;
	cout << ", mean: " << setw(10) << left << phistats.get_mean();
	cout << ", std: " << setw(10) << left << phistats.get_std() << endl;
}


map<string, double> IterEnsembleSmoother::get_phi_map(ObservationEnsemble &_oe)
{
	map<string, double> phi_map;
	ObservationInfo oinfo = pest_scenario.get_ctl_observation_info();
	Eigen::VectorXd oe_base_vec, oe_vec, q, diff;

	q.resize(act_obs_names.size());
	double w;
	//for (auto &oname : act_obs_names)
	for (int i = 0; i < act_obs_names.size(); i++)
		q(i) = oinfo.get_weight(act_obs_names[i]);
	vector<string> base_real_names = oe_base.get_real_names(),oe_real_names = _oe.get_real_names();
	vector<string>::iterator start = base_real_names.begin(), end = base_real_names.end();
	double phi;
	string rname;
	
	Eigen::MatrixXd oe_reals = _oe.get_eigen(vector<string>(), oe_base.get_var_names());
	//for (auto &rname : _oe.get_real_names())
	for (int i=0; i<_oe.shape().first;i++)
	{
		rname = oe_real_names[i];
		if (find(start, end, rname) == end)
			continue;
		oe_base_vec = oe_base.get_real_vector(rname);
		oe_vec = oe_reals.row(i);
		diff = (oe_vec - oe_base_vec).cwiseProduct(q);
		phi = (diff.cwiseProduct(diff)).sum();
		phi_map[rname] = phi;
	}
	return phi_map;
}


PhiStats IterEnsembleSmoother::report_and_save()
{
	ofstream &f_rec = file_manager.rec_ofstream();
	f_rec << endl << "  ---  IterEnsembleSmoother iteration " << iter << " report  ---  " << endl;
	f_rec << "   number of active realizations:  " << pe.shape().first << endl;
	f_rec << "   number of model runs:           " << run_mgr_ptr->get_total_runs() << endl;
	
	cout << endl << "  ---  IterEnsembleSmoother iteration " << iter << " report  ---  " << endl;
	cout << "   number of active realizations:   " << pe.shape().first << endl;
	cout << "   number of model runs:            " << run_mgr_ptr->get_total_runs() << endl;

	map<string, double> phi_info = get_phi_map(oe);
	PhiStats phistats(phi_info);
	
	phistats.rec_report(f_rec);

	stringstream ss;
	ss << file_manager.get_base_filename() << "." << iter << ".obs.csv";
	oe.to_csv(ss.str());
	f_rec << "      current obs ensemble saved to " << ss.str() << endl;
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
	pe.to_csv(ss.str());
	f_rec << "      current par ensemble saved to " << ss.str() << endl;
	f_rec << "      current ensemble phi info saved to " << fphi_name << endl;
	
	phistats.csv_report(fphi, iter, oe_org_real_names);

	
	return phistats;
}


vector<ObservationEnsemble> IterEnsembleSmoother::run_lambda_ensembles(vector<ParameterEnsemble> &pe_lams)
{
	stringstream ss;
	ss << "queuing " << pe_lams.size() << " ensembles";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	vector<int> subset_idxs;
	if ((use_subset) && (subset_size < pe_lams[0].shape().first))
	{
		for (int i = 0; i < subset_size; i++)
			subset_idxs.push_back(i);
	}
	vector<map<int, int>> real_run_ids_vec;
	for (auto &pe_lam : pe_lams)
	{
		try
		{
			real_run_ids_vec.push_back(pe_lam.add_runs(run_mgr_ptr,subset_idxs));
		}
		catch (const exception &e)
		{
			stringstream ss;
			ss << "run_ensemble() error queueing runs: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("run_ensembles() error queueing runs"));
		}
	}
	performance_log->log_event("making runs");
	try
	{
		
		run_mgr_ptr->run();
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error running ensembles: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error running ensembles"));
	}

	performance_log->log_event("processing runs");
	vector<int> failed_real_indices;
	vector<ObservationEnsemble> obs_lams;
	ObservationEnsemble _oe = oe;//copy
	if (subset_size < pe_lams[0].shape().first)
		_oe.keep_rows(subset_idxs);
	map<int, int> real_run_ids;
	//for (auto &real_run_ids : real_run_ids_vec)
	for (int i=0;i<pe_lams.size();i++)
	{
		
		real_run_ids = real_run_ids_vec[i];
		try
		{
			failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr);
		}
		catch (const exception &e)
		{
			stringstream ss;
			ss << "error processing runs: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error processing runs"));
		}

		if (failed_real_indices.size() > 0)
		{
			stringstream ss;
			vector<string> par_real_names = pe.get_real_names();
			vector<string> obs_real_names = oe.get_real_names();
			ss << "the following par:obs realization runs failed: ";
			for (auto &i : failed_real_indices)
			{
				ss << par_real_names[i] << ":" << obs_real_names[i] << ',';
			}
			performance_log->log_event(ss.str());
			performance_log->log_event("dropping failed realizations");
			_oe.drop_rows(failed_real_indices);
			
		}
		obs_lams.push_back(_oe);
	}
	return obs_lams;
}


vector<int> IterEnsembleSmoother::run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, vector<int> &real_idxs)
{
	stringstream ss;
	ss << "queuing " << _pe.shape().first << " runs";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	map<int, int> real_run_ids;
	try
	{
		real_run_ids = _pe.add_runs(run_mgr_ptr,real_idxs);
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "run_ensemble() error queueing runs: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("run_ensemble() error queueing runs"));
	}
	performance_log->log_event("making runs");
	try
	{
		run_mgr_ptr->run();
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error running ensemble: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error running ensemble"));
	}

	performance_log->log_event("processing runs");
	if (real_idxs.size() > 0)
	{
		_oe.keep_rows(real_idxs);
	}
	vector<int> failed_real_indices;
	try
	{
		failed_real_indices = _oe.update_from_runs(real_run_ids,run_mgr_ptr);
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error processing runs: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error processing runs"));
	}

	if (failed_real_indices.size() > 0)
	{
		stringstream ss;
		vector<string> par_real_names = pe.get_real_names();
		vector<string> obs_real_names = oe.get_real_names();
		ss << "the following par:obs realization runs failed: ";
		for (auto &i : failed_real_indices)
		{
			ss << par_real_names[i] << ":" << obs_real_names[i] << ',';
		}
		performance_log->log_event(ss.str());
		performance_log->log_event("dropping failed realizations");
		_pe.drop_rows(failed_real_indices);
		_oe.drop_rows(failed_real_indices);
	}
	return failed_real_indices;
}


void IterEnsembleSmoother::finalize()
{

}