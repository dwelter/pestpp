#include <random>
#include <iomanip>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "Ensemble.h"
#include "EnsembleSmoother.h"
#include "ObjectiveFunc.h"
#include "covariance.h"
#include "RedSVD-h.h";


PhiStats::PhiStats(const map<string, PhiComponets> &_phi_info)
{
	update(_phi_info);
}
void PhiStats::update(const map<string, PhiComponets> &_phi_info)
{
	phi_info = _phi_info;
	map<string, double> stats;
	Eigen::VectorXd phi_vec;
	size = phi_info.size();
	phi_vec.resize(size);
	PhiComponets pc;
	int i = 0;
	for (auto &pi : phi_info)
		phi_vec[i] = pi.second.meas;
	mean = phi_vec.mean();
	Eigen::VectorXd centered = (phi_vec - (Eigen::VectorXd::Ones(phi_vec.size()) * mean));
	double var = (centered.cwiseProduct(centered)).sum();
	std = sqrt(var / double(phi_vec.size()));
	max = phi_vec.maxCoeff();
	min = phi_vec.minCoeff();
}

void PhiStats::rec_report(ofstream &f_rec)
{
	f_rec << "ensemble size:               " << setw(20) << left << size << endl;
	f_rec << "ensemble mean:               " << setw(20) << left << mean << endl;
	f_rec << "ensemble standard deviation: " << setw(20) << left << std << endl;
	f_rec << "ensemble min:                " << setw(20) << left << min << endl;
	f_rec << "ensemble max:                " << setw(20) << left << max << endl;
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
		if (phi_info.find(rname) != phi_info.end())
			csv << phi_info[rname].meas;
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

	last_best_mean = 1.0E+30;
	last_best_std = 1.0e+30;

	lam_mults = pest_scenario.get_pestpp_options_ptr()->get_ies_lam_mults();
	if (lam_mults.size() == 0)
		lam_mults.push_back(1.0);

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
			throw_ies_error(string("error loading restat obs ensemble csv"));
		}
	}
	
	

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
}

void IterEnsembleSmoother::solve()
{
	iter++;
	stringstream ss;
	ss << "starting solve for iteration: " << iter;
	performance_log->log_event(ss.str());

	if (oe.shape().first < 2)
		throw_ies_error(string("less than 2 active realizations...cannot continue"));

	//build up some components for subsetting later
	vector<int> drop_idxs;
	int subset = 3; // eventually ++ies_subset
	if (false) //eventaully if ies_subset
	{
		for (int i = subset; i < pe.shape().first; i++)
			drop_idxs.push_back(i);
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
	//cout << obscov_inv_sqrt.rows() << ',' << obscov_inv_sqrt.cols() << endl;
	Eigen::MatrixXd obs_diff = scale * (obscov_inv_sqrt * diff);

	performance_log->log_event("calculate scaled par diff");
	diff = pe.get_eigen_mean_diff(vector<string>(), act_par_names).transpose();
	Eigen::MatrixXd par_diff = scale * diff;

#ifdef _DEBUG
	cout << "scaled_residual" << endl << scaled_residual << endl << endl;
	cout << "par_diff" << endl << par_diff << endl << endl;
	cout << "obs_diff" << endl << obs_diff << endl << endl;
#endif

	performance_log->log_event("SVD of obs diff");
	RedSVD::RedSVD<Eigen::MatrixXd> rsvd(obs_diff);
	
	Eigen::MatrixXd ivec, upgrade_1, s = rsvd.singularValues(), V = rsvd.matrixV(), Ut = rsvd.matrixU().transpose();
	Eigen::MatrixXd s2 = s.cwiseProduct(s);
	vector<ParameterEnsemble> pe_lams;
	for (auto &lam_mult : lam_mults)
	{
		ss.str("");
		double cur_lam = last_best_lam * lam_mult;
		ss << "starting calcs for lambda" << cur_lam;
		performance_log->log_event(ss.str());
		performance_log->log_event("form scaled identity matrix");
		ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();
		performance_log->log_event("calculate portion of upgrade_1 for localization");
		upgrade_1 = -1.0 * par_diff * V * s.asDiagonal() * ivec * Ut;
		
		//localization here...

		performance_log->log_event("apply residuals to upgrade_1");
		upgrade_1 = (upgrade_1 * scaled_residual).transpose();
		 
		ParameterEnsemble pe_lam = pe;
		pe_lam.set_eigen(*pe_lam.get_eigen_ptr() + upgrade_1);

#ifdef _DEBUG
		cout << "V:" << V.rows() << ',' << V.cols() << endl;
		cout << V << endl << endl;
		cout << "Ut" << Ut.rows() << ',' << Ut.cols() << endl;
		cout << Ut << endl << endl;
		cout << "s:" << s.size() << endl;
		cout << s << endl << endl;
		cout << "ivec:" << ivec.size() << endl;
		cout << ivec << endl << endl;
		cout << "par_diff:" << par_diff.rows() << ',' << par_diff.cols() << endl;
		cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;
		cout << upgrade_1 << endl << endl;
		cout << "pe_lam:" << pe_lam.shape().first << "," << pe.shape().second << endl;
		cout << *pe_lam.get_eigen_ptr() << endl << endl;
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

		pe_lams.push_back(pe_lam);
	}

	//queue up runs
	ObservationEnsemble oe_lam = oe, oe_lam_best = oe;//two copies
	//drop the last rows of oe_lam and oe_lam_best if we are subsetting
	if (false) //++ies_subset
	{
		oe_lam.drop_rows(drop_idxs);
		oe_lam_best.drop_rows(drop_idxs);
	}
	ParameterEnsemble pe_lam_best = pe; //copy
	vector<map<int, int>> real_run_ids_lams;
	PhiStats phistats;
	int best_idx = -1;
	double best_mean = 1.0e+30, best_std = 1.0e+30;
	for (auto &pe_lam : pe_lams)
	{
		if (false) //++ies_subset
		{
			//don't want to drop rows because we are going to need one of these whole later
			//prob need to pass through a vector of ints to the _pe.run_ensemble...
			//pe_lam.drop_rows(drop_idxs);
			
		}
		run_ensemble(pe_lam, oe_lam);
		phistats.update(get_phi_info(oe_lam));
		if (phistats.get_mean() < best_mean)
		{
			oe_lam_best = oe_lam;
			best_mean = phistats.get_mean();
			pe_lam_best = pe_lam;
		}
	}

	if (false)//subset stuff here
	{

	}

	//phi_stats = report_and_save();

	//make runs
	//run_mgr_ptr->run();

	//process runs to find "best" lambda
	//ObservationEnsemble oe_lam = oe; //copy
	//vector<int> failed_real_idxs;
	//for (int ilam = 0; ilam < real_run_ids_lams.size(); ilam++)
	//{
	//	map<int, int> rids = real_run_ids_lams[ilam];
	//	failed_real_idxs = oe_lam.update_from_runs(rids, run_mgr_ptr);
	//	if 
	//}


}


map<string,PhiComponets> IterEnsembleSmoother::get_phi_info(ObservationEnsemble &_oe)
{
	Observations obs = pest_scenario.get_ctl_observations();
	Parameters pars = pest_scenario.get_ctl_parameters();
	ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
	PhiComponets phi_comps;
	vector<string> par_reals = pe.get_real_names(), obs_reals = oe.get_real_names();
	vector<double> svec;
	Eigen::VectorXd evec;
	map<string, PhiComponets> phi_comps_ensemble;
	//for (auto &name : oe.get_real_names())
	for (int ireal=0;ireal<oe.shape().first;ireal++)
	{
		obs.update_without_clear(oe.get_var_names(),oe.get_real_vector(obs_reals[ireal]));
		pars.update_without_clear(pe.get_var_names(), pe.get_real_vector(par_reals[ireal]));
		phi_comps = obj_func.get_phi_comp(obs, pars, pest_scenario.get_regul_scheme_ptr());
		phi_comps_ensemble[obs_reals[ireal]] = phi_comps;
	}
	return phi_comps_ensemble;
}

//map<string, double> IterEnsembleSmoother::get_phi_vec_stats(map<string,PhiComponets> &phi_info)
//{
//	//map<string,PhiComponets> phi_info = get_phi_info();
//	map<string, double> stats;
//	Eigen::VectorXd phi_vec;
//	phi_vec.resize(oe.shape().first);
//	PhiComponets pc;
//	vector<string> obs_reals = oe.get_real_names();
//	for (int ireal = 0; ireal < oe.shape().first; ireal++)
//		phi_vec[ireal] = phi_info[obs_reals[ireal]].meas;
//	stats["mean"] = phi_vec.mean();
//	Eigen::VectorXd centered = (phi_vec - (Eigen::VectorXd::Ones(phi_vec.size()) * stats["mean"]));
//	double var = (centered.cwiseProduct(centered)).sum();
//
//	stats["standard deviation"] = sqrt(var / double(phi_vec.size()));
//	stats["max"] = phi_vec.maxCoeff();
//	stats["min"] = phi_vec.minCoeff();
//	return stats;
//}

PhiStats IterEnsembleSmoother::report_and_save()
{
	ofstream &f_rec = file_manager.rec_ofstream();
	f_rec << endl << "  ---  IterEnsembleSmoother iteration " << iter << " report  ---  " << endl;
	f_rec << "        number of active realizations: " << pe.shape().first << endl;
	f_rec << "        number of model runs: " << run_mgr_ptr->get_total_runs() << endl;
	
	map<string, PhiComponets> phi_info = get_phi_info(oe);
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

//EnsemblePair IterEnsembleSmoother::run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe)
//{
//
//	EnsemblePair epair(pe, oe);
//	stringstream ss;
//	ss << "queuing " << _pe.shape().first << " runs";
//	performance_log->log_event(ss.str());
//	try
//	{
//		epair.queue_runs(run_mgr_ptr);
//	}
//	catch (const exception &e)
//	{
//		stringstream ss;
//		ss << "error processing obs csv: " << e.what();
//		throw_ies_error(ss.str());
//	}
//	catch (...)
//	{
//		throw_ies_error(string("error queuing runs"));
//	}
//	performance_log->log_event("making runs");
//	try
//	{
//		epair.run(run_mgr_ptr);
//	}
//	catch (const exception &e)
//	{
//		stringstream ss;
//		ss << "error running ensemble: " << e.what();
//		throw_ies_error(ss.str());
//	}
//	catch (...)
//	{
//		throw_ies_error(string("error running ensemble"));
//	}
//	
//	performance_log->log_event("processing runs");
//	vector<int> failed_real_indices;
//	try
//	{
//		failed_real_indices = epair.process_runs(run_mgr_ptr);
//	}
//	catch (const exception &e)
//	{
//		stringstream ss;
//		ss << "error processing runs: " << e.what();
//		throw_ies_error(ss.str());
//	}
//	catch (...)
//	{
//		throw_ies_error(string("error processing runs"));
//	}
//
//	if (failed_real_indices.size() > 0)
//	{
//		stringstream ss;
//		vector<string> par_real_names = pe.get_real_names();
//		vector<string> obs_real_names = oe.get_real_names();
//		ss << "the following par:obs realization runs failed: ";
//		for (auto &i : failed_real_indices)
//		{
//			ss << par_real_names[i] << ":" << obs_real_names[i] << ',';
//		}
//		performance_log->log_event(ss.str());
//	}
//	return epair;
//}

vector<int> IterEnsembleSmoother::run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe)
{
	stringstream ss;
	ss << "queuing " << _pe.shape().first << " runs";
	performance_log->log_event(ss.str());
	
	map<int, int> real_run_ids;
	try
	{
		real_run_ids = _pe.add_runs(run_mgr_ptr);
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