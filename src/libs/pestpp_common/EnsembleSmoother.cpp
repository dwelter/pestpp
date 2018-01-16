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


PhiHandler::PhiHandler(Pest *_pest_scenario, FileManager *_file_manager,
	ObservationEnsemble *_oe_base, ParameterEnsemble *_pe_base,
	Covariance *_parcov_inv, double *_reg_factor)	
{
	pest_scenario = _pest_scenario;
	file_manager = _file_manager;
	oe_base = _oe_base;
	pe_base = _pe_base;
	reg_factor = _reg_factor;
	parcov_inv = _parcov_inv;
	oreal_names = oe_base->get_real_names();
	preal_names = pe_base->get_real_names();
	prepare_csv(file_manager->open_ofile_ext("phi.actual.csv"),oreal_names);
	prepare_csv(file_manager->open_ofile_ext("phi.meas.csv"),oreal_names);
	prepare_csv(file_manager->open_ofile_ext("phi.composite.csv"),oreal_names);
	prepare_csv(file_manager->open_ofile_ext("phi.regul.csv"),preal_names);

}


void PhiHandler::update(ObservationEnsemble & oe, ParameterEnsemble & pe)
{
	meas.clear();
	meas = calc_meas(oe);
	regul.clear();
	regul = calc_regul(pe);
	actual.clear();
	actual = calc_actual(oe);
	composite.clear();
	composite = calc_composite(meas,regul);

}

map<string, double>* PhiHandler::get_phi_map(PhiHandler::phiType &pt)
{
	switch (pt)
	{
	case PhiHandler::phiType::ACTUAL:
		return &actual;
	case PhiHandler::phiType::COMPOSITE:
		return &composite;
	case PhiHandler::phiType::MEAS:
		return &meas;
	case PhiHandler::phiType::REGUL:
		return &regul;
	}
}


double PhiHandler::get_mean(phiType pt)
{
	double mean = 0.0;
	map<string, double>* phi_map = get_phi_map(pt);
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (pi;pi != end; ++pi)
		mean = mean + pi->second;
	return mean / phi_map->size();
}

double PhiHandler::get_std(phiType pt)
{
	double mean = get_mean(pt);
	double var = 0.0;
	map<string, double>* phi_map = get_phi_map(pt);
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (pi; pi != end; ++pi)
		var = var + (pow(pi->second - mean,2));
	return sqrt(var/(phi_map->size()));
}

double PhiHandler::get_max(phiType pt)
{
	double mx = -1.0e+30;
	map<string, double>* phi_map = get_phi_map(pt);
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (pi; pi != end; ++pi)
		mx = (pi->second > mx) ? pi->second : mx;
	return mx;
}

double PhiHandler::get_min(phiType pt)
{
	double mn = 1.0e+30;
	map<string, double>* phi_map = get_phi_map(pt);
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (pi; pi != end; ++pi)
		mn = (pi->second < mn) ? pi->second : mn;
	return mn;
}

map<string, double> PhiHandler::get_summary_stats(PhiHandler::phiType pt)
{
	map<string, double> stats;
	stats["mean"] = get_mean(pt);
	stats["std"] = get_std(pt);
	stats["max"] = get_max(pt);
	stats["min"] = get_min(pt);
	return stats;
}

string PhiHandler::get_summary_string(PhiHandler::phiType pt)
{
	map<string, double> stats = get_summary_stats(pt);
	stringstream ss;
	string typ;
	switch (pt)
	{
	case PhiHandler::phiType::ACTUAL:
		typ = "actual";
		break;
	case PhiHandler::phiType::MEAS:
		typ = "measured";
		break;
	case PhiHandler::phiType::REGUL:
		typ = "regularization";
		break;
	case PhiHandler::phiType::COMPOSITE:
		typ = "composite";
		break;
	}

	ss << setw(15) << typ << setw(15) << stats["mean"] << setw(15) << stats["std"] << setw(15) << stats["min"] << setw(15) << stats["max"] << endl;
	return ss.str();
}


string PhiHandler::get_summary_header()
{
	stringstream ss;
	ss << setw(15) << "phi type" << setw(15) << "mean" << setw(15) << "std" << setw(15) << "min" << setw(15) << "max" << endl;
	return ss.str();

}


void PhiHandler::report()
{
	ofstream &f = file_manager->rec_ofstream();
	f << get_summary_header();
	cout << get_summary_header();
	string s = get_summary_string(PhiHandler::phiType::COMPOSITE);
	f << s;
	cout << s;
	s = get_summary_string(PhiHandler::phiType::MEAS);
	f << s;
	cout << s;
	s = get_summary_string(PhiHandler::phiType::REGUL);
	f << s;
	cout << s;
	s = get_summary_string(PhiHandler::phiType::ACTUAL);
	f << s;
	cout << s;
	f << endl << endl;
	f.flush();

}

void PhiHandler::write(int iter_num, int total_runs)
{
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.actual.csv"), phiType::ACTUAL,oreal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.meas.csv"), phiType::MEAS, oreal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.regul.csv"), phiType::REGUL, preal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.composite.csv"), phiType::COMPOSITE, oreal_names);
}

void PhiHandler::write_csv(int iter_num, int total_runs, ofstream &csv, phiType pt, vector<string> &names)
{
	map<string, double>* phi_map = get_phi_map(pt);
	map<string, double>::iterator pmi = phi_map->end();
	csv << iter_num << ',' << total_runs;
	map<string, double> stats = get_summary_stats(pt);
	csv << ',' << stats["mean"] << ',' << stats["std"] << ',' << stats["min"] << ',' << stats["max"];
	for (auto &name : names)
	{
		csv << ',';
		if (phi_map->find(name) != pmi)
			csv << phi_map->at(name);
	}
	csv << endl;
	csv.flush();
}

void PhiHandler::prepare_csv(ofstream & csv,vector<string> &names)
{
	csv << "iteration,total_runs,mean,standard_deviation,min,max";
	for (auto &name : names)
		csv << ',' << name;
	csv << endl;
}

map<string, double> PhiHandler::calc_meas(ObservationEnsemble & oe)
{
	map<string, double> phi_map;
	ObservationInfo oinfo = pest_scenario->get_ctl_observation_info();
	Eigen::VectorXd oe_base_vec, oe_vec, q, diff;
	vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	q.resize(act_obs_names.size());
	double w;
	//for (auto &oname : act_obs_names)
	for (int i = 0; i < act_obs_names.size(); i++)
		q(i) = oinfo.get_weight(act_obs_names[i]);
	vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
	vector<string>::iterator start = base_real_names.begin(), end = base_real_names.end();
	double phi;
	string rname;

	Eigen::MatrixXd oe_reals = oe.get_eigen(vector<string>(), oe_base->get_var_names());
	//for (auto &rname : _oe.get_real_names())
	//meas.clear();
	for (int i = 0; i<oe.shape().first; i++)
	{
		rname = oe_real_names[i];
		if (find(start, end, rname) == end)
			continue;
		oe_base_vec = oe_base->get_real_vector(rname);
		oe_vec = oe_reals.row(i);
		diff = (oe_vec - oe_base_vec).cwiseProduct(q);
		phi = (diff.cwiseProduct(diff)).sum();
		phi_map[rname] = phi;
	}
	return phi_map;
}

map<string, double> PhiHandler::calc_regul(ParameterEnsemble & pe)
{	
	map<string, double> phi_map;
	vector<string> real_names = pe.get_real_names();
	pe_base->transform_ip(ParameterEnsemble::transStatus::NUM);
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	Eigen::MatrixXd diff = pe.get_eigen() - pe_base->get_eigen(real_names, vector<string>());
	//cout << diff << endl;
	//cout << diff.rows() << ' ' << diff.cols() << endl;
	//cout << parcov_inv->e_ptr()->rows() << ' ' << parcov_inv->e_ptr()->cols() << endl;
	//cout << *parcov_inv->e_ptr() << endl;
	diff = diff * *parcov_inv->e_ptr() * diff.transpose();
	//cout << diff.rows() << ' ' << diff.cols() << endl;

	//regul.clear();
	for (int i = 0; i < real_names.size(); i++)
		phi_map[real_names[i]] = diff(i,i);
	return phi_map;
}

map<string, double> PhiHandler::calc_actual(ObservationEnsemble & oe)
{
	map<string, double> phi_map;
	ObservationInfo oinfo = pest_scenario->get_ctl_observation_info();
	Observations obs = pest_scenario->get_ctl_observations();
	Eigen::VectorXd obs_val_vec, oe_vec, q, diff;
	vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	q.resize(act_obs_names.size());
	obs_val_vec.resize(act_obs_names.size());
	double w;
	//for (auto &oname : act_obs_names)
	for (int i = 0; i < act_obs_names.size(); i++)
	{
		q(i) = oinfo.get_weight(act_obs_names[i]);
		obs_val_vec(i) = obs[act_obs_names[i]];
	}
		
	vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
	vector<string>::iterator start = base_real_names.begin(), end = base_real_names.end();
	double phi;
	string rname;

	Eigen::MatrixXd oe_reals = oe.get_eigen(vector<string>(), oe_base->get_var_names());
	//for (auto &rname : _oe.get_real_names())
	//meas.clear();
	for (int i = 0; i<oe.shape().first; i++)
	{
		rname = oe_real_names[i];
		if (find(start, end, rname) == end)
			continue;
		oe_vec = oe_reals.row(i);
		diff = (oe_vec - obs_val_vec).cwiseProduct(q);
		phi = (diff.cwiseProduct(diff)).sum();
		phi_map[rname] = phi;
	}
	return phi_map;
}


map<string, double> PhiHandler::calc_composite(map<string, double> &_meas, map<string, double> &_regul)
{
	map<string, double> phi_map;
	//composite.clear();
	string prn, orn;
	double reg, mea;
	map<string, double>::iterator meas_end = _meas.end(), regul_end = _regul.end();
	for (int i = 0; i < oe_base->shape().first; i++)
	{
		prn = preal_names[i];
		orn = oreal_names[i];
		if (meas.find(orn) != meas_end)
		{
			mea = _meas[orn];
			reg = _regul[prn];
			phi_map[orn] = mea + (reg * *reg_factor);

		}
	}
	return phi_map;
}	




IterEnsembleSmoother::IterEnsembleSmoother(Pest &_pest_scenario, FileManager &_file_manager,
	OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log,
	RunManagerAbstract* _run_mgr_ptr) : pest_scenario(_pest_scenario), file_manager(_file_manager),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr)
{
	pe.set_pest_scenario(&pest_scenario);
	oe.set_pest_scenario(&pest_scenario);
	//fphi_name = file_manager.get_base_filename() + ".phi.csv";
	//fphi.open(fphi_name);
	//if (!fphi.good())
	//	throw_ies_error("initialize: unable to open phi csv file " + fphi_name);
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
	cout << "  ---  initializing ---  " << endl << endl;
	ies_save_mat = true;
	iter = 0;
	ofstream &frec = file_manager.rec_ofstream();
	last_best_mean = 1.0E+30;
	last_best_std = 1.0e+30;
	lambda_max = 1.0E+10;
	lambda_min = 1.0E-10;
	lam_mults = pest_scenario.get_pestpp_options_ptr()->get_ies_lam_mults();
	if (lam_mults.size() == 0)
		lam_mults.push_back(1.0);
	cout << "...using lambda multipliers: ";
	for (auto &lm : lam_mults)
		cout << lm << ' , ';
	cout << endl;

	subset_size = pest_scenario.get_pestpp_options().get_ies_subset_size();
	reg_factor = pest_scenario.get_pestpp_options().get_ies_reg_factor();
	cout << "...using reg_factor:" << reg_factor << endl;

	string par_csv = pest_scenario.get_pestpp_options().get_ies_par_csv();
	performance_log->log_event("processing par csv "+par_csv);
	cout << "  ---  loading par ensemble from csv file " << par_csv << endl << endl;
	stringstream ss;
	try
	{
		pe.from_csv(par_csv);
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
	string obs_csv = pest_scenario.get_pestpp_options().get_ies_obs_csv();
	performance_log->log_event("processing obs csv "+obs_csv);
	cout << "  ---  loading obs ensemble from csv file " << obs_csv << endl << endl;

	try
	{
		oe.from_csv(obs_csv);
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

	//PhiStats::initialize_csv(fphi, oe_org_real_names);

	oe_base = oe; //copy
	//reorder this for later...
	oe_base.reorder(vector<string>(), act_obs_names);

	pe_base = pe;
	//reorder this for later
	pe_base.reorder(vector<string>(), act_par_names);

	performance_log->log_event("load obscov");
	cout << " --- initializing obs cov from obnservation weights  ---  " << endl << endl;
	Covariance obscov;
    obscov.from_observation_weights(pest_scenario);
	obscov = obscov.get(act_obs_names);
	obscov.inv_ip();
	//NEED SQRT operator to matrix
	obscov_inv_sqrt = obscov.get_matrix().diagonal().cwiseSqrt().asDiagonal();
	
	performance_log->log_event("load parcov");
	 
	string parcov_filename = pest_scenario.get_pestpp_options().get_parcov_filename();
	if (parcov_filename.size() == 0)
	{
		cout << "  ---  initializing par cov from parameter bounds  ---  " << endl << endl;
		parcov_inv.from_parameter_bounds(pest_scenario);
	}
	else
	{
		cout << "  ---  initializing par cov from file " << parcov_filename << endl << endl;
		string extension = parcov_filename.substr(parcov_filename.size() - 3, 3);
		pest_utils::upper_ip(extension);
		if (extension.compare("COV") == 0)
			parcov_inv.from_ascii(parcov_filename);
		else if ((extension.compare("JCO") == 0) || (extension.compare("JCB") == 0))
			parcov_inv.from_binary(parcov_filename);
		else if (extension.compare("UNC") == 0)
			parcov_inv.from_uncertainty_file(parcov_filename);
		else
			throw_ies_error("unrecognized parcov_filename extension: " + extension);
	}
	performance_log->log_event("inverting parcov");
	parcov_inv = parcov_inv.get(act_par_names);
	parcov_inv.inv_ip();
	//need this here for Am calcs...
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	pe_base.transform_ip(ParameterEnsemble::transStatus::NUM);
	
	
	if (!pest_scenario.get_pestpp_options().get_ies_use_approx()) //eventually ies_use_approx
	{
		cout << "  ---  using full (MAP) update solution  ---  " << endl << endl;
		performance_log->log_event("calculating 'Am' matrix for full solution");
		cout << "  forming Am matrix" << endl;
		double scale = (1.0 / (sqrt(double(pe.shape().first - 1))));
		Eigen::MatrixXd par_diff = scale * pe.get_eigen_mean_diff();
		par_diff.transposeInPlace();
		if (ies_save_mat)
			save_mat("prior_par_diff.dat",par_diff);
		//RedSVD::RedSVD<Eigen::MatrixXd> rsvd(par_diff);
		//Am = rsvd.matrixU() * rsvd.singularValues().inverse();
		
		Eigen::MatrixXd ivec, upgrade_1, s, V, U, st;
		SVD_REDSVD rsvd;
		rsvd.set_performance_log(performance_log);

		rsvd.solve_ip(par_diff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
		//rsvd.solve_ip(par_diff, s, U, V);
		//cout << "U " << U.rows() << ',' << U.cols() << endl;
		//cout << "s " << s.rows() << ',' << s.cols() << endl;
		//cout << s << endl;
		Eigen::MatrixXd temp = s.asDiagonal();
		//cout << temp << endl;
		Eigen::MatrixXd temp2 = temp.inverse();
		//cout << "U " << U.rows() << ',' << U.cols() << endl;
		//cout << "s " << s.rows() << ',' << s.cols() << endl;
		//cout << temp2 << endl;
		//cout << U << endl;

		Am = U * temp;
		//cout << Am << endl;
		//cout << "Am " << Am.rows() << ',' << Am.cols() << endl;
		if (ies_save_mat)
		{
			cout << "Am:" << Am.rows() << ',' << Am.cols() << endl;
			/*file_manager.open_ofile_ext("am.dat") << Am << endl;
			file_manager.open_ofile_ext("Am_U.dat") << U << endl;
			file_manager.open_ofile_ext("Am_s_inv.dat") << temp2 << endl;*/
			save_mat("am.dat", Am);
			save_mat("am_u.dat", U);
			save_mat("am_v.dat", V);
			save_mat("am_s_inv.dat", temp2);



		}
			
	}

	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();	
	//no restart
	if (obs_restart_csv.size() == 0)
	{
		performance_log->log_event("running initial ensemble");
		cout << "  ---  running initial ensemble of " << oe.shape().first << " realizations" << endl << endl;
		run_ensemble(pe, oe);
	}
	else
	{
		performance_log->log_event("restart IES with existing obs ensemble csv: " + obs_restart_csv);
		cout << "  ---  restarting with existing obs csv " << obs_restart_csv << endl;
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
			cout << "  ---  shape mismatch detected with restart obs ensemble...checking for compatibility..." << endl;
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
	performance_log->log_event("calc initial phi");
	ph = PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov_inv, &reg_factor);
	ph.update(oe, pe);
	frec << endl <<endl << "  ---  initial phi summary ---  " << endl;
	cout << endl <<endl << "  ---  initial phi summary ---  " << endl;
	ph.report();
	ph.write(0, run_mgr_ptr->get_total_runs());
	
	last_best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
	last_best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
	last_best_lam = pest_scenario.get_pestpp_options().get_ies_init_lam();
	if (last_best_lam <= 0.0)
	{
		double x = last_best_mean / (2.0 * double(oe.shape().second));
		last_best_lam = pow(10.0,(floor(log10(x))));
	}
	
	file_manager.rec_ofstream() << "  ---  current lambda: " << setw(20) << left << last_best_lam << endl;
	cout << "  ---  current lambda: " << setw(20) << left << last_best_lam << endl;
	
	//pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	cout << "  ---  initialization complete  --- " << endl << endl << endl;
}


void IterEnsembleSmoother::save_mat(string prefix, Eigen::MatrixXd &mat)
{
	stringstream ss;
	ss << iter << '.' << prefix;
	ofstream &f = file_manager.open_ofile_ext(ss.str());
	f << mat << endl;
	f.close();

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
	cout << "...calculating residual matrix" << endl;
	Eigen::MatrixXd scaled_residual = obscov_inv_sqrt * (oe.get_eigen(vector<string>(), act_obs_names) - 
		oe_base.get_eigen(oe.get_real_names(), vector<string>())).transpose();
	if (ies_save_mat)
	{
		cout << "scaled_residual:" << scaled_residual.rows() << ',' << scaled_residual.cols() << endl;
		save_mat("scaled_residual", scaled_residual);
		//file_manager.open_ofile_ext("scaled_residual.dat") << scaled_residual << endl;
	}
		
	performance_log->log_event("calculate scaled obs diff");
	cout << "...calculating obs diff matrix" << endl;
	Eigen::MatrixXd diff = oe.get_eigen_mean_diff(vector<string>(),act_obs_names).transpose();
	//cout << diff.rows() << ',' << diff.cols() << endl;
	//cout << diff << endl;
	//cout << obscov_inv_sqrt.rows() << ',' << obscov_inv_sqrt.cols() << endl;
	Eigen::MatrixXd obs_diff = scale * (obscov_inv_sqrt * diff);
	if (ies_save_mat)
	{
		cout << "obs_diff" << obs_diff.rows() << ',' << obs_diff.cols() << endl;
		//file_manager.open_ofile_ext("obs_diff.dat") << obs_diff << endl;
		save_mat("obs_diff.dat", obs_diff);
	}

	performance_log->log_event("calculate scaled par diff");
	cout << "...calculating par diff matrix" << endl;

	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	diff = pe.get_eigen_mean_diff(vector<string>(), act_par_names).transpose();
	Eigen::MatrixXd par_diff = scale * diff;
	if (ies_save_mat)
	{
		cout << "par_diff:" << par_diff.rows() << ',' << par_diff.cols() << endl;
		save_mat("par_diff.dat", par_diff);
		//file_manager.open_ofile_ext("par_diff.dat") << par_diff << endl;
	}

//#ifdef _DEBUG
	//cout << "scaled_residual" << endl << scaled_residual << endl << endl;
	//cout << "par_diff" << endl << par_diff << endl << endl;
	//cout << "obs_diff" << endl << obs_diff << endl << endl;
//#endif

	performance_log->log_event("SVD of obs diff");
	cout << "...calculating SVD of obs diff matrix" << endl;

	//RedSVD::RedSVD<Eigen::MatrixXd> rsvd1(obs_diff);	
	//Eigen::MatrixXd s1 = rsvd1.singularValues(), V1 = rsvd1.matrixV(), Ut1 = rsvd1.matrixU().transpose();
	Eigen::MatrixXd ivec, upgrade_1, s,V,Ut;

	SVD_REDSVD rsvd;
	rsvd.set_performance_log(performance_log);
	rsvd.solve_ip(obs_diff, s, Ut, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
	Ut.transposeInPlace();
	//cout << "V" << V << endl;
	//cout << "Ut" << endl << Ut << endl;
	Eigen::MatrixXd s2 = s.cwiseProduct(s);
	if (ies_save_mat)
	{
		cout << "s2:" << s2.rows() << ',' << s2.cols() << endl;
		//file_manager.open_ofile_ext("s2.dat") << s2 << endl;
		cout << "Ut:" << Ut.rows() << ',' << Ut.cols() << endl;
		cout << "V:" << V.rows() << ',' << V.cols() << endl;
		//file_manager.open_ofile_ext("Ut.dat") << Ut << endl;
		save_mat("ut.dat", Ut);
		save_mat("s2.dat", s2);
	}

	vector<ParameterEnsemble> pe_lams;
	vector<double> lam_vals;
	for (auto &lam_mult : lam_mults)
	{
		ss.str("");
		double cur_lam = last_best_lam * lam_mult;
		ss << "starting calcs for lambda" << cur_lam;
		//cout << "   ...starting calcs for lambda: " << cur_lam << endl;
		frec << "   ...starting calcs for lambda: " << cur_lam << endl;

		performance_log->log_event(ss.str());
		performance_log->log_event("form scaled identity matrix");
		cout << "...calculating scaled identity matrix" << endl;

		ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();
		if (ies_save_mat)
		{
			cout << "ivec:" << ivec.rows() << ',' << ivec.cols() << endl;
			//file_manager.open_ofile_ext("ivec.dat") << ivec << endl;
			save_mat("ivec.dat", ivec);
		}
		
		performance_log->log_event("calculate portion of upgrade_1 for localization");
		cout << "...calculating portion of upgrade_1 matrix for localization" << endl;

		upgrade_1 = -1.0 * par_diff;
		cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;

		cout << "..times V" << endl;
		upgrade_1 = upgrade_1 *  V;
		cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;

		cout << "..times s" << endl;
		upgrade_1 = upgrade_1 * s.asDiagonal();
		cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;

		cout << "..times ivec" << endl;
		upgrade_1 = upgrade_1 * ivec;
		cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;

		cout << "..times Ut" << endl;
		upgrade_1 = upgrade_1 * Ut;

		
		//localization here...

		performance_log->log_event("apply residuals to upgrade_1");
		upgrade_1 = (upgrade_1 * scaled_residual).transpose();
		if (ies_save_mat)
		{
			cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;
			//file_manager.open_ofile_ext("upgrade_1.dat") << upgrade_1 << endl;
			save_mat("upgrade_1.dat", upgrade_1);
		}
		//cout << "upgrade_1" << endl << upgrade_1 << endl;
		ParameterEnsemble pe_lam = pe;//copy
		pe_lam.add_to_cols(upgrade_1, pe_base.get_var_names());


		if ((!pest_scenario.get_pestpp_options().get_ies_use_approx()) && (iter > 1))
		{
			performance_log->log_event("calculating parameter correction (full solution)");

			performance_log->log_event("forming scaled par resid");
			Eigen::MatrixXd scaled_par_resid = pe.get_eigen(vector<string>(), act_par_names) - 
				pe_base.get_eigen(pe.get_real_names(), vector<string>());
			scaled_par_resid.transposeInPlace();
			
			performance_log->log_event("forming x4");
			if (ies_save_mat)
			{
				cout << "scaled_par_resid:" << scaled_par_resid.rows() << ',' << scaled_par_resid.cols() << endl;
				//file_manager.open_ofile_ext("scaled_par_resid.dat") << scaled_par_resid << endl;
				save_mat("scaled_par_resid.dat", scaled_par_resid);
			}
			Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
			if (ies_save_mat)
			{
				cout << "x4:" << x4.rows() << ',' << x4.cols() << endl;
				//file_manager.open_ofile_ext("x4.dat") << x4 << endl;
				save_mat("x4.dat", x4);
			}


			performance_log->log_event("forming x5");
			Eigen::MatrixXd x5 = Am * x4;
			if (ies_save_mat)
			{
				cout << "x5:" << x5.rows() << ',' << x5.cols() << endl;
				//file_manager.open_ofile_ext("x5.dat") << x5 << endl;
				save_mat("x5.dat", x5);
			}

			
			performance_log->log_event("forming x6");
			Eigen::MatrixXd x6 = par_diff.transpose() * x5;
			if (ies_save_mat)
			{
				cout << "x6:" << x6.rows() << ',' << x6.cols() << endl;
				//file_manager.open_ofile_ext("x6.dat") << x6 << endl;
				save_mat("x6.dat", x6);
			}
			
			
			performance_log->log_event("forming x7");
			if (ies_save_mat)
			{
				cout << "V:" << V.rows() << ',' << V.cols() << endl;
				//file_manager.open_ofile_ext("V.dat") << V << endl;
				save_mat("V.dat", V);
			}
			Eigen::MatrixXd x7 = V * ivec *V.transpose() * x6;
			if (ies_save_mat)
			{
				cout << "x7:" << x7.rows() << ',' << x7.cols() << endl;
				//file_manager.open_ofile_ext("x7.dat") << x7 << endl;
				save_mat("x7.dat", x7);
			}


			performance_log->log_event("forming upgrade_2");
			Eigen::MatrixXd upgrade_2 = -1.0 * (par_diff * x7);
			if (ies_save_mat)
			{
				cout << "upgrade_2:" << upgrade_2.rows() << ',' << upgrade_2.cols() << endl;
				//file_manager.open_ofile_ext("upgrade_2.dat") << upgrade_2 << endl;
				save_mat("upgrade_2", upgrade_2);
			}


			/*performance_log->log_event("forming x4");
			cout << "par_diff " << par_diff.rows() << ',' << par_diff.cols() << endl;
			cout << "Am " << Am.rows() << ',' << Am.cols() << endl;
			Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
			performance_log->log_event("forming x5");
			cout << "x4 " << x4.rows() << ',' << x4.cols() << endl;
			Eigen::MatrixXd x5 = Am * x4;
			performance_log->log_event("forming x6");
			cout <<"x5 " <<  x5.rows() << ',' << x5.cols() << endl;
			Eigen::MatrixXd x6 = par_diff * x5;
			performance_log->log_event("forming x7");
			cout << "V: " << V.rows() << "," << V.cols() << endl;
			cout << "ivec: " << ivec.rows() << ',' << ivec.cols() << endl;
			cout << "x6: " << x6.rows() << ',' << x6.cols() << endl;
			Eigen::MatrixXd x7 = V * ivec *V.transpose() * x6;
			cout << "x7: " << x7.rows() << x7.cols() << endl;
			performance_log->log_event("forming upgrade_2");
			Eigen::MatrixXd upgrade_2 = -1.0 * (par_diff * x7);
			cout << "upgrade_2" << endl << upgrade_2 << endl;
			return;*/

			//cout << "upgrade_2: " << upgrade_2.rows() << ',' << upgrade_2.cols() << endl;
			//cout << "oe_lam: " << pe_lam.shape().first << ',' << pe_lam.shape().second << endl;
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

	vector<map<int, int>> real_run_ids_lams;
	int best_idx = -1;
	double best_mean = 1.0e+30, best_std = 1.0e+30;
	double mean, std;
	frec << "  ---  running lambda ensembles...  ---  " << endl;
	cout << "  ---  running lambda ensembles...  ---  " << endl;
	vector<ObservationEnsemble> oe_lams = run_lambda_ensembles(pe_lams);
	frec << "  ---  evaluting lambda ensemble results  --  " << endl;
	cout << "  ---  evaluting lambda ensemble results  --  " << endl;

	frec << "  ---  last mean: " << setw(15) << last_best_mean << ", last stdev: " << setw(15) << last_best_std << endl;
	cout << "  ---  last mean: " << setw(15) << last_best_mean << ", last stdev: " << setw(15) << last_best_std << endl;

	ObservationEnsemble oe_lam_best;
	for (int i=0;i<pe_lams.size();i++)
	{	
		frec << "  --- lambda value " << lam_vals[i] << " phi summary ---  " << endl;
		cout << "  --- lambda value " << lam_vals[i] << " phi summary ---  " << endl;
		ph.update(oe_lams[i], pe_lams[i]);
		ph.report();
		mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
		std = ph.get_std(PhiHandler::phiType::COMPOSITE);
		if (mean < best_mean)
		{
			oe_lam_best = oe_lams[i];
			best_mean = mean;
			best_std = std;
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
	ph.update(oe_lam_best, pe_lams[best_idx]);
	best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
	best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
	if (best_mean < last_best_mean * 1.1)
	{
		frec << "  ---  updating parameter ensemble  ---  " << endl;
		cout << "  ---  updating parameter ensemble  ---  " << endl;
		performance_log->log_event("updating parameter ensemble");
		last_best_mean = best_mean;

		pe = pe_lams[best_idx];
		oe = oe_lam_best;
		
		if (best_std < last_best_std * 1.1)
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
		last_best_std = best_std;
	}

	else
	{
		frec << "  ---  not updating parameter ensemble  ---  " << endl;
		cout << "  ---  not updating parameter ensemble  ---  " << endl;
		ph.update(oe, pe);
		//double new_lam = last_best_mean * *max_element(lam_mults.begin(), lam_mults.end()) * 10.0;
		double new_lam = last_best_lam * 10.0;
		new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
		frec << " ---  increasing lambda from " << last_best_lam << " to " << new_lam << endl;
		cout << " ---  increasing lambda from " << last_best_lam << " to " << new_lam << endl;
		last_best_lam = new_lam;
	}
	report_and_save();
}

void IterEnsembleSmoother::report_and_save()
{
	ofstream &frec = file_manager.rec_ofstream();
	frec << endl << "  ---  IterEnsembleSmoother iteration " << iter << " report  ---  " << endl;
	frec << "   number of active realizations:  " << pe.shape().first << endl;
	frec << "   number of model runs:           " << run_mgr_ptr->get_total_runs() << endl;
	
	cout << endl << "  ---  IterEnsembleSmoother iteration " << iter << " report  ---  " << endl;
	cout << "   number of active realizations:   " << pe.shape().first << endl;
	cout << "   number of model runs:            " << run_mgr_ptr->get_total_runs() << endl;

	stringstream ss;
	ss << file_manager.get_base_filename() << "." << iter << ".obs.csv";
	oe.to_csv(ss.str());
	frec << "      current obs ensemble saved to " << ss.str() << endl;
	cout << "      current obs ensemble saved to " << ss.str() << endl;
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
	pe.to_csv(ss.str());
	frec << "      current par ensemble saved to " << ss.str() << endl;
	cout << "      current par ensemble saved to " << ss.str() << endl;
	
	ph.report();
	ph.write(iter,run_mgr_ptr->get_total_runs());
	
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
