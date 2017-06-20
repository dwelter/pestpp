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

	last_best_lam = 1.0;

	phi_stat_names.push_back("mean");
	phi_stat_names.push_back("standard deviation");
	phi_stat_names.push_back("min");
	phi_stat_names.push_back("max");
	

	performance_log->log_event("processing par csv");
	try
	{
		pe.from_csv(pest_scenario.get_pestpp_options().get_ies_par_csv());
	}
	catch (const exception &e)
	{
		stringstream ss;
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
		stringstream ss;
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

	fphi << "iteration,num_active,";
	for (auto &s : phi_stat_names)
		fphi << s << ',';
	for (auto &rname : oe.get_real_names())
		fphi << rname << ',';
	fphi << endl;
	oe_base = oe; //copy
	//reorder this for later...
	oe_base.reorder(vector<string>(), act_obs_names);

	if (true) //eventually the restart conditional
	{
		performance_log->log_event("running initial ensemble");
		EnsemblePair epair = run_ensemble(pe,oe);
		if (pe.shape().first != epair.num_active())
		{
			performance_log->log_event("resizing initial par and obs ensembles to remove failed realizations");
			pe.reorder(epair.get_pe_active_names(), vector<string>());
			oe.reorder(epair.get_oe_active_names(), vector<string>());
		}
		report_and_save();
	}
	else
	{

	}
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
	prior_pe_diff = pe.get_eigen_mean_diff();

}

void IterEnsembleSmoother::solve()
{
	iter++;
	stringstream ss;
	ss << "starting solve for iteration: " << iter;
	performance_log->log_event(ss.str());

	if (oe.shape().first < 2)
		throw_ies_error(string("less than 2 active realizations...cannot continue"));

	double scale = (1.0 / (sqrt(double(oe.shape().first - 1))));
	
	performance_log->log_event("calculate residual matrix");
	//oe_base var_names should be ordered by act_obs_names, so only reorder real_names
	//oe should only include active realizations, so only reorder var_names
	Eigen::MatrixXd scaled_residual = obscov_inv_sqrt * (oe.get_eigen(vector<string>(), act_obs_names) - oe_base.get_eigen(oe.get_real_names(), vector<string>())).transpose();

	
	performance_log->log_event("calculate scaled obs diff");
	Eigen::MatrixXd diff = oe.get_eigen_mean_diff(vector<string>(),act_obs_names).transpose();
	//cout << diff.rows() << ',' << diff.cols() << endl;
	//cout << obscov_inv_sqrt.rows() << ',' << obscov_inv_sqrt.cols() << endl;
	Eigen::MatrixXd obs_diff = scale * (obscov_inv_sqrt * diff);

	performance_log->log_event("calculate scaled par diff");
	diff = pe.get_eigen_mean_diff(vector<string>(), act_par_names).transpose();
	Eigen::MatrixXd par_diff = scale * diff;

	performance_log->log_event("SVD of obs diff");
	RedSVD::RedSVD<Eigen::MatrixXd> rsvd(obs_diff);
	
	
	vector<double> lam_mults;
	lam_mults.push_back(1.0);
	Eigen::MatrixXd ivec, upgrade_1, s = rsvd.singularValues(), V = rsvd.matrixV(), Ut = rsvd.matrixU().transpose();

	for (auto &lam_mult : lam_mults)
	{
		ss.str("");
		double cur_lam = last_best_lam * lam_mult;
		ss << "starting calcs for lambda" << cur_lam;
		performance_log->log_event("form scaled identity matrix");
		ivec = (1.0 / ((cur_lam + 1.0) * s.cwiseProduct(s)).array());
		cout << "V:" << V.rows() << ',' << V.cols() << endl;
		cout << "Ut" << Ut.rows() << ',' << Ut.cols() << endl;
		cout << "s:" << s.size() << endl;
		cout << "ivec:" << ivec.size() << endl;
		cout << "par_diff:" << par_diff.rows() << ',' << par_diff.cols() << endl;
		performance_log->log_event("calculate portion of upgrade_1 for localization");
		upgrade_1 = -1.0 * par_diff * V * s.asDiagonal() * ivec.asDiagonal() * Ut;
		//localization here...

		performance_log->log_event("apply residuals to upgrade_1");
		upgrade_1 = upgrade_1 * scaled_residual;
		cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;




	}



}


map<string,PhiComponets> IterEnsembleSmoother::get_phi_info()
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

map<string, double> IterEnsembleSmoother::get_phi_vec_stats()
{
	map<string,PhiComponets> phi_info = get_phi_info();
	map<string, double> stats;
	Eigen::VectorXd phi_vec;
	phi_vec.resize(oe.shape().first);
	PhiComponets pc;
	vector<string> obs_reals = oe.get_real_names();
	for (int ireal = 0; ireal < oe.shape().first; ireal++)
		phi_vec[ireal] = phi_info[obs_reals[ireal]].meas;
	stats["mean"] = phi_vec.mean();
	Eigen::VectorXd centered = (phi_vec - (Eigen::VectorXd::Ones(phi_vec.size()) * stats["mean"]));
	double var = (centered.cwiseProduct(centered)).sum();

	stats["standard deviation"] = sqrt(var / double(phi_vec.size()));
	stats["max"] = phi_vec.maxCoeff();
	stats["min"] = phi_vec.minCoeff();
	return stats;
}

void IterEnsembleSmoother::report_and_save()
{
	ofstream &f_rec = file_manager.rec_ofstream();
	f_rec << endl << "  ---  IterEnsembleSmoother iteration " << iter << " report  ---  " << endl;
	f_rec << "        number of active realizations: " << pe.shape().first << endl;
	f_rec << "        number of model runs: " << run_mgr_ptr->get_total_runs() << endl;
	map<string, double> phi_stats = get_phi_vec_stats();
	for (auto phi_info : phi_stats)
		f_rec << "        phi vector " << phi_info.first << ": " << setw(10) << phi_info.second << endl;
	stringstream ss;
	ss << file_manager.get_base_filename() << "." << iter << ".obs.csv";
	oe.to_csv(ss.str());
	f_rec << "      current obs ensemble saved to " << ss.str() << endl;
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
	pe.to_csv(ss.str());
	f_rec << "      current par ensemble saved to " << ss.str() << endl;
	f_rec << "      current ensemble phi info saved to " << fphi_name << endl;
	map<string, PhiComponets> phi_info = get_phi_info();

	fphi << iter << ',';
	fphi << oe.shape().first << ',';

	for (auto &s : phi_stat_names)
		fphi << phi_stats[s] << ',';
	vector<string> oe_real_names = oe.get_real_names();
	vector<string>::iterator start = oe_real_names.begin(), end = oe_real_names.end();
	for (auto &rname : oe_org_real_names)
	{
		if (find(start, end, rname) != end)
			fphi << phi_info[rname].meas;
		fphi << ',';
	}
}

EnsemblePair IterEnsembleSmoother::run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe)
{

	EnsemblePair epair(pe, oe);
	stringstream ss;
	ss << "queuing " << _pe.shape().first << " runs";
	performance_log->log_event(ss.str());
	try
	{
		epair.queue_runs(run_mgr_ptr);
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error processing obs csv: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error queuing runs"));
	}
	performance_log->log_event("making runs");
	try
	{
		epair.run(run_mgr_ptr);
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
		failed_real_indices = epair.process_runs(run_mgr_ptr);
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
	}
	return epair;
}





void IterEnsembleSmoother::finalize()
{

}