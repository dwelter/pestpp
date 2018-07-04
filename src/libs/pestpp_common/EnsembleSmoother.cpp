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
	Covariance *_parcov, double *_reg_factor, ObservationEnsemble *_weights)	
{
	pest_scenario = _pest_scenario;
	file_manager = _file_manager;
	oe_base = _oe_base;
	pe_base = _pe_base;
	weights = _weights;

	//check for inequality constraints
	//for (auto &og : pest_scenario.get_ctl_ordered_obs_group_names())
	string og;
	double weight; 
	ObservationInfo oi = pest_scenario->get_ctl_observation_info();
	for (auto &oname : pest_scenario->get_ctl_ordered_nz_obs_names())
	{
		og = oi.get_group(oname);
		weight = oi.get_weight(oname);
		if (weight == 0)
			continue;
		if ((og.compare(0, 2, "L_") == 0) || (og.compare(0, 4, "LESS")==0))
		{
			lt_obs_names.push_back(oname);
		}
		else if ((og.compare(0, 2, "G_")==0) || (og.compare(0, 7, "GREATER")==0))
		{
			gt_obs_names.push_back(oname);
		}
	}


	//build up obs group and par group idx maps for group reporting
	vector<string> nnz_obs = oe_base->get_var_names();
	ObservationInfo oinfo = pest_scenario->get_ctl_observation_info();
	vector<int> idx;
	for (auto &og : pest_scenario->get_ctl_ordered_obs_group_names())
	{
		idx.clear();
		for (int i = 0; i < nnz_obs.size(); i++)
		{
			if (oinfo.get_group(nnz_obs[i]) == og)
			{
				idx.push_back(i);
			}
		}
		if (idx.size() > 0)
			obs_group_idx_map[og] = idx;
	}

	vector<string> pars = pe_base->get_var_names();
	ParameterInfo pi = pest_scenario->get_ctl_parameter_info();
	for (auto &pg : pest_scenario->get_ctl_ordered_par_group_names())
	{
		idx.clear();
		for (int i = 0; i < pars.size(); i++)
		{
			if (pi.get_parameter_rec_ptr(pars[i])->group == pg)
				idx.push_back(i);
		}
		if (idx.size() > 0)
			par_group_idx_map[pg] = idx;
	}


	reg_factor = _reg_factor;
	//save the org reg factor and org q vector
	org_reg_factor = *_reg_factor;
	org_q_vec = get_q_vector();
	//Eigen::VectorXd parcov_inv_diag = parcov_inv.e_ptr()->diagonal();
	parcov_inv_diag = _parcov->e_ptr()->diagonal();
	for (int i = 0; i < parcov_inv_diag.size(); i++)
		parcov_inv_diag(i) = 1.0 / parcov_inv_diag(i);

	//parcov_inv = _parcov->inv();
	//parcov.inv_ip();
	oreal_names = oe_base->get_real_names();
	preal_names = pe_base->get_real_names();
	prepare_csv(file_manager->open_ofile_ext("phi.actual.csv"),oreal_names);
	prepare_csv(file_manager->open_ofile_ext("phi.meas.csv"),oreal_names);
	prepare_csv(file_manager->open_ofile_ext("phi.composite.csv"),oreal_names);
	prepare_csv(file_manager->open_ofile_ext("phi.regul.csv"),preal_names);
	prepare_group_csv(file_manager->open_ofile_ext("phi.group.csv"));

}

Eigen::MatrixXd PhiHandler::get_obs_resid(ObservationEnsemble &oe)
{
	Eigen::MatrixXd resid = oe.get_eigen(vector<string>(), oe_base->get_var_names()) -
		oe_base->get_eigen(oe.get_real_names(), vector<string>());
	
	apply_ineq_constraints(resid);
	return resid;
}


Eigen::MatrixXd PhiHandler::get_obs_resid_subset(ObservationEnsemble &oe)
{
	
	Eigen::MatrixXd resid = oe.get_eigen() - oe_base->get_eigen(oe.get_real_names(), oe.get_var_names());
	return resid;
}

Eigen::MatrixXd PhiHandler::get_par_resid(ParameterEnsemble &pe)
{
	Eigen::MatrixXd resid = pe.get_eigen(vector<string>(), pe_base->get_var_names()) -
		pe_base->get_eigen(pe.get_real_names(), vector<string>());
	return resid;
}

Eigen::MatrixXd PhiHandler::get_par_resid_subset(ParameterEnsemble &pe)
{
	Eigen::MatrixXd resid = pe.get_eigen() - pe_base->get_eigen(pe.get_real_names(),pe.get_var_names());
	return resid;
}

Eigen::MatrixXd PhiHandler::get_actual_obs_resid(ObservationEnsemble &oe)
{
	//vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	vector<string> act_obs_names = oe_base->get_var_names();
	Eigen::MatrixXd resid(oe.shape().first, act_obs_names.size());
	resid.setZero();
	Observations obs = pest_scenario->get_ctl_observations();
	Eigen::MatrixXd oe_vals = oe.get_eigen(vector<string>(), act_obs_names);
	Eigen::MatrixXd ovals = obs.get_data_eigen_vec(act_obs_names);
	ovals.transposeInPlace();
	for (int i = 0; i < resid.rows(); i++)
		resid.row(i) = oe_vals.row(i) - ovals;
	apply_ineq_constraints(resid);
	return resid;
}

Eigen::VectorXd PhiHandler::get_q_vector()
{
	ObservationInfo oinfo = pest_scenario->get_ctl_observation_info();
	Eigen::VectorXd q;
	//vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	vector<string> act_obs_names = oe_base->get_var_names();

	/*if (act_obs_names.size() == 0)
		act_obs_names = oe_base->get_var_names();*/
	q.resize(act_obs_names.size());
	double w;
	for (int i = 0; i < act_obs_names.size(); i++)
	{
		q(i) = oinfo.get_weight(act_obs_names[i]);
	}
	return q;
}

map<string, double> PhiHandler::get_obs_group_contrib(Eigen::VectorXd &phi_vec)
{
	map<string, double> group_phi_map;
	double sum;
	for (auto &og : obs_group_idx_map)
	{
		sum = 0.0;
		for (auto i : og.second)
			sum += phi_vec[i];
		group_phi_map[og.first] = sum;
	}

	return group_phi_map;
}

map<string, double> PhiHandler::get_par_group_contrib(Eigen::VectorXd &phi_vec)
{
	map<string, double> group_phi_map;
	double sum;
	for (auto &pg : par_group_idx_map)
	{
		sum = 0.0;
		for (auto i : pg.second)
			sum += phi_vec[i];
		group_phi_map[pg.first] = sum;
	}
	return group_phi_map;
}

void PhiHandler::update(ObservationEnsemble & oe, ParameterEnsemble & pe)
{
	//update the various phi component vectors
	meas.clear();
	obs_group_phi_map.clear();
	Eigen::VectorXd q = get_q_vector();
	map<string, Eigen::VectorXd> meas_map = calc_meas(oe, q);
	for (auto &pv : meas_map)
	{
		meas[pv.first] = pv.second.sum();
		
	}
	/*if (pest_scenario->get_control_info().pestmode == ControlInfo::PestMode::PARETO)
	{
		meas_map.clear();
		meas_map = calc_meas(oe, org_q_vec);
		
	}
	for (auto &pv : meas_map)
	{
		obs_group_phi_map[pv.first] = get_obs_group_contrib(pv.second);
	}*/

	regul.clear();
	map<string, Eigen::VectorXd> reg_map = calc_regul(pe);//, *reg_factor);
	//for (auto &pv : calc_regul(pe))
	string name;
	//big assumption - if oe is a diff shape, then this 
	//must be a subset, so just use the first X rows of pe
	for (int i=0;i<oe.shape().first;i++)
	{
		//name = preal_names[i];
		name = pe.get_real_names()[i];
		//cout << name << endl;
		regul[name] = reg_map[name].sum();
		par_group_phi_map[name] = get_par_group_contrib(reg_map[name]);
	}

	/*if (pest_scenario->get_control_info().pestmode == ControlInfo::PestMode::PARETO)
	{
		reg_map.clear();
		reg_map = calc_regul(pe, org_reg_factor);
	}
	
	for (int i = 0; i<oe.shape().first; i++)
	{
		name = preal_names[i];
		par_group_phi_map[name] = get_par_group_contrib(reg_map[name]);
	}*/
	
	actual.clear();
	for (auto &pv : calc_actual(oe, q))
	{
		actual[pv.first] = pv.second.sum();
		obs_group_phi_map[pv.first] = get_obs_group_contrib(pv.second);
	}
 	composite.clear();
	composite = calc_composite(meas, regul);
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
	throw runtime_error("PhiHandler::get_phi_map() didn't find a phi map...");
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
	if (var == 0.0)
		return 0.0;
	return sqrt(var/(phi_map->size()-1));
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
	/*if (*reg_factor == 0.0)
	{
		f << "    (note: reg_factor is zero; regularization phi reported but not used)" << endl;
		cout  << "    (note: reg_factor is zero; regularization phi reported but not used)" << endl;
	}*/
	f << "     current reg_factor: " << *reg_factor << endl;
	cout << "     current reg_factor: " << *reg_factor << endl;
	if (*reg_factor != 0.0)
	{
		
		f << "     note: regularization phi reported above does not " << endl;
		f << "           include the effects of reg_factor, " << endl;
		f << "           but composite phi does." << endl;
		cout << "     note: regularization phi reported above does not " << endl;
		cout << "           include the effects of reg_factor, " << endl;
		cout << "           but composite phi does." << endl;
	}
	f << endl << endl;
	f.flush();
}

void PhiHandler::write(int iter_num, int total_runs, bool write_group)
{
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.actual.csv"), phiType::ACTUAL,oreal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.meas.csv"), phiType::MEAS, oreal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.regul.csv"), phiType::REGUL, preal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.composite.csv"), phiType::COMPOSITE, oreal_names);
	if (write_group)
		write_group_csv(iter_num, total_runs, file_manager->get_ofstream("phi.group.csv"));
}

void PhiHandler::write_group(int iter_num, int total_runs, vector<double> extra)
{
	write_group_csv(iter_num, total_runs, file_manager->get_ofstream("phi.group.csv"),extra);
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



void PhiHandler::write_group_csv(int iter_num, int total_runs, ofstream &csv, vector<double> extra)
{
	//csv << "iteration,total_runs,realiation";
	string oreal, preal;
	for (int ireal = 0; ireal < oreal_names.size(); ireal++)
	{
		oreal = oreal_names[ireal];
		preal = preal_names[ireal];
		if (obs_group_phi_map.find(oreal) == obs_group_phi_map.end())
			continue;

		csv << iter_num << ',' << total_runs << ',' << oreal << ',' << preal;
		for (auto &e : extra)
			csv << ',' << e;

		for (auto &name : pest_scenario->get_ctl_ordered_obs_group_names())
			if (obs_group_phi_map[oreal].find(name) == obs_group_phi_map[oreal].end())
				csv << ',' << 0.0;
			else
				csv  << ',' << obs_group_phi_map[oreal][name];

		for (auto &name : pest_scenario->get_ctl_ordered_par_group_names())
			if (par_group_phi_map[preal].find(name) == par_group_phi_map[preal].end())
				csv << ',' << 0.0;
			else
				csv << ',' << par_group_phi_map[preal][name];
		csv << endl;;
		csv.flush();
	}
}

void PhiHandler::prepare_group_csv(ofstream &csv, vector<string> extra)
{
	csv << "iteration,total_runs,obs_realization,par_realization";
	for (auto &name : extra)
		csv << ',' << name;
	for (auto &name : pest_scenario->get_ctl_ordered_obs_group_names())
		csv << ',' << name;
	for (auto &name : pest_scenario->get_ctl_ordered_par_group_names())
		csv << ',' << name;
	csv << endl;
}

vector<int> PhiHandler::get_idxs_greater_than(double bad_phi, ObservationEnsemble &oe)
{
	map<string, double> _meas;
	Eigen::VectorXd q = get_q_vector();
	for (auto &pv : calc_meas(oe, q))
		_meas[pv.first] = pv.second.sum();
	vector<int> idxs;
	vector<string> names = oe.get_real_names();
	for (int i=0;i<names.size();i++)
		if (_meas[names[i]] > bad_phi)
			idxs.push_back(i);
	return idxs;
}

map<string, Eigen::VectorXd> PhiHandler::calc_meas(ObservationEnsemble & oe, Eigen::VectorXd &q_vec)
{
	map<string, Eigen::VectorXd> phi_map;
	Eigen::VectorXd oe_base_vec, oe_vec, diff,w_vec;
	//Eigen::VectorXd q = get_q_vector();
	vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
	vector<string>::iterator start = base_real_names.begin(), end = base_real_names.end();

	Eigen::MatrixXd w_mat;
	if (weights->shape().first > 0)
		w_mat = weights->get_eigen(vector<string>(), oe_base->get_var_names());
	double phi;
	string rname;

	Eigen::MatrixXd resid = get_obs_resid(oe);
	assert(oe_real_names.size() == resid.rows());
	for (int i = 0; i<resid.rows(); i++)
	{
		rname = oe_real_names[i];
		if (find(start, end, rname) == end)
			continue;
		diff = resid.row(i);
		
		if (weights->shape().first == 0)
			diff = diff.cwiseProduct(q_vec);
		else
		{
			w_vec = w_mat.row(i);
			diff = diff.cwiseProduct(w_vec);
		}
		phi = (diff.cwiseProduct(diff)).sum();
		phi_map[rname] = diff.cwiseProduct(diff);
	}
	return phi_map;
}

map<string, Eigen::VectorXd> PhiHandler::calc_regul(ParameterEnsemble & pe)
{	
	map<string, Eigen::VectorXd> phi_map;
	vector<string> real_names = pe.get_real_names();
	pe_base->transform_ip(ParameterEnsemble::transStatus::NUM);
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	Eigen::MatrixXd diff_mat = get_par_resid(pe);

	
	Eigen::VectorXd diff;
	for (int i = 0; i < real_names.size(); i++)
	{	
		diff = diff_mat.row(i);
		diff = diff.cwiseProduct(diff);
		diff = diff.cwiseProduct(parcov_inv_diag);
		//phi_map[real_names[i]] = _reg_fac * diff;
		phi_map[real_names[i]] = diff;
	}
	return phi_map;
}


void PhiHandler::apply_ineq_constraints(Eigen::MatrixXd &resid)
{
	vector<string> names = oe_base->get_var_names();
	vector<string> lt_names = get_lt_obs_names(), gt_names = get_gt_obs_names();
	vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	assert(names.size() == resid.cols());
	
	map<string, double> lt_vals,gt_vals;
	Observations obs = pest_scenario->get_ctl_observations();
	for (auto &n : lt_names)
		lt_vals[n] = obs.get_rec(n);
	for (auto &n : gt_names)
		gt_vals[n] = obs.get_rec(n);
	if ((lt_vals.size() == 0) && (gt_vals.size() == 0))
		return;
	map<string, int> idxs;
	//for (int i = 0; i < act_obs_names.size(); i++)
	//	idxs[act_obs_names[i]] = i;
	for (int i = 0; i < names.size(); i++)
		idxs[names[i]] = i;
	int idx;
	double val;
	Eigen::VectorXd col;
	
	for (auto iv : lt_vals)
	{
		idx = idxs[iv.first];
		col = resid.col(idx);
		val = iv.second;
		//cout << resid.col(idx) << endl;
		for (int i = 0; i < resid.rows(); i++)
			col(i) = (col(i) < 0.0) ? 0.0 : col(i);
		//cout << resid.col(idx) << endl;
		resid.col(idx) = col;
		//cout << resid.col(idx) << endl;
	}

	for (auto iv : gt_vals)
	{
		idx = idxs[iv.first];
		col = resid.col(idx);
		val = iv.second;
		for (int i = 0; i < resid.rows(); i++)
			col(i) = (col(i) > 0.0) ? 0.0 : col(i);
		resid.col(idx) = col;
	}
}


map<string, Eigen::VectorXd> PhiHandler::calc_actual(ObservationEnsemble & oe, Eigen::VectorXd &q_vec)
{
	map<string, Eigen::VectorXd> phi_map;
	Eigen::MatrixXd resid = get_actual_obs_resid(oe);
	vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
	vector<string>::iterator start = base_real_names.begin(), end = base_real_names.end();
	double phi;
	string rname;

	Eigen::MatrixXd oe_reals = oe.get_eigen(vector<string>(), oe_base->get_var_names());
	Eigen::VectorXd diff;
	for (int i = 0; i<oe.shape().first; i++)
	{
		rname = oe_real_names[i];
		if (find(start, end, rname) == end)
			continue;
		//diff = (oe_vec - obs_val_vec).cwiseProduct(q);
		diff = resid.row(i);
		diff = diff.cwiseProduct(q_vec);
		//phi = (diff.cwiseProduct(diff)).sum();
		phi_map[rname] = diff.cwiseProduct(diff);
	}
	return phi_map;
}


map<string, double> PhiHandler::calc_composite(map<string, double> &_meas, map<string, double> &_regul)
{
	map<string, double> phi_map;
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
	weights.set_pest_scenario(&pest_scenario);
	localizer.set_pest_scenario(&pest_scenario);
}

void IterEnsembleSmoother::throw_ies_error(string message)
{
	performance_log->log_event("IterEnsembleSmoother error: " + message);
	cout << endl << "   ************   " << endl << "    IterEnsembleSmoother error: " << message << endl << endl;
	file_manager.rec_ofstream() << endl << "   ************   " << endl << "    IterEnsembleSmoother error: " << message << endl << endl;
	file_manager.close_file("rec");
	performance_log->~PerformanceLog();
	throw runtime_error("IterEnsembleSmoother error: " + message);
}

bool IterEnsembleSmoother::initialize_pe(Covariance &cov)
{
	stringstream ss;
	int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
	string par_csv = pest_scenario.get_pestpp_options().get_ies_par_csv();
	bool drawn = false;
	if (par_csv.size() == 0)
	{
		message(1, "drawing parameter realizations: ", num_reals);
		pe.draw(num_reals, cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level());
		// stringstream ss;
		// ss << file_manager.get_base_filename() << ".0.par.csv";
		// message(1, "saving initial parameter ensemble to ", ss.str());
		// pe.to_csv(ss.str());
		drawn = true;
	}
	else
	{
		string par_ext = pest_utils::lower_cp(par_csv).substr(par_csv.size() - 3, par_csv.size());
		performance_log->log_event("processing par csv " + par_csv);
		if (par_ext.compare("csv") == 0)
		{
			message(1, "loading par ensemble from csv file", par_csv);
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
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
		{
			message(1, "loading par ensemble from binary file", par_csv);
			try
			{
				pe.from_binary(par_csv);
			}
			catch (const exception &e)
			{
				ss << "error processing par jcb: " << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error processing par jcb"));
			}
		}
		else
		{
			ss << "unrecognized par csv extension " << par_ext << ", looking for csv, jcb, or jco";
			throw_ies_error(ss.str());
		}

		pe.transform_ip(ParameterEnsemble::transStatus::NUM);

		if (pp_args.find("IES_NUM_REALS") != pp_args.end())
		{
			int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
			if (num_reals < pe.shape().first)
			{
				message(1,"ies_num_reals arg passed, truncated parameter ensemble to ",num_reals);
				vector<string> keep_names,real_names=pe.get_real_names();
				for (int i=0;i<num_reals;i++)
				{
					keep_names.push_back(real_names[i]);
				}
				pe.keep_rows(keep_names);
			}
		}


		if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior())
		{
			message(1, "initializing prior parameter covariance matrix from ensemble (using diagonal matrix)");
			parcov = pe.get_diagonal_cov_matrix();
			if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1)
			{
				
				if (pest_scenario.get_pestpp_options().get_ies_save_binary())
				{
					string filename = file_manager.get_base_filename() + ".prior.jcb";
					message(1, "saving emprirical parameter covariance matrix to binary file: ", filename);
					parcov.to_binary(filename);
				}
				else
				{
					string filename = file_manager.get_base_filename() + ".prior.cov";
					message(1, "saving emprirical parameter covariance matrix to ASCII file: ", filename);
					parcov.to_ascii(filename);
				}
			}
		}
		
	}
	return drawn;
	
}

void IterEnsembleSmoother::add_bases()
{
	//check that 'base' isn't already in ensemble
	vector<string> rnames = pe.get_real_names();
	bool inpar = false;
	if (find(rnames.begin(), rnames.end(), "base") != rnames.end())
	{
		message(1, "'base' realization already in parameter ensemble, ignoring '++ies_include_base'");
		inpar = true;
	}
	else
	{
		message(1, "adding 'base' parameter values to ensemble");
		Parameters pars = pest_scenario.get_ctl_parameters();
		pe.get_par_transform().active_ctl2numeric_ip(pars);
		pe.append("base", pars);
	}
	
	//check that 'base' isn't already in ensemble
	rnames = oe.get_real_names();
	if (find(rnames.begin(), rnames.end(), "base") != rnames.end())
	{		
		message(1, "'base' realization already in observation ensemble, ignoring '++ies_include_base'");
	}
	else
	{
		Observations obs = pest_scenario.get_ctl_observations();
		if (inpar)
		{
			vector<string> prnames = pe.get_real_names();

			int idx = find(prnames.begin(), prnames.end(), "base") - prnames.begin();
			cout << idx << "," << rnames.size() << endl;
			string oreal = rnames[idx];
			stringstream ss;
			ss << "warning: 'base' realization in par ensenmble but not in obs ensemble," << endl;
			ss << "         replacing obs realization '" << oreal << "' with 'base'";
			string mess = ss.str();
			message(1, mess);
			vector<string> drop;
			drop.push_back(oreal);
			oe.drop_rows(drop);
			oe.append("base", obs);
			//rnames.insert(rnames.begin() + idx, string("base"));
			rnames[idx] = "base";
			oe.reorder(rnames, vector<string>());
		}
		else
		{
			message(1, "adding 'base' observation values to ensemble");
			oe.append("base", obs);
		}
	}

	//check that 'base' isn't already in ensemble
	rnames = weights.get_real_names();
	if (rnames.size() == 0)
		return;
	if (find(rnames.begin(), rnames.end(), "base") != rnames.end())
	{
		message(1, "'base' realization already in weights ensemble, ignoring '++ies_include_base'");
	}
	else
	{
		//Observations obs = pest_scenario.get_ctl_observations();
		ObservationInfo oinfo = pest_scenario.get_ctl_observation_info();
		Eigen::VectorXd q;
		//vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
		vector<string> vnames = weights.get_var_names();
		q.resize(vnames.size());
		double w;
		for (int i = 0; i < vnames.size(); i++)
		{
			q(i) = oinfo.get_weight(vnames[i]);
		}
		
		Observations wobs(vnames, q);
		if (inpar)
		{
			vector<string> prnames = pe.get_real_names();

			int idx = find(prnames.begin(), prnames.end(), "base") - prnames.begin();
			//cout << idx << "," << rnames.size() << endl;
			string oreal = rnames[idx];
			stringstream ss;
			ss << "warning: 'base' realization in par ensenmble but not in weights ensemble," << endl;
			ss << "         replacing weights realization '" << oreal << "' with 'base'";
			string mess = ss.str();
			message(1, mess);
			vector<string> drop;
			drop.push_back(oreal);
			weights.drop_rows(drop);
			weights.append("base", wobs);
			//rnames.insert(rnames.begin() + idx, string("base"));
			rnames[idx] = "base";
			weights.reorder(rnames, vector<string>());
		}
		else
		{
			message(1, "adding 'base' weight values to weights");
			
			
			weights.append("base", wobs);
		}
	}
}

bool IterEnsembleSmoother::initialize_oe(Covariance &cov)
{
	stringstream ss;
	int num_reals = pe.shape().first;

	performance_log->log_event("load obscov");
	string obs_csv = pest_scenario.get_pestpp_options().get_ies_obs_csv();
	bool drawn = false;
	if (obs_csv.size() == 0)
	{
		message(1, "drawing observation noise realizations: ", num_reals);
		oe.draw(num_reals, cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level());
		// stringstream ss;
		// ss << file_manager.get_base_filename() << ".base.obs.csv";
		// message(1, "saving initial observation ensemble to ", ss.str());
		// oe.to_csv(ss.str());
		drawn = true;
	}
	else
	{
		string obs_ext = pest_utils::lower_cp(obs_csv).substr(obs_csv.size() - 3, obs_csv.size());
		performance_log->log_event("processing obs csv " + obs_csv);
		if (obs_ext.compare("csv") == 0)
		{
			message(1, "loading obs ensemble from csv file", obs_csv);
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
		}
		else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0))
		{
			message(1, "loading obs ensemble from binary file", obs_csv);
			try
			{
				oe.from_binary(obs_csv);
			}
			catch (const exception &e)
			{
				stringstream ss;
				ss << "error processing obs binary file: " << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error processing obs binary file"));
			}
		}
		else
		{
			ss << "unrecognized obs ensemble extension " << obs_ext << ", looing for csv, jcb, or jco";
			throw_ies_error(ss.str());
		}
		if (pp_args.find("IES_NUM_REALS") != pp_args.end())
		{
			int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
			if (num_reals < oe.shape().first)
			{
				message(1,"ies_num_reals arg passed, truncated observation ensemble to ",num_reals);
				vector<string> keep_names,real_names=oe.get_real_names();
				for (int i=0;i<num_reals;i++)
				{
					keep_names.push_back(real_names[i]);
				}
				oe.keep_rows(keep_names);
			}
		}
	}
	return drawn;
	
}

template<typename T, typename A>
void IterEnsembleSmoother::message(int level, char* _message, vector<T, A> _extras)
{ 
	string s(_message);
	message(level, s, _extras); 
}

void IterEnsembleSmoother::message(int level, char* _message)
{ 
	string s(_message);
	message(level, s); 
}

template<typename T>
void IterEnsembleSmoother::message(int level, char* _message, T extra)
{
	string s(_message);
	message(level, s, extra);

}

template<typename T, typename A>
void IterEnsembleSmoother::message(int level, string _message, vector<T, A> _extras)
{
	stringstream ss;
	if (level == 0)
		ss << endl << "  ---  ";
	else if (level == 1)
		ss << "...";
	ss << _message;
	if (_extras.size() > 0)
	{
		
		for (auto &e : _extras)
			ss << e << " , ";
		
	}
	if (level == 0)
		ss << "  ---  ";

	cout << ss.str() << endl;
	file_manager.rec_ofstream() <<ss.str() << endl;
	performance_log->log_event(ss.str());

}

void IterEnsembleSmoother::message(int level, string _message)
{
	message(level, _message, vector<string>());
}

template<typename T>
void IterEnsembleSmoother::message(int level, string _message, T extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

void IterEnsembleSmoother::sanity_checks()
{
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
	vector<string> errors;
	vector<string> warnings;
	stringstream ss;
	string par_csv = ppo->get_ies_par_csv();
	string obs_csv = ppo->get_ies_obs_csv();
	string restart = ppo->get_ies_obs_restart_csv();
	

	double acc_phi = ppo->get_ies_accept_phi_fac();
	if (acc_phi < 1.0)
		errors.push_back("ies_accept_phi_fac < 1.0, not good!");
	if (acc_phi > 10.0)
		warnings.push_back("ies_accept_phi_fac > 10.0, this is prob too big, typical values 1.05 to 1.3");

	double lam_inc = ppo->get_ies_lambda_inc_fac();
	if (lam_inc < 1.0)
		errors.push_back("ies_lambda_inc_fac < 1.0, nope! how can lambda increase if this is less than 1.0???");
	
	double lam_dec = ppo->get_ies_lambda_dec_fac();
	if (lam_dec > 1.0)
		errors.push_back("ies_lambda_dec_fac > 1.0, nope!  how can lambda decrease if this is greater than 1.0???");
	
	if ((ppo->get_ies_par_csv().size() == 0) && (ppo->get_ies_use_empirical_prior()))
	{
		warnings.push_back("no point in using an empirical prior if we are drawing the par ensemble...resetting ies_use_empirical_prior to false");
		ppo->set_ies_use_empirical_prior(false);
	}
	if ((par_csv.size() == 0) && (restart.size() > 0))
		errors.push_back("ies_par_csv is empty but ies_restart_obs_csv is not - how can this work?");
	if (ppo->get_ies_bad_phi() <= 0.0)
		errors.push_back("ies_bad_phi <= 0.0, really?");
	if ((ppo->get_ies_num_reals() < error_min_reals) && (par_csv.size() == 0))
	{
		ss.str("");
		ss << "ies_num_reals < " << error_min_reals << ", this is redic, increaing to " << warn_min_reals;
		warnings.push_back(ss.str());
		ppo->set_ies_num_reals(warn_min_reals);
	}
	if ((ppo->get_ies_num_reals() < warn_min_reals) && (par_csv.size() == 0))
	{
		ss.str("");
		ss << "ies_num_reals < " << warn_min_reals << ", this is prob too few";
		warnings.push_back(ss.str());
	}
	if (ppo->get_ies_reg_factor() < 0.0)
		errors.push_back("ies_reg_factor < 0.0 - WRONG!");
	//if (ppo->get_ies_reg_factor() > 1.0)
	//	errors.push_back("ies_reg_factor > 1.0 - nope");
	if ((par_csv.size() == 0) && (ppo->get_ies_subset_size() < 10000000) && (ppo->get_ies_num_reals() < ppo->get_ies_subset_size() * 2))
		warnings.push_back("ies_num_reals < 2*ies_subset_size: you not gaining that much using subset here");
	//if ((ppo->get_ies_subset_size() < 100000001) && (ppo->get_ies_lam_mults().size() == 1))
	//{
	//	warnings.push_back("only one lambda mult to test, no point in using a subset");
	//	//ppo->set_ies_subset_size(100000000);
	//}

	string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
	if ((how != "FIRST") && (how != "LAST") && (how != "RANDOM") && (how != "PHI_BASED"))
	{
		ss.str("");
		ss << "'subset_how' is '" << how << "' but should be 'FIRST','LAST','RANDOM','PHI_BASED'";
		errors.push_back(ss.str());
	}

	if ((ppo->get_ies_verbose_level() < 0) || (ppo->get_ies_verbose_level() > 3))
	{
		warnings.push_back("ies_verbose_level must be between 0 and 3, resetting to 3");
		ppo->set_ies_verbose_level(3);
	}
	
	/*if (!ppo->get_ies_use_prior_scaling())
	{
		warnings.push_back("not using prior scaling - this is really a dev option, you should always use prior scaling...");
	}*/

	if (warnings.size() > 0)
	{
		message(0, "sanity_check warnings");
		for (auto &w : warnings)
			message(1, w);
	}
	if (errors.size() > 0)
	{
		message(0, "sanity_check errors - uh oh");
		for (auto &e : errors)
			message(1, e);
		throw_ies_error(string("sanity_check() found some problems - please review rec file"));
	}
	//cout << endl << endl;
}

void IterEnsembleSmoother::initialize_restart_oe()
{
	stringstream ss;
	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();
	performance_log->log_event("restart IES with existing obs ensemble: " + obs_restart_csv);
	message(1, "restarting with existing obs ensemble", obs_restart_csv);
	string obs_ext = pest_utils::lower_cp(obs_restart_csv).substr(obs_restart_csv.size() - 3, obs_restart_csv.size());
	if (obs_ext.compare("csv") == 0)
	{
		message(1, "loading restart obs ensemble from csv file", obs_restart_csv);
		try
		{
			oe.from_csv(obs_restart_csv);
		}
		catch (const exception &e)
		{
			ss << "error processing restart obs csv: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error processing restart obs csv"));
		}
	}
	else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0))
	{
		message(1, "loading restart obs ensemble from binary file", obs_restart_csv);
		try
		{
			oe.from_binary(obs_restart_csv);
		}
		catch (const exception &e)
		{
			ss << "error processing restart obs binary file: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error processing restart obs binary file"));
		}
	}
	else
	{
		ss << "unrecognized restart obs ensemble extension " << obs_ext << ", looing for csv, jcb, or jco";
		throw_ies_error(ss.str());
	}

	if (pp_args.find("IES_NUM_REALS") != pp_args.end())
	{
		int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
		if (num_reals < oe.shape().first)
		{
			message(1, "ies_num_reals arg passed, truncated restart obs ensemble to ", num_reals);
			vector<string> keep_names, real_names = oe.get_real_names();
			for (int i = 0; i<num_reals; i++)
			{
				keep_names.push_back(real_names[i]);
			}
			oe.keep_rows(keep_names);
		}
	}

	//check that restart oe is in sync
	vector<string> oe_real_names = oe.get_real_names(), oe_base_real_names = oe_base.get_real_names();
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
		message(1, "shape mismatch detected with restart obs ensemble...checking for compatibility");
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


void IterEnsembleSmoother::initialize_weights()
{
	string weights_csv = pest_scenario.get_pestpp_options().get_ies_weight_csv();
	//message(1, "loading weights ensemble from ", weights_csv);

	stringstream ss;
	string obs_ext = pest_utils::lower_cp(weights_csv).substr(weights_csv.size() - 3, weights_csv.size());
	if (obs_ext.compare("csv") == 0)
	{
		message(1, "loading weights ensemble from csv file", weights_csv);
		try
		{
			weights.from_csv(weights_csv);
		}
		catch (const exception &e)
		{
			ss << "error processing weights csv: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error processing weights csv"));
		}
	}
	else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0))
	{
		message(1, "loading weights ensemble from binary file", weights_csv);
		try
		{
			weights.from_binary(weights_csv);
		}
		catch (const exception &e)
		{
			ss << "error processing weights binary file: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error processing weights binary file"));
		}
	}
	else
	{
		ss << "unrecognized weights ensemble extension " << obs_ext << ", looking for csv, jcb, or jco";
		throw_ies_error(ss.str());
	}

	if (pp_args.find("IES_NUM_REALS") != pp_args.end())
	{
		int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
		if (num_reals < oe.shape().first)
		{
			message(1, "ies_num_reals arg passed, truncated weights ensemble to ", num_reals);
			vector<string> keep_names, real_names = weights.get_real_names();
			for (int i = 0; i<num_reals; i++)
			{
				keep_names.push_back(real_names[i]);
			}
			weights.keep_rows(keep_names);
		}
	}

	//check that restart oe is in sync
	vector<string> weights_real_names = weights.get_real_names(), oe_real_names = oe.get_real_names();
	vector<string>::const_iterator start, end;
	vector<string> missing;
	start = oe_real_names.begin();
	end = oe_real_names.end();
	for (auto &rname : weights_real_names)
		if (find(start, end, rname) == end)
			missing.push_back(rname);
	if (missing.size() > 0)
	{
		ss << "the following realization names were not found in the weight csv:";
		for (auto &m : missing)
			ss << m << ",";
		throw_ies_error(ss.str());

	}


	if (weights.shape().first > oe.shape().first) //something is wrong
	{
		ss << "weight ensemble has too many rows: " << weights.shape().first << " compared to oe: " << oe.shape().first;
		throw_ies_error(ss.str());
	}
	vector<string> weights_names = weights.get_var_names();
	set<string> weights_set(weights_names.begin(), weights_names.end());
	for (auto &name : oe.get_var_names())
	{
		if (weights_set.find(name) == weights_set.end())
			missing.push_back(name);
	}

	if (missing.size() > 0)
	{
		ss << "weights ensemble is missing the following observations: ";
		for (auto m : missing)
			ss << ',' << m;
		throw_ies_error(ss.str());
	}

}

void IterEnsembleSmoother::initialize_parcov()
{
	stringstream ss;
	performance_log->log_event("initializing parcov");

	if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior())
		return;
	string how = parcov.try_from(pest_scenario, file_manager);
	message(1, "parcov loaded ", how);
	//if (parcov.e_ptr()->rows() > 0)
	parcov = parcov.get(act_par_names);

}


void IterEnsembleSmoother::initialize_obscov()
{
	message(1, "initializing observation noise covariance matrix");
	string obscov_filename = pest_scenario.get_pestpp_options().get_obscov_filename();
	
	string how = obscov.try_from(pest_scenario, file_manager, false);
	message(1, "obscov loaded ", how);
	obscov = obscov.get(act_obs_names);



}


void IterEnsembleSmoother::initialize()
{
	message(0, "initializing");
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();
	stringstream ss;

	//set some defaults
	PestppOptions *ppo = pest_scenario.get_pestpp_options_ptr();
	
	if (pp_args.find("IES_LAMBDA_MULTS") == pp_args.end())
		ppo->set_ies_lam_mults(vector<double>{0.1, 0.5, 1.0, 2.0, 5.0});
	if (pp_args.find("IES_SUBSET_SIZE") == pp_args.end())
		ppo->set_ies_subset_size(4);
	if (pp_args.find("LAMBDA_SCALE_FAC") == pp_args.end())
		ppo->set_lambda_scale_vec(vector<double>{0.5, 0.75, 0.95, 1.0, 1.1});


	verbose_level = pest_scenario.get_pestpp_options_ptr()->get_ies_verbose_level();
	if (pest_scenario.get_n_adj_par() >= 1e6)
	{
		message(0, "You are a god among mere mortals!");
	}

	use_localizer = localizer.initialize(performance_log);
	if (use_localizer)
	{
		ss.str("");
		ss << "using localized solution with " << localizer.get_localizer_map().size() << " sequential upgrade steps";
		message(1, ss.str());
		ss.str("");
	}
	iter = 0;
	//ofstream &frec = file_manager.rec_ofstream();
	last_best_mean = 1.0E+30;
	last_best_std = 1.0e+30;
	lambda_max = 1.0E+30;
	lambda_min = 1.0E-30;
	warn_min_reals = 10;
	error_min_reals = 2;
	consec_bad_lambda_cycles = 0;

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::REGUL)
	{
		message(1, "WARNING: 'pestmode' == 'regularization', in pestpp-ies, this is controlled with the ++ies_reg_factor argument, resetting to 'estimation'");
		//throw_ies_error("'pestmode' == 'regularization', please reset to 'estimation'");
	}
	else if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::UNKNOWN)
	{
		message(1,"WARNING: unrecognized 'pestmode', using 'estimation'");
	}
	else if ((pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::PARETO))
	{
		message(1, "using pestpp-ies 'pareto' mode");
		string pobs_group = pest_scenario.get_pareto_info().obsgroup;
		message(1, "pareto obs group: ", pobs_group);

		if (pobs_group.substr(0, 5) == "REGUL")
		{
			message(1, "'regul' detected in pareto obs group name");
			//if (pest_scenario.get_pestpp_options().get_ies_reg_factor() == 0.0)
			//	throw_ies_error("pareto model problem: pareto obs group is 'regul'-ish but ies_reg_fac is zero");
		
		}
		else
		{
			
			ObservationInfo *oi = pest_scenario.get_observation_info_ptr();
			Observations obs = pest_scenario.get_ctl_observations();
			vector<string> zero;
			for (auto &oname : pest_scenario.get_ctl_ordered_obs_names())
				if (oi->get_group(oname) == pobs_group)
				{
					if (oi->get_weight(oname) == 0.0)
						zero.push_back(oname);
					pareto_obs[oname] = obs[oname];
					pareto_weights[oname] = oi->get_weight(oname);
				}

			if (pareto_obs.size() == 0)
				throw_ies_error("no observations found for pareto obs group");
			
			
			if (zero.size() > 0)
			{
				message(0, "error: the following pareto obs have zero weight:", zero);
				throw_ies_error("atleast one pareto obs has zero weight");
			}


		}
		//throw_ies_error("pareto mode not finished");
	}

	if (pest_scenario.get_ctl_ordered_pi_names().size() > 0)
	{
		message(1, "WARNING: prior information equations not supported in pestpp-ies, ignoring...");
	}

	lam_mults = pest_scenario.get_pestpp_options().get_ies_lam_mults();
	if (lam_mults.size() == 0)
		lam_mults.push_back(1.0);
	message(1, "using lambda multipliers: ", lam_mults);
	vector<double> scale_facs = pest_scenario.get_pestpp_options().get_lambda_scale_vec();
	message(1, "usnig lambda scaling factors: ", scale_facs);
	double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
	message(1, "acceptable phi factor: ", acc_fac);
	double inc_fac = pest_scenario.get_pestpp_options().get_ies_lambda_inc_fac();
	message(1, "lambda increase factor: ", inc_fac);
	double dec_fac = pest_scenario.get_pestpp_options().get_ies_lambda_dec_fac();
	message(1, "lambda decrease factor: ", dec_fac);
	message(1, "max run fail: ", ppo->get_max_run_fail());
	string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
	message(1, "subset how: ", how);
	sanity_checks();
	
	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	initialize_parcov();
	initialize_obscov();

	subset_size = pest_scenario.get_pestpp_options().get_ies_subset_size();
	reg_factor = pest_scenario.get_pestpp_options().get_ies_reg_factor();
	message(1,"using reg_factor: ", reg_factor);
	double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
	if (bad_phi < 1.0e+30)
		message(1, "using bad_phi: ", bad_phi);

	

	int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
	
	
	
	bool pe_drawn = initialize_pe(parcov);

	if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
	{
		message(1, "forming inverse sqrt of prior parameter covariance matrix");
		
		if (parcov.isdiagonal())
			parcov_inv_sqrt = parcov.inv(echo).get_matrix().diagonal().cwiseSqrt().asDiagonal();
		else
		{
			message(1, "first extracting diagonal from prior parameter covariance matrix");
			Covariance parcov_diag;
			parcov_diag.from_diagonal(parcov);
			parcov_inv_sqrt = parcov_diag.inv(echo).get_matrix().diagonal().cwiseSqrt().asDiagonal();
		}
	}
	else {
		message(1, "not using prior parameter covariance matrix scaling");
	}

	bool oe_drawn = initialize_oe(obscov);

	try
	{
		pe.check_for_dups();
	}
	catch (const exception &e)
	{
		string message = e.what();
		throw_ies_error("error in parameter ensemble: " + message);
	}

	try
	{
		oe.check_for_dups();
	}
	catch (const exception &e)
	{
		string message = e.what();
		throw_ies_error("error in observation ensemble: " + message);
	}

	if (pe.shape().first != oe.shape().first)
	{
		ss.str("");
		ss << "parameter ensemble rows (" << pe.shape().first << ") not equal to observation ensemble rows (" << oe.shape().first << ")";
		throw_ies_error(ss.str());
	}

	if (ppo->get_ies_weight_csv().size() > 0)
	{
		initialize_weights();
	}

	//if pareto mode, reset the stochastic obs vals for the pareto obs to the value in the control file
	if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::PARETO)
	{
		string oname;
		vector<string> oe_names = oe.get_var_names();
		double val;
		Eigen::VectorXd col;
		Eigen::MatrixXd oe_reals = *oe.get_eigen_ptr();
		for (int i = 0; i < oe.shape().second; i++)
		{
			if (pareto_obs.find(oe_names[i]) != pareto_obs.end())
			{
				val = pareto_obs[oe_names[i]];
				col = Eigen::VectorXd::Zero(oe.shape().first);
				col.setConstant(val);
				oe_reals.col(i) = col;
			}
		}
		oe.set_eigen(oe_reals);
	}
	
	//need this here for Am calcs...
	message(1, "transforming parameter ensembles to numeric");
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);

	if (pest_scenario.get_pestpp_options().get_ies_include_base())
		add_bases();

	ss.str("");
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << ".0.par.jcb";
		pe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".0.par.csv";
		pe.to_csv(ss.str());
	}
	message(1, "saved initial parameter ensemble to ", ss.str());
	
	ss.str("");
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << ".base.obs.jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".base.obs.csv";
		oe.to_csv(ss.str());
	}
	message(1, "saved base observation ensemble (obsval+noise) to ", ss.str());
	
	if (pest_scenario.get_control_info().noptmax == 0)
	{
		message(0, "'noptmax'=0, running mean parameter ensemble values and quitting");
		message(1, "calculating mean parameter values");
		Parameters pars;
		vector<double> mv = pe.get_mean_stl_vector();
		pars.update(pe.get_var_names(), pe.get_mean_stl_vector());
		ParamTransformSeq pts = pe.get_par_transform();

		ParameterEnsemble _pe(&pest_scenario);
		_pe.reserve(vector<string>(), pe.get_var_names());
		_pe.set_trans_status(pe.get_trans_status());
		_pe.append("mean", pars);
		string par_csv = file_manager.get_base_filename() + ".mean.par.csv";
		message(1, "saving mean parameter values to ", par_csv);
		_pe.to_csv(par_csv);
		pe_base = _pe;
		pe_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario);
		_oe.reserve(vector<string>(), oe.get_var_names());
		_oe.append("mean", pest_scenario.get_ctl_observations());
		oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);
		//initialize the phi handler
		ph = PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov, &reg_factor, &weights);
		if (ph.get_lt_obs_names().size() > 0)
		{
			message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
		}
		if (ph.get_gt_obs_names().size())
		{
			message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
		}
		message(1, "running mean parameter values");

		vector<int> failed_idxs = run_ensemble(_pe, _oe);
		if (failed_idxs.size() != 0)
		{
			message(0, "mean parmeter value run failed...bummer");
			return;
		}
		string obs_csv = file_manager.get_base_filename() + ".mean.obs.csv";
		message(1, "saving results from mean parameter value run to ", obs_csv);
		_oe.to_csv(obs_csv);

		ph.update(_oe, _pe);
		message(0, "mean parameter phi report:");
		ph.report();

		return;
	}

	if (subset_size > pe.shape().first)
	{
		use_subset = false;
	}
	else
	{
		message(1, "using subset in lambda testing, only first ", subset_size);
		use_subset = true;
	}
	
	oe_org_real_names = oe.get_real_names();
	pe_org_real_names = pe.get_real_names();

	oe_base = oe; //copy
	//reorder this for later...
	oe_base.reorder(vector<string>(), act_obs_names);

	pe_base = pe; //copy
	//reorder this for later
	pe_base.reorder(vector<string>(), act_par_names);
	
	message(1, "forming inverse sqrt obscov");
	//obscov.inv_ip(echo);
	/*obscov_inv_sqrt = obscov.inv().get_matrix().diagonal().cwiseSqrt().asDiagonal();
	if (verbose_level > 2)
	{
		Eigen::MatrixXd oc = obscov_inv_sqrt.toDenseMatrix();
		save_mat("obscov_inv_sqrt.dat", oc);
	}*/
	

	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();	
	//no restart
	if (obs_restart_csv.size() == 0)
	{
		performance_log->log_event("running initial ensemble");
		message(1, "running initial ensemble of size", oe.shape().first);
		vector<int> failed = run_ensemble(pe, oe);
		if (pe.shape().first == 0)
			throw_ies_error("all realizations failed during initial evaluation");
		
		ss.str("");
		if (pest_scenario.get_pestpp_options().get_ies_save_binary())
		{
			ss << file_manager.get_base_filename() << ".0.obs.jcb";
			oe.to_binary(ss.str());
		}
		else
		{
			ss << file_manager.get_base_filename() << ".0.obs.csv";
			oe.to_csv(ss.str());
		}
		message(1, "saving results of initial ensemble run to", ss.str());

		//string obs_csv = file_manager.get_base_filename() + ".0.obs.csv";
		//message(1, "saving results of initial ensemble run to", obs_csv);
		//oe.to_csv(obs_csv);
		pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	}
	else
	{
		initialize_restart_oe();
	}



	performance_log->log_event("calc initial phi");
	//initialize the phi handler
	ph = PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov, &reg_factor, &weights);

	if (ph.get_lt_obs_names().size() > 0)
	{
		message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
	}
	if (ph.get_gt_obs_names().size())
	{
		message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
	}

	ph.update(oe, pe);
	message(0, "pre-drop initial phi summary");
	ph.report();
	drop_bad_phi(pe, oe);
	if (oe.shape().first == 0)
	{
		throw_ies_error(string("all realizations dropped as 'bad'"));
	}
	if (oe.shape().first <= error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need at least ", error_min_reals);
		throw_ies_error(string("too few active realizations, cannot continue"));
	}
	if (oe.shape().first <= warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(0, s);
	}
	performance_log->log_event("calc initial phi");
	ph.update(oe, pe);
	message(0, "initial phi summary");
	ph.report();
	ph.write(0, run_mgr_ptr->get_total_runs());
	
	if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
	{
		message(1, "using full (MAP) update solution");
		/*performance_log->log_event("calculating 'Am' matrix for full solution");
		message(1, "forming Am matrix");
		double scale = (1.0 / (sqrt(double(pe.shape().first - 1))));
		Eigen::MatrixXd par_diff = scale * pe.get_eigen_mean_diff();
		par_diff.transposeInPlace();
		if (verbose_level > 1)
		{
			cout << "prior_par_diff: " << par_diff.rows() << ',' << par_diff.cols() << endl;
			if (verbose_level > 2)
				save_mat("prior_par_diff.dat", par_diff);
		}

		Eigen::MatrixXd ivec, upgrade_1, s, V, U, st;
		SVD_REDSVD rsvd;
		rsvd.set_performance_log(performance_log);

		rsvd.solve_ip(par_diff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
		par_diff.resize(0, 0);
		Eigen::MatrixXd temp = s.asDiagonal();
		Eigen::MatrixXd temp2 = temp.inverse();
		Am = U * temp;
		if (verbose_level > 1)
		{
			cout << "Am:" << Am.rows() << ',' << Am.cols() << endl;
			if (verbose_level > 2)
			{
				save_mat("am.dat", Am);
				save_mat("am_u.dat", U);
				save_mat("am_v.dat", V);
				save_mat("am_s_inv.dat", temp2);
			}
		}*/
	}

	last_best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
	last_best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
	last_best_lam = pest_scenario.get_pestpp_options().get_ies_init_lam();
	if (last_best_lam <= 0.0)
	{
		double x = last_best_mean / (2.0 * double(oe.shape().second));
		last_best_lam = pow(10.0,(floor(log10(x))));
	}
	
	message(1, "current lambda:", last_best_lam);
	message(0, "initialization complete");
}

Eigen::MatrixXd IterEnsembleSmoother::get_Am(const vector<string> &real_names, const vector<string> &par_names)
{
	
	double scale = (1.0 / (sqrt(double(real_names.size() - 1))));
	Eigen::MatrixXd par_diff = scale * pe_base.get_eigen_mean_diff(real_names,par_names);
	par_diff.transposeInPlace();
	if (verbose_level > 1)
	{
		cout << "prior_par_diff: " << par_diff.rows() << ',' << par_diff.cols() << endl;
		if (verbose_level > 2)
			save_mat("prior_par_diff.dat", par_diff);
	}

	Eigen::MatrixXd ivec, upgrade_1, s, V, U, st;
	SVD_REDSVD rsvd;
	rsvd.set_performance_log(performance_log);

	rsvd.solve_ip(par_diff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
	par_diff.resize(0, 0);
	Eigen::MatrixXd temp = s.asDiagonal();
	Eigen::MatrixXd temp2 = temp.inverse();
	Eigen::MatrixXd Am = U * temp;
	return Am;
}

void IterEnsembleSmoother::drop_bad_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe)
{
	//don't use this assert because _pe maybe full size, but _oe might be subset size
	//assert(_pe.shape().first == _oe.shape().first);
	vector<int> idxs = ph.get_idxs_greater_than(pest_scenario.get_pestpp_options().get_ies_bad_phi(), _oe);
	
	//for testing
	//idxs.push_back(0);

	if (idxs.size() > 0)
	{

		message(0, "droppping realizations as bad: ", idxs.size());

		vector<string> par_real_names = _pe.get_real_names(), obs_real_names = _oe.get_real_names();
		stringstream ss;
		string pname;
		int pidx;
		vector<string> full_onames;
		// if a subset drop, then use the full oe index, otherwise, just use _oe index
		if (_oe.shape().first != _pe.shape().first)
		{
			full_onames = oe.get_real_names();
		}
		else
		{
			full_onames = _oe.get_real_names();
		}
		vector<string> pdrop, odrop;
		for (auto i : idxs)
		{
			//find the index of this i in the full set of obs names - for the case of randomized subsets
			pidx = find(full_onames.begin(), full_onames.end(), obs_real_names[i]) - full_onames.begin();
			pname = par_real_names[pidx];
			ss << pname << " : " << obs_real_names[i] << " , ";
			pdrop.push_back(pname);
			odrop.push_back(obs_real_names[i]);
		}

		string s = "dropping par:obs realizations: "+ ss.str();
		message(1, s);
		try
		{
			_pe.drop_rows(pdrop);
			_oe.drop_rows(odrop);
		}
		catch (const exception &e)
		{
			stringstream ss;
			ss << "drop_bad_phi() error : " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("drop_bad_phi() error"));
		}
	}
}

void IterEnsembleSmoother::save_mat(string prefix, Eigen::MatrixXd &mat)
{
	stringstream ss;
	ss << iter << '.' << prefix;
	try
	{	
		ofstream &f = file_manager.open_ofile_ext(ss.str());
		f << mat << endl;
		f.close();
		file_manager.close_file(ss.str());
	}
	catch (...)
	{
		message(1, "error saving matrix", ss.str());
	}
}

void IterEnsembleSmoother::adjust_pareto_weight(string &obsgroup, double wfac)

{
	if (obsgroup.substr(0, 5) == "REGUL")

	{
		message(1, "resetting reg fac for pareto to ", wfac);
		reg_factor = wfac;
	}
	else
	{
		message(1, "applying weight factor for pareto obs group of ", wfac);

		//pest_scenario.get_observation_info_ptr()->scale_group_weights(obsgroup, wfac);
		ObservationInfo *oi = pest_scenario.get_observation_info_ptr();
		double val;
		for (auto &pw : pareto_weights)
		{
			val = wfac * pw.second;
			oi->set_weight(pw.first, val);
		}
			

		Covariance obscov;
		obscov.from_observation_weights(pest_scenario);
		obscov = obscov.get(act_obs_names);
		cout << obscov_inv_sqrt.diagonal() << endl;
		obscov_inv_sqrt = obscov.inv().get_matrix().diagonal().cwiseSqrt().asDiagonal();
	    cout << obscov_inv_sqrt.diagonal() << endl;
		cout << endl;
	}
}

void IterEnsembleSmoother::pareto_iterate_2_solution()
{
	//todo:
	//get initial obsgroup weight and use wf mults intead of vals
	ParetoInfo pi = pest_scenario.get_pareto_info();
	stringstream ss;
	//double init_lam = last_best_lam, init_mean = last_best_mean, init_std = last_best_std;
	double init_lam = last_best_lam, init_mean = 1.0e+30, init_std = 1.0e+30;

	message(0, "starting pareto analysis");
	message(0, "initial pareto wfac", pi.wf_start);
	message(1, "starting initial pareto iterations", pi.niter_start);
	adjust_pareto_weight(pi.obsgroup, pi.wf_start);
	last_best_lam = init_lam, last_best_mean = init_mean, last_best_std = init_std;
	for (int i = 0; i < pi.niter_start; i++)
	{
		iter++;
		message(1, "starting solve for iteration:", iter);
		ss << "starting solve for iteration: " << iter;
		performance_log->log_event(ss.str());
		solve_new();
		report_and_save();
		ph.update(oe, pe);
		last_best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
		last_best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
		ph.report();
		ph.write(iter, run_mgr_ptr->get_total_runs(),false);

	}
	ph.write_group(iter, run_mgr_ptr->get_total_runs(),vector<double>());
	double wfac = pi.wf_start + pi.wf_inc;
	vector<double> wfacs;
	wfacs.push_back(pi.wf_start + pi.wf_inc);
	while (true)
	{
		if (pi.wf_inc < 0.0)
		{
			if (wfac < pi.wf_fin)
				break;
			//wfac = wfac - pi.wf_inc;
		}
		else
		{
			if (wfac > pi.wf_fin)
				break;
			
		}
		wfac = wfac + pi.wf_inc;
		wfacs.push_back(wfac);
	}
	//while (wfac < pi.wf_fin)
	//pe = pe_base;
	//oe = oe_base;
	for (auto &wfac : wfacs)
	{
		last_best_lam = init_lam, last_best_mean = init_mean, last_best_std = init_std;
		message(0, "using pareto wfac", wfac);
		message(0, "starting pareto iterations", pi.niter_gen);
		adjust_pareto_weight(pi.obsgroup, wfac);
		for (int i = 0; i < pi.niter_gen; i++)
		{
			iter++;
			message(0, "starting solve for iteration:", iter);
			ss << "starting solve for iteration: " << iter;
			performance_log->log_event(ss.str());
			solve_new();
			report_and_save();
			ph.update(oe, pe);
			last_best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
			last_best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
			ph.report();
			ph.write(iter, run_mgr_ptr->get_total_runs(),false);
		}
		ph.write_group(iter, run_mgr_ptr->get_total_runs(), vector<double>());
		//pe = pe_base;
		//oe = oe_base;
	}
	message(1, "final pareto wfac", pi.niter_fin);
	message(0, "starting final pareto iterations", pi.niter_fin);
	adjust_pareto_weight(pi.obsgroup, pi.wf_fin);
	last_best_lam = init_lam, last_best_mean = init_mean, last_best_std = init_std;
	//pe = pe_base;
	//oe = oe_base;
	for (int i = 0; i < pi.niter_fin; i++)
	{
		iter++;
		message(0, "starting solve for iteration:", iter);
		ss << "starting solve for iteration: " << iter;
		performance_log->log_event(ss.str());
		solve_new();
		report_and_save();
		ph.update(oe,pe);
		last_best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
		last_best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
		ph.report();
		ph.write(iter, run_mgr_ptr->get_total_runs(),false);
	}
	ph.write_group(iter, run_mgr_ptr->get_total_runs(), vector<double>());

	
}

void IterEnsembleSmoother::iterate_2_solution()
{
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();
	if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::PARETO)
		pareto_iterate_2_solution();
	else
	{
		bool accept;
		for (int i = 0; i < pest_scenario.get_control_info().noptmax; i++)
		{
			iter++;
			message(0, "starting solve for iteration:", iter);
			ss << "starting solve for iteration: " << iter;
			performance_log->log_event(ss.str());
			accept = solve_new();
			report_and_save();
			ph.update(oe,pe);
			last_best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
			last_best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
			ph.report();
			ph.write(iter, run_mgr_ptr->get_total_runs());
			if (accept)
				consec_bad_lambda_cycles = 0;
			else
				consec_bad_lambda_cycles++;
			
			if (should_terminate())
				break;
		}
	}
}

bool IterEnsembleSmoother::should_terminate()
{
	//todo: use ies accept fac here?
	double phiredstp = pest_scenario.get_control_info().phiredstp;
	int nphistp = pest_scenario.get_control_info().nphistp;
	int nphinored = pest_scenario.get_control_info().nphinored;
	bool phiredstp_sat = false, nphinored_sat = false, consec_sat = false;
	double phi, ratio;
	int count = 0;
	int nphired = 0;
	//best_mean_phis = vector<double>{ 1.0,0.8,0.81,0.755,1.1,0.75,0.75,1.2 };

	

	/*if ((!consec_sat )&& (best_mean_phis.size() == 0))
		return false;*/
	message(0, "phi-based termination criteria check");
	message(1, "phiredstp: ", phiredstp);
	message(1, "nphistp: ", nphistp);
	message(1, "nphinored (also used for consecutive bad lamdba cycles): ", nphinored);
	if (best_mean_phis.size() > 0)
	{
		vector<double>::iterator idx = min_element(best_mean_phis.begin(), best_mean_phis.end());
		nphired = (best_mean_phis.end() - idx) - 1;
		best_phi_yet = best_mean_phis[idx - best_mean_phis.begin()];// *pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
		message(1, "best mean phi sequence: ", best_mean_phis);
		message(1, "best phi yet: ", best_phi_yet);
	}
	message(1, "number of consecutive bad lambda testing cycles: ", consec_bad_lambda_cycles);
	if (consec_bad_lambda_cycles >= nphinored)
	{
		message(1, "number of consecutive bad lambda testing cycles > nphinored");
		consec_sat = true;
	}
			
	for (auto &phi : best_mean_phis)
	{
		ratio = (phi - best_phi_yet) / phi;
		if (ratio <= phiredstp)
			count++;
	}
	message(1, "number of iterations satisfying phiredstp criteria: ", count);
	if (count >= nphistp)
	{
		message(1, "number iterations satisfying phiredstp criteria > nphistp");
		phiredstp_sat = true;
	}
	
	message(1, "number of iterations since best yet mean phi: ", nphired);
	if (nphired >= nphinored)
	{
		message(1, "number of iterations since best yet mean phi > nphinored");
		nphinored_sat = true;
	}
	
	if ((nphinored_sat) || (phiredstp_sat) || (consec_sat))
	{
		message(1, "phi-based termination criteria satisfied, all done");
		return true;
	}
	return false;
}

ParameterEnsemble IterEnsembleSmoother::calc_upgrade(vector<string> &obs_names, vector<string> &par_names,double cur_lam, int num_reals)
{

	stringstream ss;

	//make sure par_names are in pe var_names and obs_names in oe var_names

	//use num_reals to select out a block from eigen, rather than getting and reordering ensembles

	ObservationEnsemble oe_upgrade = oe;
	oe_upgrade.reorder(vector<string>(), obs_names);
	ParameterEnsemble pe_upgrade = pe;
	pe_upgrade.reorder(vector<string>(), par_names);

	Eigen::DiagonalMatrix<double, Eigen::Dynamic> parcov_inv;// = parcov.get(par_names).inv().e_ptr()->toDense().cwiseSqrt().asDiagonal();
	if (parcov.isdiagonal())
		parcov_inv = parcov.get(par_names).inv().get_matrix().diagonal().cwiseSqrt().asDiagonal();
	else
	{
		message(1, "first extracting diagonal from prior parameter covariance matrix");
		Covariance parcov_diag;
		Covariance parcov_local = parcov.get(par_names);
		parcov_diag.from_diagonal(parcov_local);
		parcov_inv = parcov_diag.inv().get_matrix().diagonal().cwiseSqrt().asDiagonal();
	}


	double scale = (1.0 / (sqrt(double(oe_upgrade.shape().first - 1))));
	Eigen::VectorXd weight_vec(obs_names.size());
	for (int i = 0; i < obs_names.size(); i++)
	{
		weight_vec[i] = pest_scenario.get_observation_info_ptr()->get_weight(obs_names[i]);
	}
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> weights = weight_vec.asDiagonal();

	performance_log->log_event("calculate residual matrix");
	//oe_base var_names should be ordered by act_obs_names, so only reorder real_names
	//oe should only include active realizations, so only reorder var_names
	message(1, "calculating residual matrix");
	//Eigen::MatrixXd scaled_residual = obscov_inv_sqrt * ph.get_obs_resid(oe).transpose();
	Eigen::MatrixXd scaled_residual = weights * ph.get_obs_resid_subset(oe_upgrade).transpose();

	if (verbose_level > 1)
	{
		cout << "scaled_residual: " << scaled_residual.rows() << ',' << scaled_residual.cols() << endl;
		if (verbose_level > 2)
		{
			save_mat("scaled_residual.dat", scaled_residual);
			Eigen::MatrixXd residual = ph.get_obs_resid_subset(oe_upgrade).transpose();
			save_mat("residual.dat", residual);
		}
	}

	performance_log->log_event("calculate scaled obs diff");
	message(1, "calculating obs diff matrix");
	Eigen::MatrixXd diff = oe_upgrade.get_eigen_mean_diff().transpose();
	Eigen::MatrixXd obs_diff = scale * (weights * diff);
	if (verbose_level > 1)
	{
		cout << "obs_diff: " << obs_diff.rows() << ',' << obs_diff.cols() << endl;
		if (verbose_level > 2)
			save_mat("obs_diff.dat", obs_diff);
	}

	performance_log->log_event("calculate scaled par diff");
	message(1, "calculating par diff matrix");
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	diff = pe_upgrade.get_eigen_mean_diff().transpose();
	Eigen::MatrixXd par_diff;
	if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
	{
		//throw runtime_error("parcov scale not implemented for localization");
		cout << "...applying prior par cov scaling to par diff matrix" << endl;
		par_diff = scale * parcov_inv * diff;
	}
	else
		par_diff = scale * diff;

	if (verbose_level > 1)
	{
		cout << "par_diff:" << par_diff.rows() << ',' << par_diff.cols() << endl;
		if (verbose_level > 2)
		{
			save_mat("scaled_par_diff.dat", par_diff);
			save_mat("par_diff.dat", diff);
		}
	}
	performance_log->log_event("SVD of obs diff");
	message(1, "calculating SVD of obs diff matrix");
	Eigen::MatrixXd ivec, upgrade_1, s, V, Ut;

	
	SVD_REDSVD rsvd;
	rsvd.set_performance_log(performance_log);
	rsvd.solve_ip(obs_diff, s, Ut, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);

	//SVD_EIGEN esvd;
	//esvd.set_performance_log(performance_log);
	//esvd.solve_ip(obs_diff, s, Ut, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);

	Ut.transposeInPlace();
	obs_diff.resize(0, 0);

	Eigen::MatrixXd s2 = s.cwiseProduct(s);
	if (verbose_level > 1)
	{
		cout << "s2: " << s2.rows() << ',' << s2.cols() << endl;
		cout << "Ut: " << Ut.rows() << ',' << Ut.cols() << endl;
		cout << "V:" << V.rows() << ',' << V.cols() << endl;
		if (verbose_level > 2)
		{
			save_mat("ut.dat", Ut);
			save_mat("s2.dat", s2);
		}
	}

	//vector<ParameterEnsemble> pe_lams;
	//vector<double> lam_vals, scale_vals;
	
	/*ss.str("");
	ss << "starting calcs for lambda" << cur_lam;
	message(0, "starting lambda calcs for lambda", cur_lam);\
	performance_log->log_event(ss.str());*/

	performance_log->log_event("form scaled identity matrix");
	message(1, "calculating scaled identity matrix");
	ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();
	if (verbose_level > 1)
	{
		cout << "ivec:" << ivec.rows() << ',' << ivec.cols() << endl;
		if (verbose_level > 2)
			save_mat("ivec.dat", ivec);
	}


	message(1, "forming X1");
	Eigen::MatrixXd X1 = Ut * scaled_residual;
	if (verbose_level > 1)
	{
		cout << "X1: " << X1.rows() << ',' << X1.cols() << endl;
		if (verbose_level > 2)
			save_mat("X1.dat", X1);
	}
	//scaled_residual.resize(0, 0);
	//Ut.resize(0, 0);

	message(1, "forming X2");
	Eigen::MatrixXd X2 = ivec * X1;
	if (verbose_level > 1)
	{
		cout << "X2: " << X2.rows() << ',' << X2.cols() << endl;
		if (verbose_level > 2)
			save_mat("X2.dat", X2);
	}
	X1.resize(0, 0);

	message(1, "forming X3");
	Eigen::MatrixXd X3 = V * s.asDiagonal() * X2;
	if (verbose_level > 1)
	{
		cout << "X3: " << X3.rows() << ',' << X3.cols() << endl;
		if (verbose_level > 2)
			save_mat("X3.dat", X3);
	}
	X2.resize(0, 0);

	message(1, "forming upgrade_1");
	if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
	{
		upgrade_1 = -1.0 * parcov_inv * par_diff * X3;
	}
	else
	{
		upgrade_1 = -1.0 * par_diff * X3;
	}
	upgrade_1.transposeInPlace();
	if (verbose_level > 1)
	{
		cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;
		if (verbose_level > 2)
			save_mat("upgrade_1.dat", upgrade_1);
	}
	X3.resize(0, 0);

	/*ParameterEnsemble pe_lam = pe;
	pe_lam.set_eigen(*pe_lam.get_eigen_ptr() + upgrade_1);
	upgrade_1.resize(0, 0);*/

	Eigen::MatrixXd upgrade_2;
	if ((!pest_scenario.get_pestpp_options().get_ies_use_approx()))// && (iter > 1))
	{
		performance_log->log_event("calculating parameter correction (full solution)");
		message(1, "calculating parameter correction (full solution, MAP)");
		performance_log->log_event("forming scaled par resid");
		Eigen::MatrixXd scaled_par_resid;
		if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
		{
			//throw runtime_error("parcov scaling not implemented for localization");
			scaled_par_resid = parcov_inv * ph.get_par_resid(pe_upgrade).transpose();
		}
		else
		{
			scaled_par_resid = ph.get_par_resid_subset(pe_upgrade).transpose();
		}

		//scaled_par_resid.transposeInPlace();

		performance_log->log_event("forming x4");
		message(1, "forming X4");
		if (verbose_level > 1)
		{
			cout << "scaled_par_resid: " << scaled_par_resid.rows() << ',' << scaled_par_resid.cols() << endl;
			if (verbose_level > 2)
				save_mat("scaled_par_resid.dat", scaled_par_resid);
		}
		Eigen::MatrixXd Am = get_Am(pe.get_real_names(), par_names);
		Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
		if (verbose_level > 1)
		{
			cout << "x4: " << x4.rows() << ',' << x4.cols() << endl;
			if (verbose_level > 2)
				save_mat("x4.dat", x4);
		}

		performance_log->log_event("forming x5");
		message(1, "forming X5");
		Eigen::MatrixXd x5 = Am * x4;
		if (verbose_level > 1)
		{
			cout << "x5: " << x5.rows() << ',' << x5.cols() << endl;
			if (verbose_level > 2)
				save_mat("x5.dat", x5);
		}

		message(1, "forming X6");
		performance_log->log_event("forming x6");
		Eigen::MatrixXd x6 = par_diff.transpose() * x5;
		if (verbose_level > 1)
		{
			cout << "x6: " << x6.rows() << ',' << x6.cols() << endl;
			if (verbose_level > 2)
				save_mat("x6.dat", x6);
		}

		message(1, "forming X7");
		performance_log->log_event("forming x7");
		if (verbose_level > 1)
		{
			cout << "V: " << V.rows() << ',' << V.cols() << endl;
			if (verbose_level > 2)
				save_mat("V.dat", V);
		}
		Eigen::MatrixXd x7 = V * ivec *V.transpose() * x6;
		if (verbose_level > 1)
		{
			cout << "x7: " << x7.rows() << ',' << x7.cols() << endl;
			if (verbose_level > 2)
				save_mat("x7.dat", x7);
		}

		performance_log->log_event("forming upgrade_2");
		message(1, "forming upgrade_2");

		if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
		{
			upgrade_2 = -1.0 * parcov_inv * par_diff * x7;
		}
		else
		{
			upgrade_2 = -1.0 * (par_diff * x7);
		}

		if (verbose_level > 1)
		{
			cout << "upgrade_2: " << upgrade_2.rows() << ',' << upgrade_2.cols() << endl;
			if (verbose_level > 2)
				save_mat("upgrade_2", upgrade_2);
		}
		upgrade_1 = upgrade_1 + upgrade_2.transpose();
		//pe_lam.set_eigen(*pe_lam.get_eigen_ptr() + upgrade_2.transpose());

	}
	ParameterEnsemble upgrade_pe(&pest_scenario);
	upgrade_pe.from_eigen_mat(upgrade_1,pe.get_real_names(), par_names);
	return upgrade_pe;
}


bool IterEnsembleSmoother::solve_new()
{
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();
	if (pe.shape().first <= error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need at least ", error_min_reals);
		throw_ies_error(string("too few active realizations, cannot continue"));
	}
	if (pe.shape().first <= warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(1, s);
	}

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

	
	vector<ParameterEnsemble> pe_lams;
	vector<double> lam_vals, scale_vals;
	for (auto &lam_mult : lam_mults)
	{
		ss.str("");
		double cur_lam = last_best_lam * lam_mult;
		ss << "starting calcs for lambda" << cur_lam;
		message(0, "starting lambda calcs for lambda", cur_lam);
		ParameterEnsemble pe_upgrade = pe.zeros_like();
		if (use_localizer)
		{
			int i = 0;
			int lsize = localizer.get_localizer_map().size();

			for (auto local_pair : localizer.get_localizer_map())
			{
				ss.str("");
				ss << "localized upgrade part " << i + 1 << " of " << lsize;
				message(1, ss.str());
				ParameterEnsemble pe_local = calc_upgrade(local_pair.first, local_pair.second, cur_lam, pe.shape().first);
				pe_upgrade.add_2_cols_ip(pe_local);
				i++;
			}
		}
		else
			pe_upgrade = calc_upgrade(act_obs_names, act_par_names, cur_lam, pe.shape().first);

		for (auto sf : pest_scenario.get_pestpp_options().get_lambda_scale_vec())
		{

			ParameterEnsemble pe_lam_scale = pe;
			pe_lam_scale.set_eigen(*pe_lam_scale.get_eigen_ptr() + (*pe_upgrade.get_eigen_ptr() * sf));
			if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds())
				pe_lam_scale.enforce_bounds();
			pe_lams.push_back(pe_lam_scale);
			lam_vals.push_back(cur_lam);
			scale_vals.push_back(sf);
			if (!pest_scenario.get_pestpp_options().get_ies_save_lambda_en())
				continue;
			ss.str("");
			ss << file_manager.get_base_filename() << "." << iter << "." << cur_lam << ".lambda." << sf << ".scale.par";

			if (pest_scenario.get_pestpp_options().get_ies_save_binary())
			{
				ss << ".jcb";
				pe_lam_scale.to_binary(ss.str());
			}
			else
			{
				ss << ".csv";
				pe_lam_scale.to_csv(ss.str());
			}
			frec << "lambda, scale value " << cur_lam << ',' << sf << " pars saved to " << ss.str() << endl;

		}

		ss.str("");
		message(1, "finished calcs for lambda:", cur_lam);

	}
	//return;
	vector<map<int, int>> real_run_ids_lams;
	int best_idx = -1;
	double best_mean = 1.0e+30, best_std = 1.0e+30;
	double mean, std;

	message(0, "running lambda ensembles");
	vector<ObservationEnsemble> oe_lams = run_lambda_ensembles(pe_lams, lam_vals, scale_vals);

	message(0, "evaluting lambda ensembles");
	message(1, "last mean: ", last_best_mean);
	message(1, "last stdev: ", last_best_std);

	ObservationEnsemble oe_lam_best;
	for (int i = 0; i<pe_lams.size(); i++)
	{
		if (oe_lams[i].shape().first == 0)
			continue;
		vector<double> vals({ lam_vals[i],scale_vals[i] });
		drop_bad_phi(pe_lams[i], oe_lams[i]);
		if (oe_lams[i].shape().first == 0)
		{
			message(1, "all realizations dropped as 'bad' for lambda, scale fac ", vals);
			continue;
		}
		message(0, "phi summary for lambda, scale fac:", vals);
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
	if (best_idx == -1)
	{
		message(0, "WARNING:  unsuccessful lambda testing, resetting lambda to 10000.0");
		last_best_lam = 10000.0;
		return false;

	}
	double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
	double lam_inc = pest_scenario.get_pestpp_options().get_ies_lambda_inc_fac();
	double lam_dec = pest_scenario.get_pestpp_options().get_ies_lambda_dec_fac();


	//subset stuff here
	if ((best_idx != -1) && (use_subset) && (subset_size < pe.shape().first))
	{

		double acc_phi = last_best_mean * acc_fac;
		if (best_mean > acc_phi)
		{
			double new_lam = last_best_lam * lam_inc;
			new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
			last_best_lam = new_lam;
			ss.str("");
			ss << "best subset mean phi  (" << best_mean << ") greater than acceptable phi : " << acc_phi;
			string m = ss.str();
			message(0, m);
			message(1, "abandoning current lambda ensembles, increasing lambda to ", new_lam);
			message(1, "returing to lambda calculations...");
			return false;
		}


		//release the memory of the unneeded pe_lams
		for (int i = 0; i < pe_lams.size(); i++)
		{
			if (i == best_idx)
				continue;
			pe_lams[i] = ParameterEnsemble();
		}
		//need to work out which par and obs en real names to run - some may have failed during subset testing...
		ObservationEnsemble remaining_oe_lam = oe;//copy
		ParameterEnsemble remaining_pe_lam = pe_lams[best_idx];
		vector<string> pe_keep_names, oe_keep_names;
		vector<string> pe_names = pe.get_real_names(), oe_names = oe.get_real_names();

		vector<string> org_pe_idxs,org_oe_idxs;
		set<string> ssub;
		for (auto &i : subset_idxs)
			ssub.emplace(pe_names[i]);
		for (int i=0;i<pe_names.size();i++)
			if (ssub.find(pe_names[i]) == ssub.end())
			{
				pe_keep_names.push_back(pe_names[i]);
				oe_keep_names.push_back(oe_names[i]);
			}
		message(0, "running remaining realizations for best lambda, scale:", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));

		//pe_keep_names and oe_keep_names are names of the remaining reals to eval
		remaining_pe_lam.keep_rows(pe_keep_names);
		remaining_oe_lam.keep_rows(oe_keep_names);
		//save these names for later
		org_pe_idxs = remaining_pe_lam.get_real_names();
		org_oe_idxs = remaining_oe_lam.get_real_names();
		///run
		vector<int> fails = run_ensemble(remaining_pe_lam, remaining_oe_lam);

		//for testing
		//fails.push_back(0);

		//if any of the remaining runs failed
		if (fails.size() == org_pe_idxs.size())
			throw_ies_error(string("all remaining realizations failed...something is prob wrong"));
		if (fails.size() > 0)
		{

			vector<string> new_pe_idxs, new_oe_idxs;
			vector<int>::iterator start = fails.begin(), end = fails.end();
			stringstream ss;
			ss << "the following par:obs realizations failed during evaluation of the remaining ensemble";
			for (int i = 0; i < org_pe_idxs.size(); i++)
				if (find(start, end, i) == end)
				{
					new_pe_idxs.push_back(org_pe_idxs[i]);
					new_oe_idxs.push_back(org_oe_idxs[i]);
				}
				else
				{
					ss << org_pe_idxs[i] << ":" << org_oe_idxs[i] << " , ";
				}
			string s = ss.str();
			message(1, s);
			remaining_oe_lam.keep_rows(new_oe_idxs);
			remaining_pe_lam.keep_rows(new_pe_idxs);

		}
		//drop the remaining runs from the par en then append the remaining par runs (in case some failed)
		performance_log->log_event("assembling ensembles");
		pe_lams[best_idx].drop_rows(pe_keep_names);
		pe_lams[best_idx].append_other_rows(remaining_pe_lam);
		//append the remaining obs en
		oe_lam_best.append_other_rows(remaining_oe_lam);
		assert(pe_lams[best_idx].shape().first == oe_lam_best.shape().first);
		drop_bad_phi(pe_lams[best_idx], oe_lam_best);
		if (oe_lam_best.shape().first == 0)
		{
			throw_ies_error(string("all realization dropped after finishing subset runs...something might be wrong..."));
		}
		performance_log->log_event("updating phi");
		ph.update(oe_lam_best, pe_lams[best_idx]);
		best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
		best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
		message(1, "phi summary for entire ensemble using lambda,scale_fac ", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));
		ph.report();
	}
	else
	{
		ph.update(oe_lam_best, pe_lams[best_idx]);
		best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
		best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
	}

	ph.update(oe_lam_best, pe_lams[best_idx]);
	best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
	best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
	message(1, "last best mean phi * acceptable phi factor: ", last_best_mean * acc_fac);
	message(1, "current best mean phi: ", best_mean);

	//track this here for phi-based termination check
	best_mean_phis.push_back(best_mean);

	if (best_mean < last_best_mean * acc_fac)
	{
		message(0, "updating parameter ensemble");
		performance_log->log_event("updating parameter ensemble");
		last_best_mean = best_mean;

		pe = pe_lams[best_idx];
		oe = oe_lam_best;
		if (best_std < last_best_std * acc_fac)
		{
			double new_lam = lam_vals[best_idx] * lam_dec;
			new_lam = (new_lam < lambda_min) ? lambda_min : new_lam;
			message(0, "updating lambda to ", new_lam);
			last_best_lam = new_lam;
		}
		else
		{
			message(0, "not updating lambda");
		}
		last_best_std = best_std;
	}

	else
	{
		message(0, "not updating parameter ensemble");
		ph.update(oe, pe);
		double new_lam = last_best_lam * lam_inc;
		new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
		message(0, "incresing lambda to: ", new_lam);
		last_best_lam = new_lam;
	}
	//report_and_save();
	return true;
}

//
//bool IterEnsembleSmoother::solve_old()
//{
//	stringstream ss;
//	ofstream &frec = file_manager.rec_ofstream();
//	if (pe.shape().first <= error_min_reals)
//	{
//		message(0, "too few active realizations:", oe.shape().first);
//		message(1, "need at least ", error_min_reals);
//		throw_ies_error(string("too few active realizations, cannot continue"));
//	}
//	if (pe.shape().first <= warn_min_reals)
//	{
//		ss.str("");
//		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
//		string s = ss.str();
//		message(1, s);
//	}
//
//	if ((use_subset) && (subset_size > pe.shape().first))
//	{
//		ss.str("");
//		ss << "++ies_subset size (" << subset_size << ") greater than ensemble size (" << pe.shape().first << ")";
//		frec << "  ---  " << ss.str() << endl;
//		cout << "  ---  " << ss.str() << endl;
//		frec << "  ...reducing ++ies_subset_size to " << pe.shape().first << endl;
//		cout << "  ...reducing ++ies_subset_size to " << pe.shape().first << endl;
//		subset_size = pe.shape().first;
//	}
//
//	double scale = (1.0 / (sqrt(double(oe.shape().first - 1))));
//	
//	performance_log->log_event("calculate residual matrix");
//	//oe_base var_names should be ordered by act_obs_names, so only reorder real_names
//	//oe should only include active realizations, so only reorder var_names
//	message(1, "calculating residual matrix");
//	Eigen::MatrixXd scaled_residual = obscov_inv_sqrt * ph.get_obs_resid(oe).transpose();
//	if (verbose_level > 1)
//	{
//		cout << "scaled_residual: " << scaled_residual.rows() << ',' << scaled_residual.cols() << endl;
//		if (verbose_level > 2)
//		{
//			save_mat("scaled_residual.dat", scaled_residual);
//			Eigen::MatrixXd residual = ph.get_obs_resid(oe).transpose();
//			save_mat("residual.dat",residual);
//		}
//	}
//		
//	performance_log->log_event("calculate scaled obs diff");
//	message(1, "calculating obs diff matrix");
//	Eigen::MatrixXd diff = oe.get_eigen_mean_diff(vector<string>(),act_obs_names).transpose();
//	Eigen::MatrixXd obs_diff = scale * (obscov_inv_sqrt * diff);
//	if (verbose_level > 1)
//	{
//		cout << "obs_diff: " << obs_diff.rows() << ',' << obs_diff.cols() << endl;
//		if (verbose_level > 2)
//			save_mat("obs_diff.dat", obs_diff);
//	}
//
//	performance_log->log_event("calculate scaled par diff");
//	message(1, "calculating par diff matrix");
//	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
//	diff = pe.get_eigen_mean_diff(vector<string>(), act_par_names).transpose();
//	Eigen::MatrixXd par_diff;
//	if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
//	{
//		cout << "...applying prior par cov scaling to par diff matrix" << endl;
//		par_diff = scale * parcov_inv_sqrt * diff;
//	}
//	else
//		par_diff = scale * diff;
//		
//	if (verbose_level > 1)
//	{
//		cout << "par_diff:" << par_diff.rows() << ',' << par_diff.cols() << endl;
//		if (verbose_level > 2)
//		{
//			save_mat("scaled_par_diff.dat", par_diff);
//			save_mat("par_diff.dat", diff);
//		}
//	}
//
//	performance_log->log_event("SVD of obs diff");
//	message(1, "calculating SVD of obs diff matrix");
//	Eigen::MatrixXd ivec, upgrade_1, s,V,Ut;
//
//	SVD_REDSVD rsvd;
//	rsvd.set_performance_log(performance_log);
//	rsvd.solve_ip(obs_diff, s, Ut, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
//	
//	//SVD_EIGEN esvd;
//	//esvd.set_performance_log(performance_log);
//	//esvd.solve_ip(obs_diff, s, Ut, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
//
//	Ut.transposeInPlace();
//	obs_diff.resize(0, 0);
//	
//	Eigen::MatrixXd s2 = s.cwiseProduct(s);
//	if (verbose_level > 1)
//	{
//		cout << "s2: " << s2.rows() << ',' << s2.cols() << endl;
//		cout << "Ut: " << Ut.rows() << ',' << Ut.cols() << endl;
//		cout << "V:" << V.rows() << ',' << V.cols() << endl;
//		if (verbose_level > 2)
//		{
//			save_mat("ut.dat", Ut);
//			save_mat("s2.dat", s2);
//		}
//	}
//
//	vector<ParameterEnsemble> pe_lams;
//	vector<double> lam_vals, scale_vals;
//	for (auto &lam_mult : lam_mults)
//	{
//		ss.str("");
//		double cur_lam = last_best_lam * lam_mult;
//		ss << "starting calcs for lambda" << cur_lam;
//		message(0, "starting lambda calcs for lambda", cur_lam);
//
//		performance_log->log_event(ss.str());
//		performance_log->log_event("form scaled identity matrix");
//		message(1, "calculating scaled identity matrix");
//		ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();
//		if (verbose_level > 1)
//		{
//			cout << "ivec:" << ivec.rows() << ',' << ivec.cols() << endl;
//			if (verbose_level > 2)
//				save_mat("ivec.dat", ivec);
//		}
//
//
//		message(1, "forming X1");
//		Eigen::MatrixXd X1 = Ut * scaled_residual;
//		if (verbose_level > 1)
//		{
//			cout << "X1: " << X1.rows() << ',' << X1.cols() << endl;
//			if (verbose_level > 2)
//				save_mat("X1.dat", X1);
//		}
//		//scaled_residual.resize(0, 0);
//		//Ut.resize(0, 0);
//
//		message(1, "forming X2");
//		Eigen::MatrixXd X2 = ivec * X1;
//		if (verbose_level > 1)
//		{
//			cout << "X2: " << X2.rows() << ',' << X2.cols() << endl;
//			if (verbose_level > 2)
//				save_mat("X2.dat", X2);
//		}
//		X1.resize(0, 0);
//
//		message(1, "forming X3");
//		Eigen::MatrixXd X3 = V * s.asDiagonal() * X2;
//		if (verbose_level > 1)
//		{
//			cout << "X3: " << X3.rows() << ',' << X3.cols() << endl;
//			if (verbose_level > 2)
//				save_mat("X3.dat", X3);
//		}
//		X2.resize(0, 0);
//
//		message(1, "forming upgrade_1");
//		if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
//		{
//			upgrade_1 = -1.0 * parcov_inv_sqrt * par_diff * X3;
//		}
//		else
//		{
//			upgrade_1 = -1.0 * par_diff * X3;
//		}
//		upgrade_1.transposeInPlace();
//		if (verbose_level > 1)
//		{
//			cout << "upgrade_1:" << upgrade_1.rows() << ',' << upgrade_1.cols() << endl;
//			if (verbose_level > 2)
//				save_mat("upgrade_1.dat", upgrade_1);
//		}
//		X3.resize(0,0);
//
//		/*ParameterEnsemble pe_lam = pe;
//		pe_lam.set_eigen(*pe_lam.get_eigen_ptr() + upgrade_1);
//		upgrade_1.resize(0, 0);*/
//		
//		Eigen::MatrixXd upgrade_2;
//		if ((!pest_scenario.get_pestpp_options().get_ies_use_approx()) && (iter > 1))
//		{
//			performance_log->log_event("calculating parameter correction (full solution)");
//			message(1, "calculating parameter correction (full solution, MAP)");
//			performance_log->log_event("forming scaled par resid");
//			Eigen::MatrixXd scaled_par_resid;
//			if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
//			{
//				scaled_par_resid = parcov_inv_sqrt *  ph.get_par_resid(pe).transpose();
//			}
//			else
//			{
//				scaled_par_resid = ph.get_par_resid(pe).transpose();
//			}
//
//			//scaled_par_resid.transposeInPlace();
//			
//			performance_log->log_event("forming x4");
//			message(1, "forming X4");
//			if (verbose_level > 1)
//			{
//				cout << "scaled_par_resid: " << scaled_par_resid.rows() << ',' << scaled_par_resid.cols() << endl;
//				if (verbose_level > 2)
//					save_mat("scaled_par_resid.dat", scaled_par_resid);
//			}
//			Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
//			if (verbose_level > 1)
//			{
//				cout << "x4: " << x4.rows() << ',' << x4.cols() << endl;
//				if (verbose_level > 2)
//					save_mat("x4.dat", x4);
//			}
//
//			performance_log->log_event("forming x5");
//			message(1, "forming X5");
//			Eigen::MatrixXd x5 = Am * x4;
//			if (verbose_level > 1)
//			{
//				cout << "x5: " << x5.rows() << ',' << x5.cols() << endl;
//				if (verbose_level > 2)
//					save_mat("x5.dat", x5);
//			}
//			
//			message(1, "forming X6");
//			performance_log->log_event("forming x6");
//			Eigen::MatrixXd x6 = par_diff.transpose() * x5;
//			if (verbose_level > 1)
//			{
//				cout << "x6: " << x6.rows() << ',' << x6.cols() << endl;
//				if (verbose_level > 2)
//					save_mat("x6.dat", x6);
//			}
//				
//			message(1, "forming X7");
//			performance_log->log_event("forming x7");
//			if (verbose_level > 1)
//			{
//				cout << "V: " << V.rows() << ',' << V.cols() << endl;
//				if (verbose_level > 2)
//					save_mat("V.dat", V);
//			}	
//			Eigen::MatrixXd x7 = V * ivec *V.transpose() * x6;
//			if (verbose_level > 1)
//			{
//				cout << "x7: " << x7.rows() << ',' << x7.cols() << endl;
//				if (verbose_level > 2)
//					save_mat("x7.dat", x7);
//			}
//
//			performance_log->log_event("forming upgrade_2");
//			message(1, "forming upgrade_2");
//			
//			if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
//			{
//				upgrade_2 = -1.0 * parcov_inv_sqrt * par_diff * x7;
//			}
//			else
//			{
//				upgrade_2 = -1.0 * (par_diff * x7);
//			}
//			
//			if (verbose_level > 1)
//			{
//				cout << "upgrade_2: " << upgrade_2.rows() << ',' << upgrade_2.cols() << endl;
//				if (verbose_level > 2)
//					save_mat("upgrade_2", upgrade_2);
//			}
//
//			//pe_lam.set_eigen(*pe_lam.get_eigen_ptr() + upgrade_2.transpose());
//		}
//
//		for (auto sf : pest_scenario.get_pestpp_options().get_lambda_scale_vec())
//		{
//
//			ParameterEnsemble pe_lam_scale = pe;
//			if (upgrade_2.rows() > 0)
//				upgrade_1 = upgrade_1 + upgrade_2.transpose();
//			pe_lam_scale.set_eigen(*pe_lam_scale.get_eigen_ptr() + (upgrade_1 * sf));
//			if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds())
//				pe_lam_scale.enforce_bounds();
//			pe_lams.push_back(pe_lam_scale);
//			lam_vals.push_back(cur_lam);
//			scale_vals.push_back(sf);
//			if (!pest_scenario.get_pestpp_options().get_ies_save_lambda_en())
//				continue;
//			ss.str("");
//			ss << file_manager.get_base_filename() << "." << iter << "." << cur_lam << ".lambda." << sf << ".scale.par";
//
//			if (pest_scenario.get_pestpp_options().get_ies_save_binary())
//			{
//				ss << ".jcb";
//				pe_lam_scale.to_binary(ss.str());
//			}
//			else
//			{
//				ss << ".csv";
//				pe_lam_scale.to_csv(ss.str());
//			}
//			frec << "lambda, scale value " << cur_lam << ',' << sf << " pars saved to " << ss.str() << endl;
//
//		}
//
//		
//		ss.str("");
//		message(1, "finished calcs for lambda:", cur_lam);
//
//	}
//	//return;
//	vector<map<int, int>> real_run_ids_lams;
//	int best_idx = -1;
//	double best_mean = 1.0e+30, best_std = 1.0e+30;
//	double mean, std;
//	
//	message(0, "running lambda ensembles");
//	vector<ObservationEnsemble> oe_lams = run_lambda_ensembles(pe_lams, lam_vals,scale_vals);
//		
//	message(0, "evaluting lambda ensembles");
//	message(1, "last mean: ", last_best_mean);
//	message(1, "last stdev: ", last_best_std);
//
//	ObservationEnsemble oe_lam_best;
//	for (int i=0;i<pe_lams.size();i++)
//	{	
//		//for testing...
//		//pest_scenario.get_pestpp_options_ptr()->set_ies_bad_phi(0.0);
//		//if all runs failed...
//		if (oe_lams[i].shape().first == 0)
//			continue;
//		vector<double> vals({ lam_vals[i],scale_vals[i] });
//		drop_bad_phi(pe_lams[i], oe_lams[i]);
//		if (oe_lams[i].shape().first == 0)
//		{
//			message(1, "all realizations dropped as 'bad' for lambda, scale fac ",vals);
//			continue;
//		}
//		message(0, "phi summary for lambda, scale fac:", vals);
//		ph.update(oe_lams[i], pe_lams[i]);
//		ph.report();
//		mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
//		std = ph.get_std(PhiHandler::phiType::COMPOSITE);
//		if (mean < best_mean)
//		{
//			oe_lam_best = oe_lams[i];
//			best_mean = mean;
//			best_std = std;
//			best_idx = i;
//		}
//	}
//	if (best_idx == -1)
//	{
//		message(0, "WARNING:  unsuccessful lambda testing, resetting lambda to 10000.0");
//		last_best_lam = 10000.0;
//		return false;
//
//	}
//	double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
//	double lam_inc = pest_scenario.get_pestpp_options().get_ies_lambda_inc_fac();
//	double lam_dec = pest_scenario.get_pestpp_options().get_ies_lambda_dec_fac();
//
//	
//	//subset stuff here
//	if ((best_idx != -1) && (use_subset) && (subset_size < pe.shape().first))
//	{ 
//		
//		double acc_phi = last_best_mean * acc_fac;
//		if (best_mean > acc_phi)
//		{
//			double new_lam = last_best_lam * lam_inc;
//			new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
//			last_best_lam = new_lam;
//			ss.str("");
//			ss << "best subset mean phi  (" << best_mean << ") greater than acceptable phi : " << acc_phi;
//			string m = ss.str();
//			message(0, m);
//			message(1, "abandoning current lambda ensembles, increasing lambda to ", new_lam);
//			message(1,"returing to lambda calculations...");
//			return false;
//		}
//
//
//		//release the memory of the unneeded pe_lams
//		for (int i = 0; i < pe_lams.size(); i++)
//		{
//			if (i == best_idx)
//				continue;
//			pe_lams[i] = ParameterEnsemble();
//		}
//		//need to work out which par and obs en real names to run - some may have failed during subset testing...
//		ObservationEnsemble remaining_oe_lam = oe;//copy
//		ParameterEnsemble remaining_pe_lam = pe_lams[best_idx];
//		vector<string> pe_keep_names, oe_keep_names;
//		vector<string> pe_names = pe.get_real_names(), oe_names = oe.get_real_names();
//		vector<string> org_pe_idxs,org_oe_idxs;
//		for (int i = subset_size; i <pe.shape().first; i++)
//		{	
//			pe_keep_names.push_back(pe_names[i]);
//			oe_keep_names.push_back(oe_names[i]);
//		}
//		message(0, "running remaining realizations for best lambda, scale:", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));
//
//		//pe_keep_names and oe_keep_names are names of the remaining reals to eval
//		remaining_pe_lam.keep_rows(pe_keep_names);
//		remaining_oe_lam.keep_rows(oe_keep_names);
//		//save these names for later
//		org_pe_idxs = remaining_pe_lam.get_real_names();
//		org_oe_idxs = remaining_oe_lam.get_real_names();
//		///run
//		vector<int> fails = run_ensemble(remaining_pe_lam, remaining_oe_lam);
//		//if any of the remaining runs failed
//		if (fails.size() ==org_pe_idxs.size())
//			throw_ies_error(string("all remaining realizations failed...something is prob wrong"));
//		if (fails.size() > 0)
//		{
//
//			vector<string> new_pe_idxs, new_oe_idxs;
//			vector<int>::iterator start = fails.begin(), end = fails.end();
//			stringstream ss;
//			ss << "the following par:obs realizations failed during evaluation of the remaining ensemble";
//			for (int i = 0; i < org_pe_idxs.size(); i++)
//				if (find(start,end,i) == end)
//				{
//					new_pe_idxs.push_back(org_pe_idxs[i]);
//					new_oe_idxs.push_back(org_oe_idxs[i]);
//				}
//				else
//				{
//					ss << org_pe_idxs[i] << ":" << org_oe_idxs[i] << " , ";
//				}
//			string s = ss.str();
//			message(1, s);
//			remaining_oe_lam.keep_rows(new_oe_idxs);
//			remaining_pe_lam.keep_rows(new_pe_idxs);
//
//		}
//		//drop the remaining runs from the par en then append the remaining par runs (in case some failed)
//		performance_log->log_event("assembling ensembles");
//		pe_lams[best_idx].drop_rows(pe_keep_names); 
//		pe_lams[best_idx].append_other_rows(remaining_pe_lam);
//		//append the remaining obs en
//		oe_lam_best.append_other_rows(remaining_oe_lam);
//		assert(pe_lams[best_idx].shape().first == oe_lam_best.shape().first);
//		drop_bad_phi(pe_lams[best_idx], oe_lam_best);
//		if (oe_lam_best.shape().first == 0)
//		{
//			throw_ies_error(string("all realization dropped after finishing subset runs...something might be wrong..."));
//		}
//		performance_log->log_event("updating phi");
//		ph.update(oe_lam_best, pe_lams[best_idx]);
//		best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
//		best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
//		message(1, "phi summary for entire ensemble using lambda,scale_fac ", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));
//		ph.report();
//	}
//	else
//	{ 
//		ph.update(oe_lam_best, pe_lams[best_idx]);
//		best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
//		best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
//	}
//	
//	ph.update(oe_lam_best, pe_lams[best_idx]);
//	best_mean = ph.get_mean(PhiHandler::phiType::COMPOSITE);
//	best_std = ph.get_std(PhiHandler::phiType::COMPOSITE);
//	message(1, "last best mean phi * acceptable phi factor: ", last_best_mean * acc_fac);
//	message(1, "current best mean phi: ", best_mean);
//
//	//track this here for phi-based termination check
//	best_mean_phis.push_back(best_mean);
//
//	if (best_mean < last_best_mean * acc_fac)
//	{
//		message(0,"updating parameter ensemble");
//		performance_log->log_event("updating parameter ensemble");
//		last_best_mean = best_mean;
//
//		pe = pe_lams[best_idx];
//		oe = oe_lam_best;		
//		if (best_std < last_best_std * acc_fac)
//		{
//			double new_lam = lam_vals[best_idx] * lam_dec;
//			new_lam = (new_lam < lambda_min) ? lambda_min : new_lam;
//			message(0, "updating lambda to ", new_lam);
//			last_best_lam = new_lam;
//		}
//		else
//		{
//			message(0, "not updating lambda");
//		}
//		last_best_std = best_std;
//	}
//
//	else
//	{
//		message(0, "not updating parameter ensemble");
//		ph.update(oe, pe);
//		double new_lam = last_best_lam * lam_inc;
//		new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
//		message(0, "incresing lambda to: ", new_lam);
//		last_best_lam = new_lam;
//	}
//	//report_and_save();
//	return true;
//}

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
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << "." << iter << ".obs.jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << "." << iter << ".obs.csv";
		oe.to_csv(ss.str());
	}
	frec << "      current obs ensemble saved to " << ss.str() << endl;
	cout << "      current obs ensemble saved to " << ss.str() << endl;
	ss.str("");
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << "." << iter << ".par.jcb";
		pe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
		pe.to_csv(ss.str());
	}
	//ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
	//pe.to_csv(ss.str());
	frec << "      current par ensemble saved to " << ss.str() << endl;
	cout << "      current par ensemble saved to " << ss.str() << endl;
	
	
	
}


void IterEnsembleSmoother::set_subset_idx(int size)
{
	//map<int,int> subset_idx_map;
	subset_idxs.clear();
	int nreal_subset = pest_scenario.get_pestpp_options().get_ies_subset_size();
	if ((!use_subset) || (nreal_subset >= size))
		return;
	//int size = pe.shape().first;
	string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
	if (how == "FIRST")
	{
		for (int i = 0; i < nreal_subset; i++)
			subset_idxs.push_back(i);
	}
	else if (how == "LAST")
	{
		
		for (int i = size - nreal_subset; i < size; i++)
			subset_idxs.push_back(i);
	}
	
	else if (how == "RANDOM")
	{
		std::uniform_int_distribution<int> uni(0, size);
		for (int i = 0; i < nreal_subset; i++)
			subset_idxs.push_back(uni(Ensemble::rand_engine));
	}
	else if (how == "PHI_BASED")
	{
		//sidx needs to be index of realization, not realization number
		vector<pair<double, int>> phis;
		//vector<int> sidx;
		int step;
		PhiHandler::phiType pt = PhiHandler::phiType::COMPOSITE;
		map<string, double>* phi_map = ph.get_phi_map(pt);
		map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();

		int i = 0;
		for (pi; pi != end; ++pi)
		{
			phis.push_back(make_pair(pi->second, i)); //phival,idx?
			++i;
		}
		sort(phis.begin(), phis.end());

		//include idx for lowest and highest phi reals
		subset_idxs.push_back(phis[0].second);
		subset_idxs.push_back(phis[phis.size() - 1].second);//can't get .back() to work with vector<pair<>>

													 //start near middle, work out
													 //start with mid, moving up and down gives more than enough
		int mid = phis.size() / 2;
		step = phis.size() / nreal_subset;
		for (i = 0; i < nreal_subset; ++i)
		{
			//add higher phis first
			if ((subset_idxs.size() < nreal_subset) && (find(subset_idxs.begin(), subset_idxs.end(), phis[mid + i * step].second) == subset_idxs.end()))
			{
				// slopy but otherwise hard to ensure nsubsets and good spread
				subset_idxs.push_back(phis[mid + i * step].second);
			}
			if ((subset_idxs.size() < nreal_subset) && (find(subset_idxs.begin(), subset_idxs.end(), phis[mid - i * step].second) == subset_idxs.end()))
			{
				// slopy but otherwise hard to ensure nsubsets and good spread
				subset_idxs.push_back(phis[mid - i * step].second);
			}
		}
	}
	else
	{
		//throw runtime_error("unkonwn 'subset_how'");
		throw_ies_error("unknown 'subset_how'");
	}
	message(1,"subset idxs:",subset_idxs);
	return;
	//return subset_idx_map;
}

vector<ObservationEnsemble> IterEnsembleSmoother::run_lambda_ensembles(vector<ParameterEnsemble> &pe_lams, vector<double> &lam_vals, vector<double> &scale_vals)
{
	ofstream &frec = file_manager.rec_ofstream();
	stringstream ss;
	ss << "queuing " << pe_lams.size() << " ensembles";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	/*vector<int> subset_idxs;
	if ((use_subset) && (subset_size < pe_lams[0].shape().first))
	{
		for (int i = 0; i < subset_size; i++)
			subset_idxs.push_back(i);
	}*/
	set_subset_idx(pe_lams[0].shape().first);
	vector<map<int, int>> real_run_ids_vec;
	//ParameterEnsemble pe_lam;
	//for (int i=0;i<pe_lams.size();i++)
	for (auto &pe_lam : pe_lams)
	{
		// if (pest_scenario.get_pestpp_options().get_ies_save_binary())
		// {
		// 	ss << file_manager.get_base_filename() << "." << iter << "." << i << ".lambdapars.jcb";
		// 	pe_lam.to_binary(ss.str());
		// }
		// else
		// {
		// 	ss << file_manager.get_base_filename() << "." << iter << "." << i << ".lambdapars.jcb";
		// }
		// frec << " lambda value " << lam_vals[i] << " pars saved to " << ss.str() << endl;
		// ss.str("");
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
	//ObservationEnsemble _oe = oe;//copy
	//if (subset_size < pe_lams[0].shape().first)
	//	_oe.keep_rows(subset_idxs);
	map<int, int> real_run_ids;
	//for (auto &real_run_ids : real_run_ids_vec)
	for (int i=0;i<pe_lams.size();i++)
	{
		ObservationEnsemble _oe = oe;//copy
		vector<double> rep_vals{ lam_vals[i],scale_vals[i] };
		real_run_ids = real_run_ids_vec[i];
		//if using subset, reset the real_idx in real_run_ids to be just simple counter
		if (subset_size < pe_lams[0].shape().first)
		{
			_oe.keep_rows(subset_idxs);
			int ireal = 0;
			map<int, int> temp;
			for (auto &rri : real_run_ids)
			{
				temp[ireal] = rri.second;
				ireal++;
			}

			real_run_ids = temp;
		}
		
		try
		{
			failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr);
		}
		catch (const exception &e)
		{
			stringstream ss;
			ss << "error processing runs for lambda,scale: " << lam_vals[i] << ',' << scale_vals[i] << ':' << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			stringstream ss;
			ss << "error processing runs for lambda,scale: " << lam_vals[i] << ',' << scale_vals[i];
			throw_ies_error(ss.str());
		}
		
		//for testing
		//failed_real_indices.push_back(real_run_ids.size()-1);
		
		if (failed_real_indices.size() > 0)
		{
			stringstream ss;
			vector<string> par_real_names = pe.get_real_names();
			vector<string> obs_real_names = oe.get_real_names();
			vector<string> failed_par_names, failed_obs_names;
			string oname, pname;
			ss << "the following par:obs realization runs failed for lambda,scale " << lam_vals[i] << ',' << scale_vals[i] << "-->";
			for (auto &i : failed_real_indices)
			{
				pname = par_real_names[subset_idxs[i]];
				oname = obs_real_names[subset_idxs[i]];
				failed_par_names.push_back(pname);
				failed_obs_names.push_back(oname);
				ss << pname << ":" << oname << ',';
			}
			string s = ss.str();
			message(1,s);
			if (failed_real_indices.size() == _oe.shape().first)
			{
				message(0, "WARNING: all realizations failed for lambda, scale :", rep_vals);
				_oe = ObservationEnsemble();

			}
			else
			{
				performance_log->log_event("dropping failed realizations");
				//_oe.drop_rows(failed_real_indices);
				//pe_lams[i].drop_rows(failed_real_indices);
				_oe.drop_rows(failed_obs_names);
				pe_lams[i].drop_rows(failed_par_names);
			}
			
		}
		obs_lams.push_back(_oe);
	}
	return obs_lams;
}


vector<int> IterEnsembleSmoother::run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs)
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
	//for testing
	//failed_real_indices.push_back(0);
	
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
		message(1, "failed realizations: ", failed_real_indices.size());
		string s = ss.str();
		message(1, s);
		performance_log->log_event("dropping failed realizations");
		_pe.drop_rows(failed_real_indices);
		_oe.drop_rows(failed_real_indices);
	}
	return failed_real_indices;
}


void IterEnsembleSmoother::finalize()
{

}
