/*  
    © Copyright 2012, David Welter
    
    This file is part of PEST++.
   
    PEST++ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PEST++ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/
#include <fstream>
#include <iostream>
#include <lapackpp.h>
#include <iomanip>
#include <map>
#include <algorithm>
#include "SVDSolver.h"
#include "RunManagerAbstract.h"
#include "QSqrtMatrix.h"
#include "lapack_tools.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "PriorInformation.h"
#include "Regularization.h"
#include "SVD_PROPACK.h"


using namespace std;
using namespace pest_utils;


void SVDSolver::Upgrade::print(ostream &os) 
{
	os << "Upgrade: " << n_sing_val_used << "out of " <<  tot_sing_val << "singular vales used" << endl;
	int vec_size = par_name_vec.size();
	for (int i=0; i< vec_size; ++i) {
		os << setw(21) << par_name_vec[i] << setw(20) << svd_uvec(i) << setw(20) <<grad_uvec(i) << endl;
	}
}

SVDSolver::SVDSolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_par_group_info_ptr, const ParameterInfo *_ctl_par_info_ptr,
		const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, 
		const Regularization *_regul_scheme_ptr, const string &_description)
		: ctl_info(_ctl_info), svd_info(_svd_info), par_group_info_ptr(_par_group_info_ptr), ctl_par_info_ptr(_ctl_par_info_ptr), obs_info_ptr(_obs_info), obj_func(_obj_func),
		  file_manager(_file_manager), observations_ptr(_observations), par_transform(_par_transform),
		  cur_solution(_obj_func, _par_transform, *_observations), phiredswh_flag(false), save_next_jacobian(true), prior_info_ptr(_prior_info_ptr), jacobian(_jacobian), prev_phi_percent(0.0),
		  num_no_descent(0), regul_scheme_ptr(_regul_scheme_ptr), description(_description)
{
	svd_package = new SVD_LAPACK();
}

void SVDSolver::set_svd_package(PestppOptions::SVD_PACK _svd_pack)
{
	if(_svd_pack == PestppOptions::PROPACK){
		delete svd_package;
		svd_package = new SVD_PROPACK;
	}
	else {
		delete svd_package;
		svd_package = new SVD_LAPACK;
	}
	svd_package->set_max_sing(svd_info.maxsing);
	svd_package->set_eign_thres(svd_info.eigthresh);
}


SVDSolver::~SVDSolver(void)
{
	delete svd_package;
}

ModelRun& SVDSolver::solve(RunManagerAbstract &run_manager, TerminationController &termination_ctl, int max_iter, 
	const Parameters &ctl_pars, ModelRun &optimum_run)
{
	ostream &os = file_manager.rec_ofstream();
	cur_solution.set_ctl_parameters(ctl_pars);
	// Start Solution iterations
	bool save_nextjac = false;
	for (int iter_num=1; iter_num<=max_iter;++iter_num) {
		cout << "OPTIMISATION ITERATION NUMBER: " << termination_ctl.get_iteration_number()+1 << endl;
		os   << "OPTIMISATION ITERATION NUMBER: " << termination_ctl.get_iteration_number()+1 << endl;
		cout << "  Iteration type: " << get_description() << endl;
		os   << "  Iteration type: " << get_description() << endl;
		cout << "  SVD Package: " << svd_package->description << endl;
		os   << "  SVD Package: " << svd_package->description << endl;
		os   << "    Model calls so far : " << run_manager.get_total_runs() << endl;
		iteration(run_manager, termination_ctl, false);
		if (save_nextjac) {
			jacobian.save(file_manager.jacobian_filename());
		}
		if (!optimum_run.obs_valid() || cur_solution.get_phi() < optimum_run.get_phi())
		{
			optimum_run = cur_solution;
			// save new optimum parameters to .par file
			optimum_run.get_ctl_pars().save(file_manager.par_filename(), optimum_run.get_par_tran().get_offset_ptr(), 
				optimum_run.get_par_tran().get_scale_ptr());
			jacobian.save(file_manager.jacobian_filename());
			// jacobian calculated next iteration will be at the current parameters and
			// will be more accurate than the one caluculated at the begining of this iteration
			save_nextjac = true;
		}
		os << endl << endl;
		if (termination_ctl.check_last_iteration()){
			break;
		}
	}
	return cur_solution;
}


SVDSolver::Upgrade SVDSolver::calc_upgrade_vec(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
	const LaVectorDouble &Residuals, const vector<string> &par_name_vec, const vector<string> &obs_name_vec)
{
	Upgrade upgrade;
	upgrade.par_name_vec = par_name_vec;
	LaGenMatDouble jac = jacobian.get_matrix(par_name_vec, obs_name_vec);
	LaGenMatDouble SqrtQ_J = Q_sqrt * jac;
	LaVectorDouble Sigma(min(SqrtQ_J.size(0), SqrtQ_J.size(1)));
	LaGenMatDouble U(SqrtQ_J.size(0), SqrtQ_J.size(0));
	LaGenMatDouble Vt(SqrtQ_J.size(1), SqrtQ_J.size(1));
	svd_package->solve_ip(SqrtQ_J, Sigma, U, Vt);
	//calculate the number of singluar values above the threshold
	LaGenMatDouble SqrtQ_J_inv =SVD_inv(U, Sigma, Vt, svd_info.maxsing, svd_info.eigthresh, upgrade.n_sing_val_used);
	   upgrade.svd_uvec =(SqrtQ_J_inv * Q_sqrt) * Residuals;
	upgrade.tot_sing_val = Sigma.size();
	// Calculate svd and greatest descent unit vectors
	upgrade.svd_norm = Blas_Norm2(upgrade.svd_uvec);
	if (upgrade.svd_norm != 0) upgrade.svd_uvec.scale(1.0 /upgrade.svd_norm);
	upgrade.grad_uvec = Q_sqrt.tran_q_mat_mult(jac) * Residuals;
	upgrade.grad_norm = Blas_Norm2(upgrade.grad_uvec);
	if (upgrade.grad_norm != 0) upgrade.grad_uvec.scale(1.0 /upgrade.grad_norm);
	return upgrade;
}

map<string,double> SVDSolver::freeze_parameters(ModelRun &model_run, const Upgrade &upgrade, bool use_descent,
	double scale)
{
	Parameters upgrade_pars = model_run.get_numeric_pars();
	Parameters tmp_parameters;
	int ip = 0;
	map<string,double> tmp_svd;
	map<string,double> tmp_descent;

	// tmp_parameters is used to a store a list of parameters which aare limited and need
	// to be frozen.

	// First get a list of parameters that are at their bounds in the direction computed
	//using SVD
	for(Parameters::iterator b=upgrade_pars.begin(), e=upgrade_pars.end(); b!=e; ++b, ++ip)
	{
		b->second += upgrade.svd_uvec(ip)*upgrade.svd_norm*scale;
	}
	tmp_svd = limit_parameters_ip(model_run.get_numeric_pars(), upgrade_pars);

	// If use_descent==true then only freeze parameters that are at their bounds in both the SVD 
	// and direction of steepest descent.  If descent is true remove the entries in tmp_svd
	// that are not also at their bounds in the direction of steepest descent
	if (use_descent)
	{
		// Test Greatest descent parameter upgrade
		upgrade_pars =  model_run.get_numeric_pars();
		ip = 0;
		for(Parameters::iterator b=upgrade_pars.begin(), e=upgrade_pars.end(); b!=e; ++b, ++ip)
		{
			b->second += upgrade.grad_uvec(ip)*upgrade.svd_norm*scale;
		}
		tmp_descent =  limit_parameters_ip(model_run.get_numeric_pars(), upgrade_pars);
		// remove entries from tmp_svd that are not also in tmp_descent
		map<string,double>::iterator descent_end = tmp_descent.end();
		for (map<string,double>::iterator b=tmp_svd.begin(), e=tmp_svd.end(); b!=e;) {
			if (tmp_descent.find((*b).first) != descent_end ) {
				tmp_svd.erase(b++);
			}
			else {
				++b;
			}
		}
	}
	model_run.freeze_parameters(tmp_svd);
	return tmp_svd;
}


ModelRun SVDSolver::iterative_parameter_freeze(const ModelRun &model_run, 
	Upgrade &upgrade, const QSqrtMatrix &q_sqrt_mat, const LaVectorDouble &residuals_vec, 
	const vector<string> & obs_names_vec, bool use_desent, double scale)
{
	ostream &os = file_manager.rec_ofstream();
	ModelRun base_run(model_run);
	map<string,double> freeze_pars;

	freeze_pars = freeze_parameters(base_run, upgrade, use_desent);
	while (! freeze_pars.empty() ) {
		for (map<string,double>::iterator b=freeze_pars.begin(), e=freeze_pars.end(); b!=e; ++b) {
			os << "  freezing parameter: " << (*b).first <<  " at:" << (*b).second << endl;
		}
		// If there are any frozen parameter, rebuild the LaGenMatDouble jacobian matrix so that it does not include them
		for(map<string, double>::const_iterator b_f=freeze_pars.begin(), e_f=freeze_pars.end();
			b_f!=e_f; ++b_f)
		{
			upgrade.par_name_vec.erase(std::remove(upgrade.par_name_vec.begin(), upgrade.par_name_vec.end(),
				(*b_f).first), upgrade.par_name_vec.end());
		}
		upgrade = calc_upgrade_vec(jacobian, q_sqrt_mat, residuals_vec, upgrade.par_name_vec, obs_names_vec);
		freeze_pars = freeze_parameters(base_run, upgrade, use_desent, scale);
	}
	base_run.freeze_parameters(freeze_pars);
	return base_run;
}


double SVDSolver::add_model_run(RunManagerAbstract &run_manager, const ParamTransformSeq &numeric2model_tran_seq, const Parameters &numeric_base_par, 
	const Upgrade &upgrade, double rot_fac, double scale)
{
	const double PI = 3.141592;
	double tmp_norm;
	double rot_angle = 0.0;
	LaVectorDouble upgrade_vec = ((1.0 - rot_fac) * upgrade.svd_uvec + rot_fac * upgrade.grad_uvec);
	tmp_norm = Blas_Norm2(upgrade_vec);
	// Compute unit upgrade vector;
	if (tmp_norm != 0) {
		upgrade_vec *= (1.0 / tmp_norm);
	}
	rot_angle = acos(min(1.0, upgrade.svd_uvec * upgrade_vec)) * 180.0 / PI;
	upgrade_vec *= upgrade.svd_norm*scale;
	
	// update numeric parameters in upgrade_run and test them
	{
		Parameters upgrade_numeric_pars = numeric_base_par;
		upgrade_numeric_pars.add_upgrade(upgrade.par_name_vec, upgrade_vec);
		// impose limits on parameter upgrade vector
		limit_parameters_ip(numeric_base_par, upgrade_numeric_pars);
		// add run to run manager
		run_manager.add_run(numeric2model_tran_seq.numeric2model_cp(upgrade_numeric_pars));
	}
	return rot_angle;
}

void SVDSolver::iteration(RunManagerAbstract &run_manager, TerminationController &termination_ctl, bool calc_init_obs)
{
	ostream &os = file_manager.rec_ofstream();
	const double PI = 3.141592;
	LaVectorDouble residuals_vec;
	ModelRun base_run(cur_solution);
	vector<string> obs_names_vec;

	// Calculate Jacobian
	if (!cur_solution.obs_valid() || calc_init_obs == true) {
		 calc_init_obs = true;
	 }
	cout << "    calculating jacobian... ";
	jacobian.calculate(cur_solution, cur_solution.get_numeric_pars().get_keys(),  cur_solution.get_obs_template().get_keys(),
			*par_group_info_ptr, *ctl_par_info_ptr,run_manager,  *prior_info_ptr,
			phiredswh_flag, calc_init_obs);
	cout << endl;

	// update regularization weight factor
	double tikhonov_weight = regul_scheme_ptr->get_weight(cur_solution);
	// write out report for starting phi
	obj_func->phi_report(os, cur_solution.get_obs(), cur_solution.get_ctl_pars(), tikhonov_weight);
	// populate vectors with sorted observations (standard and prior info) and parameters
	obs_names_vec = cur_solution.get_obs().get_keys();
	{
		vector<string> prior_info_names = prior_info_ptr->get_keys();
		obs_names_vec.insert(obs_names_vec.end(), prior_info_names.begin(), prior_info_names.end());
	}
	// build weights matrix sqrt(Q)
	QSqrtMatrix Q_sqrt(*obs_info_ptr, obs_names_vec, prior_info_ptr, tikhonov_weight);
	//build residuals vector
	residuals_vec = -1.0 * stlvec2LaVec(cur_solution.get_residuals_vec(obs_names_vec));

	// Freeze Parameters
	bool use_desent;
	if (prev_phi_percent > 99.0 && num_no_descent > 2) {
		use_desent = true;
		num_no_descent = 0;
	}
	else {
		use_desent = false;
		++num_no_descent;
	}

	Upgrade upgrade;
	upgrade = calc_upgrade_vec(jacobian, Q_sqrt, residuals_vec, base_run.get_numeric_pars().get_keys(), obs_names_vec);
	ModelRun upgrade_run = iterative_parameter_freeze(base_run, upgrade, Q_sqrt, residuals_vec, obs_names_vec, use_desent);
	ModelRun best_upgrade_run(base_run);

	os << endl;
	os << "      SVD information:" << endl;
	os << "        number of singular values used: " << upgrade.n_sing_val_used << "/" << upgrade.tot_sing_val << endl;
	os << "        upgrade vector magnitude (without limits or bounds) = " << upgrade.svd_norm << endl;
	os << "        angle to direction of greatest descent: ";
	os <<  acos(min(1.0, upgrade.svd_uvec * upgrade.grad_uvec)) * 180.0 / PI << " deg" << endl;
	os << endl;
	//compute rotation factor and try parameter upgrades
	Parameters::iterator b;
	Parameters::iterator e;

	//Build model runs
	int n_runs = 7;
	double rot_angle;
	vector<double> rot_fac_vec;
	vector<double> rot_angle_vec;
	run_manager.allocate_memory(base_run.get_model_pars(), base_run.get_obs_template(), n_runs);
	 
	//Make two runs with rotation factor of zero  but scaled down
	//ModelRun upgrade_run_scale(base_run);
	//Upgrade upgrade_scale;
	// Add scaled run
	//double rot_angle;
	//upgrade_scale = calc_upgrade_vec(jacobian, Q_sqrt, residuals_vec, base_run.get_numeric_pars().get_keys(), obs_names_vec);
	//upgrade_run_scale = iterative_parameter_freeze(base_run, upgrade_scale, Q_sqrt, residuals_vec, obs_names_vec, use_desent, 1.1);
	//rot_angle = add_model_run(run_manager, upgrade_run_scale.get_par_tran(), base_run.get_numeric_pars(),
	//		upgrade, 0.0, 1.1);
	//rot_angle_vec.push_back(rot_angle);
	//rot_fac_vec.push_back(0.0);

	for(int i=0; i<7; ++i) {

		double tmp_rot_fac[] =  {0.0, 0.01, 0.1, 0.2, 0.5, 0.7, 1.0};
		rot_fac_vec.push_back(tmp_rot_fac[i]);
		rot_angle = add_model_run(run_manager, upgrade_run.get_par_tran(), 
			base_run.get_numeric_pars(),upgrade, rot_fac_vec[i], 1.0);
		rot_angle_vec.push_back(rot_angle);
		rot_fac_vec.push_back(tmp_rot_fac[i]);
	}
	// process model runs
	cout << "    testing upgrade vectors... ";
	run_manager.run();
	cout << endl;
	bool best_run_updated_flag = false;
	for(int i=0; i<n_runs; ++i) {
		double rot_fac = 0;
		if (i>1) rot_fac = rot_fac_vec[i-2];
		ModelRun upgrade_run(base_run);
		run_manager.get_run(upgrade_run, i, RunManagerAbstract::FORCE_PAR_UPDATE);
		streamsize n_prec = os.precision(2);
		os << "      Rotation Factor = ";
		os << setiosflags(ios::fixed)<< setw(4) << rot_fac;
		os << " (" << rot_angle_vec[i] << " deg)";
		os.precision(n_prec);
		os.unsetf(ios_base::floatfield); // reset all flags to default
		os << ";  phi = " << upgrade_run.get_phi(tikhonov_weight); 
		os.precision(2);
		os << setiosflags(ios::fixed);
		os << " ("  << upgrade_run.get_phi(tikhonov_weight)/cur_solution.get_phi(tikhonov_weight)*100 << "% starting phi)" << endl;
		os.precision(n_prec);
		os.unsetf(ios_base::floatfield); // reset all flags to default
		if ( upgrade_run.obs_valid() &&  !best_run_updated_flag || upgrade_run.get_phi() <  best_upgrade_run.get_phi()) {
			best_run_updated_flag = true;
			best_upgrade_run = upgrade_run;
		}
	}
	// clean up run_manager memory
	run_manager.free_memory();

	// set flag to switch to central derivatives next iteration
	if(cur_solution.get_phi() != 0 && !phiredswh_flag && (cur_solution.get_phi()-best_upgrade_run.get_phi())/cur_solution.get_phi() < ctl_info->phiredswh) {
		phiredswh_flag = true;
		os << endl << "      Switching to central derivatives:" << endl;
	}

	cout << "  Starting phi = " << cur_solution.get_phi() << ";  ending phi = " << best_upgrade_run.get_phi() <<
		"  ("  << best_upgrade_run.get_phi()/cur_solution.get_phi()*100 << "% starting phi)" << endl;
	cout << endl;
	os << endl;
	iteration_update_and_report(os, best_upgrade_run, termination_ctl);
	prev_phi_percent =  best_upgrade_run.get_phi()/cur_solution.get_phi()*100;
	cur_solution = best_upgrade_run;
	// Clear Frozen Parameters
	cur_solution.thaw_parameters();
}

map<string, double> SVDSolver::limit_parameters_ip(const Parameters &init_numeric_pars, Parameters &upgrade_numeric_pars)
{
	const string *name;
	const ParameterRec *p_info;
	double p_init;
	double p_upgrade;
	double p_limit;
	double b_facorg_lim;
	map<string, double> bound_parameters;
	pair<bool, double> par_limit;
	Parameters model_parameters_at_limit;
	Parameters::const_iterator pu_iter;
	Parameters::const_iterator pu_end = upgrade_numeric_pars.end();
	Parameters init_ctl_pars = par_transform.numeric2ctl_cp(init_numeric_pars);
	Parameters upgrade_ctl_pars  = par_transform.numeric2ctl_cp(upgrade_numeric_pars);
	Parameters::const_iterator mu_iter;
	Parameters::const_iterator mu_end = upgrade_ctl_pars.end();

	for(Parameters::const_iterator b=init_ctl_pars.begin(), e=init_ctl_pars.end(); b!=e; ++b)
	{
		par_limit = pair<bool, double>(false, 0.0);
		name = &(*b).first;  // parameter name
		p_init = (*b).second; // inital parameter value
		pu_iter = upgrade_ctl_pars.find(*name);
		p_upgrade = (*pu_iter).second;  // upgrade parameter value
		p_info = ctl_par_info_ptr->get_parameter_rec_ptr(*name);

		double init_value = ctl_par_info_ptr->get_parameter_rec_ptr(*name)->init_value;
		if(init_value == 0.0)
		{
			init_value = ctl_par_info_ptr->get_parameter_rec_ptr(*name)->ubnd / 4.0;
		}
		b_facorg_lim = ctl_info->facorig * init_value;
		if (abs(p_init) >= b_facorg_lim) {
			b_facorg_lim = p_init;
		}

		// Check Relative Chanage Limit
		if(p_info->chglim == "RELATIVE" && abs((p_upgrade - p_init) / b_facorg_lim) > ctl_info->relparmax)
		{
			par_limit.first = true;
			par_limit.second = p_init + sign(p_upgrade - p_init) * ctl_info->relparmax *  abs(b_facorg_lim);
		}

		// Check Factor Change Limit
		else if(p_info->chglim == "FACTOR") {
			if (b_facorg_lim > 0 && p_upgrade < b_facorg_lim/ctl_info->facparmax ) 
			{
				par_limit.first = true;
				par_limit.second = b_facorg_lim / ctl_info->facparmax;
			}
			else if (b_facorg_lim > 0 && p_upgrade > b_facorg_lim*ctl_info->facparmax ) 
			{
				par_limit.first = true;
				par_limit.second = b_facorg_lim * ctl_info->facparmax;
			}
			else if (b_facorg_lim < 0 && p_upgrade < b_facorg_lim*ctl_info->facparmax)
			{
				par_limit.first = true;
				par_limit.second = b_facorg_lim * ctl_info->facparmax;
			}
			else if (b_facorg_lim < 0 && p_upgrade > b_facorg_lim/ctl_info->facparmax)
			{
				par_limit.first = true;
				par_limit.second = b_facorg_lim / ctl_info->facparmax;
			}
		}
		// Check parameter upper bound
		if((!par_limit.first && p_upgrade > p_info->ubnd) || 
			(par_limit.first && par_limit.second > p_info->ubnd)) {
				par_limit.first = true;
				par_limit.second = p_info->ubnd;
				bound_parameters[*name] = p_info->ubnd;
		}
		// Check parameter lower bound
		else if((!par_limit.first && p_upgrade < p_info->lbnd) || 
			(par_limit.first && par_limit.second < p_info->lbnd)) {
				par_limit.first = true;
				par_limit.second = p_info->lbnd;
				bound_parameters[*name] = p_info->lbnd;
		}
		// Add any limited parameters to model_parameters_at_limit
		if (par_limit.first) {
			model_parameters_at_limit.insert(*name, par_limit.second);
		}
	}
	// Calculate most stringent limit factor on a PEST parameter
	double limit_factor= 1.0;
	double tmp_limit;
	string limit_parameter_name = "";
	Parameters pest_parameters_at_limit = par_transform.model2numeric_cp(model_parameters_at_limit);
	for(Parameters::const_iterator b=pest_parameters_at_limit.begin(), e=pest_parameters_at_limit.end();
		b!=e; ++b) {
			name = &(*b).first;
			p_limit = (*b).second;
			p_init = init_numeric_pars.get_rec(*name);
			p_upgrade = upgrade_numeric_pars.get_rec(*name);
			tmp_limit = (p_limit - p_init) / (p_upgrade - p_init);
			if (tmp_limit < limit_factor)  {
				limit_factor = (p_limit - p_init) / (p_upgrade - p_init);
				limit_parameter_name = *name;
			}
	}
	// Apply limit factor to PEST upgrade parameters
	for(Parameters::iterator b=upgrade_numeric_pars.begin(), e=upgrade_numeric_pars.end();
		b!=e; ++b) {
		name = &(*b).first;
		p_init = init_numeric_pars.get_rec(*name);
		(*b).second = p_init + ((*b).second - p_init) *  limit_factor;
	}
	map<string,double>::iterator iter;
	iter = bound_parameters.find(limit_parameter_name);
	if (iter != bound_parameters.end()) {
		double value = bound_parameters[limit_parameter_name];
		bound_parameters.clear();
		bound_parameters[limit_parameter_name] = value;
	}
	else {
		bound_parameters.clear();
	}
	return bound_parameters;
}

void SVDSolver::param_change_stats(double p_old, double p_new, bool &have_fac, double &fac_change, bool &have_rel, double &rel_change) 
{
	have_rel = have_fac = true;
	double a = max(abs(p_new), abs(p_old));
		double b = min(abs(p_new), abs(p_old));
		// compute relative change
		if (p_old == 0) {
			have_rel = false;
			rel_change = -9999;
		}
		else 
		{
			rel_change = (p_old - p_new) / p_old;
		}
		//compute factor change
		if (p_old == 0.0 || p_new == 0.0) {
			have_fac = false;
			fac_change = -9999;
		}
		else {
			fac_change = a / b;
		}
	}


void SVDSolver::iteration_update_and_report(ostream &os, ModelRun &upgrade, TerminationController &termination_ctl)
{
	const string *p_name;
	double p_old, p_new;
	double fac_change=-9999, rel_change=-9999;
	bool have_fac=false, have_rel=false;
	const Parameters &old_ctl_pars = par_transform.numeric2ctl_cp(cur_solution.get_numeric_pars());
	const Parameters &new_ctl_pars = par_transform.numeric2ctl_cp(upgrade.get_numeric_pars());

	os << "    Parameter Upgrades (Control File Parameters)" << endl;
	os << "     Parameter      Current       Previous       Factor       Relative" << endl;
	os << "        Name         Value         Value         Change        Change" << endl;
	os << "    ------------  ------------  ------------  ------------  ------------" << endl;

	for( Parameters::const_iterator nb=new_ctl_pars.begin(), ne=new_ctl_pars.end();
		nb!=ne; ++nb) {
		p_name = &((*nb).first);
		p_new =  (*nb).second;
		p_old = old_ctl_pars.get_rec(*p_name);
		param_change_stats(p_old, p_new, have_fac, fac_change, have_rel, rel_change);
		os << left;
		os << "    " << setw(12) << *p_name;
		os << right;
		os << "  " << setw(12) << p_new;
		os << "  " << setw(12) << p_old;
		if (have_fac)
			os << "  " << setw(12) << fac_change;
		else
			os << "  " << setw(12) << "N/A";
		if (have_rel)
			os << "  " << setw(12) << rel_change;
		else
			os << "  " << setw(12) << "N/A";
		os << endl;
	}
	os << endl;

	double max_fac_change = 0;
	double max_rel_change = 0;
	const string *max_fac_par = 0;
	const string *max_rel_par = 0;
	os << "    Parameter Upgrades (Transformed Numeric Parameters)" << endl;
	os << "     Parameter      Current       Previous       Factor       Relative" << endl;
	os << "        Name         Value         Value         Change        Change" << endl;
	os << "    ------------  ------------  ------------  ------------  ------------" << endl;

	for( Parameters::const_iterator nb=upgrade.get_numeric_pars().begin(), ne=upgrade.get_numeric_pars().end();
		nb!=ne; ++nb) {
		p_name = &((*nb).first);
		p_new =  (*nb).second;
		p_old = cur_solution.get_numeric_pars().get_rec(*p_name);
		param_change_stats(p_old, p_new, have_fac, fac_change, have_rel, rel_change);
		if (fac_change >= max_fac_change) 
		{
			max_fac_change = fac_change;
			max_fac_par = p_name;
		}
		if (abs(rel_change) >= abs(max_rel_change))
		{
			max_rel_change = rel_change;
			max_rel_par = p_name;
		}
		os << left;
		os << "    " << setw(12) << *p_name;
		os << right;
		os << "  " << setw(12) << p_new;
		os << "  " << setw(12) << p_old;
		if (have_fac)
			os << "  " << setw(12) << fac_change;
		else
			os << "  " << setw(12) << "N/A";
		if (have_rel)
			os << "  " << setw(12) << rel_change;
		else
			os << "  " << setw(12) << "N/A";
		os << endl;
	}
	os << endl;
	os << "   Maximum changes in transformed numeric parameters:" << endl;
	os << "     Maximum relative change = " << max_rel_change << "   [" << *max_rel_par << "]" << endl;
	os << "     Maximum factor change = " << max_fac_change << "   [" << *max_fac_par << "]" << endl;
	termination_ctl.process_iteration(upgrade.get_phi(), max_rel_change);
}