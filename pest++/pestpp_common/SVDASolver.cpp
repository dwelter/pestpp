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
#include <iomanip>
#include <map>
#include "SVDASolver.h"
#include "ModelRunPP.h"
#include "QSqrtMatrix.h"
#include "eigen_tools.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "PriorInformation.h"

using namespace std;
using namespace pest_utils;
using namespace Eigen;

SVDASolver::SVDASolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_par_group_info_ptr, const ParameterInfo *_ctl_par_info_ptr,
		const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, const Regularization *_regul_scheme, int _max_freeze_iter)
		: SVDSolver(_ctl_info, _svd_info, _par_group_info_ptr, _ctl_par_info_ptr, _obs_info, 
		_file_manager, _observations, _obj_func, _par_transform, _prior_info_ptr, _jacobian, 
		_regul_scheme, _max_freeze_iter, "super parameter solution")
{
}


Parameters SVDASolver::limit_parameters_ip(const Parameters &init_numeric_pars, Parameters &upgrade_numeric_pars)
{
	const string *name;
	double val_init;
	double val_upgrade;
	double limit;

	pair<bool, double> par_limit;
	// superparameters are always limited by RELPARMAX
	// if any base parameters go outof bounds, they should be frozen,
	// a new SVD factorization should be computed, and then this routine should be called again

	//Assume Super parameters are independent and apply limits one at a time
	for(auto &ipar : init_numeric_pars)
	{
		par_limit = pair<bool, double>(false, 0.0);
		name = &(ipar.first);  // parameter name
		val_init = ipar.second; // inital parameter value
		auto iter_pu = upgrade_numeric_pars.find(*name);
		assert (iter_pu != upgrade_numeric_pars.end());
		val_upgrade = (*iter_pu).second;  // upgrade parameter value
		// only use relative change limits
		limit = abs((val_upgrade - val_init) / val_init);
		if(limit > ctl_info->relparmax)
		{
			(*iter_pu).second = val_init + sign(val_upgrade - val_init) * ctl_info->relparmax *  abs(val_init);
		}
	}
	//convert parameters to their deriviative form and check any any that have exceeded their bounds
	Parameters derivative_update_pars = par_transform.numeric2derivative_cp(upgrade_numeric_pars);
	//remove previously frozen parameters
	const Parameters &frz_pars = par_transform.get_svda_ptr()->get_frozen_derivative_pars();
	derivative_update_pars.erase(frz_pars);

	Parameters frozen_derivative_pars;
	const ParameterRec *p_info;
	for (auto &ipar: derivative_update_pars)
	{
		name = &(ipar.first);
		p_info = ctl_par_info_ptr->get_parameter_rec_ptr(*name);
		if(ipar.second > p_info->ubnd)
		{
			frozen_derivative_pars.insert(*name, p_info->ubnd);
		}
		// Check parameter lower bound
		else if(ipar.second < p_info->lbnd)
		{
			frozen_derivative_pars.insert(*name, p_info->lbnd);
		}
	}
	return frozen_derivative_pars;
}

//SVDASolver::Upgrade SVDSolver::calc_lambda_upgrade_vec(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
//	const Eigen::VectorXd &Residuals, const vector<string> &par_name_vec, const vector<string> &obs_name_vec,
//	const Parameters &base_numeric_pars, const Parameters &freeze_ctl_pars, int &tot_sing_val, 
//	double lambda, MarquardtMatrix marquardt_type)
//{
//	Upgrade upgrade;
//	upgrade.frozen_numeric_pars = freeze_numeric_pars;
//	upgrade.par_name_vec = par_name_vec;
//
//	//Compute effect of frozen parameters on the residuals vector
//	Parameters delta_freeze_pars = freeze_ctl_pars;
//	delta_freeze_pars -= base_ctl_pars;
//	base_run.get_par_tran().ctl2numeric_ip(delta_freeze_pars);
//	VectorXd del_residuals = calc_residual_corrections(jacobian, delta_freeze_pars, obs_name_vec);
//	MatrixXd jac = jacobian.get_matrix(obs_name_vec, par_name_vec);
//	MatrixXd SqrtQ_J = Q_sqrt * jac;
//
//	VectorXd Sigma;
//	MatrixXd U;
//	MatrixXd Vt;
//	svd_package->solve_ip(SqrtQ_J, Sigma, U, Vt);
//	if (marquardt_type == MarquardtMatrix::IDENT)
//	{
//		Sigma = Sigma.array() + lambda;
//	}
//	else
//	{
//		Sigma = Sigma.array() + (Sigma.cwiseProduct(Sigma).array() + lambda * lambda).sqrt();
//	}
//	
//	//calculate the number of singluar values above the threshold
//	int max_sing = (par_name_vec.size() < svd_info.maxsing) ? par_name_vec.size() : svd_info.maxsing;
//	MatrixXd SqrtQ_J_inv =SVD_inv(U, Sigma, Vt, max_sing, svd_info.eigthresh, max_sing);
//	VectorXd tmp_svd_uvec =(SqrtQ_J_inv * Q_sqrt) * (Residuals + del_residuals);
//
//	//build map of parameter names to index in original par_name_vec vector
//	unordered_map<string, int> par_vec_name_to_idx;
//	for(int i=0; i<par_name_vec.size(); ++i)
//	{
//		par_vec_name_to_idx[par_name_vec[i]] = i;
//	}
//
//	upgrade.uvec.resize(par_name_vec.size());
//	//tranfere newly computed componets of the ugrade vector to upgrade.svd_uvec
//	for(int i=0; i<par_name_vec.size(); ++i)
//	{
//		auto &it = par_vec_name_to_idx.find(par_name_vec[i]);
//		assert(it != par_vec_name_to_idx.end());
//		upgrade.uvec(it->second) = tmp_svd_uvec(i);
//	}
//	tot_sing_val = Sigma.size();
//	// Calculate svd unit vector
//	upgrade.norm = upgrade.uvec.norm();
//	if (upgrade.norm != 0) upgrade.uvec *= 1.0 /upgrade.norm;
//	return upgrade;
//}





void SVDASolver::iteration(RunManagerAbstract &run_manager, TerminationController &termination_ctl, bool calc_init_obs)
{
	ostream &os = file_manager.rec_ofstream();
	const double PI = 3.141592;
	VectorXd residuals_vec;
	ModelRun base_run(cur_solution);
	vector<string> obs_names_vec;
	vector<string> numeric_par_names_vec = cur_solution.get_numeric_pars().get_keys();

	// fix frozen parameters in SVDA transformation
	base_run.get_par_tran().get_svda_ptr()->update_add_frozen_pars(cur_solution.get_frozen_ctl_pars());

	// Calculate Jacobian
	if (!cur_solution.obs_valid() || calc_init_obs == true) {
		 calc_init_obs = true;
	 }
	cout << "  calculating jacobian... ";
	jacobian.calculate(cur_solution, numeric_par_names_vec,  cur_solution.get_obs_template().get_keys(),
			*par_group_info_ptr, *ctl_par_info_ptr,run_manager,  *prior_info_ptr,
			phiredswh_flag, calc_init_obs);
	cout << endl;
	cout << "  computing upgrade vectors... " << endl;

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
	residuals_vec = -1.0 * stlvec_2_egienvec(cur_solution.get_residuals_vec(obs_names_vec));

	//Build model runs
	run_manager.reinitialize(file_manager.build_filename("rnu"));
	vector<double> rot_fac_vec;
	vector<double> rot_angle_vec;
	vector<double> magnitude_vec;

	//Marquardt Lambda Update Vector
	Parameters new_numeric_pars;
	Upgrade ml_upgrade;
	int tot_sing_val;
	const Parameters &base_run_derivative_pars =  base_run.get_par_tran().ctl2derivative_cp(base_run.get_ctl_pars());
	Parameters base_numeric_pars = base_run.get_numeric_pars();
	Parameters frozen_derivative_pars;
	vector<string> par_name_vec = base_numeric_pars.get_keys();
	vector<Parameters> numeric_par_vec;

	double tmp_lambda[] = {0.1, 1.0, 10.0, 100.0, 1000.0};
	vector<double> lambda_vec(tmp_lambda, tmp_lambda+sizeof(tmp_lambda)/sizeof(double));
	lambda_vec.push_back(best_lambda);
	lambda_vec.push_back(best_lambda / 2.0);
	lambda_vec.push_back(best_lambda * 2.0);
	std::sort(lambda_vec.begin(), lambda_vec.end());
	auto iter = std::unique(lambda_vec.begin(), lambda_vec.end());
	lambda_vec.resize(std::distance(lambda_vec.begin(), iter));
	int max_freeze_iter = 1;
	for (double i_lambda : lambda_vec)
	{
		ml_upgrade = calc_lambda_upgrade_vec(jacobian, Q_sqrt, residuals_vec, par_name_vec, obs_names_vec,
			base_run.get_numeric_pars(), Parameters(), tot_sing_val, i_lambda);
		new_numeric_pars = apply_upgrade(base_numeric_pars, ml_upgrade, 1.0);
		Parameters new_frozen_pars = limit_parameters_ip(base_numeric_pars, new_numeric_pars);
		frozen_derivative_pars.insert(new_frozen_pars);
		numeric_par_vec.push_back(new_numeric_pars);
	}
	for (auto ipar : numeric_par_vec)
	{
		base_run.get_par_tran().numeric2derivative_ip(ipar);
		// impose frozen parameters
		for (auto &i : frozen_derivative_pars)
		{
		  ipar[i.first] = i.second;
		}
		base_run.get_par_tran().derivative2model_ip(ipar);
		run_manager.add_run(ipar);
		rot_fac_vec.push_back(0);
		rot_angle_vec.push_back(0);
		magnitude_vec.push_back(ml_upgrade.norm);
	}

	os << endl;
	os << "      SVD information:" << endl;
	os << endl;

	// process model runs
	cout << "  testing upgrade vectors... ";
	ofstream &fout_restart = file_manager.get_ofstream("rst");
	fout_restart << "upgrade_model_runs_begin_group_id " << run_manager.get_cur_groupid() << endl;
	run_manager.run();
	fout_restart << "upgrade_model_runs_end_group_id " << run_manager.get_cur_groupid() << endl;
	cout << endl;
	bool best_run_updated_flag = false;
	ModelRun best_upgrade_run(base_run);

	for(int i=0; i<run_manager.get_nruns(); ++i) {
		double rot_fac = rot_fac_vec[i];
		ModelRun upgrade_run(base_run);
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(i, tmp_pars, tmp_obs);
		if (success)
		{
			upgrade_run.update(tmp_pars, tmp_obs, ModelRun::FORCE_PAR_UPDATE);
			streamsize n_prec = os.precision(2);
			os << "      Marquardt Lambda = ";
			os << setiosflags(ios::fixed)<< setw(4) << lambda_vec[i];
			//os << " (" << rot_angle_vec[i] << " deg)" << "; length = " << magnitude_vec[i];
			os << "; length = " << magnitude_vec[i];
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
				best_upgrade_run.add_frozen_ctl_parameters(frozen_derivative_pars);
				best_lambda = lambda_vec[i];
			}
		}
		else
		{
			streamsize n_prec = os.precision(2);
			os << "      Rotation Factor = ";
			os << setiosflags(ios::fixed)<< setw(4) << rot_fac;
			os << " (" << rot_angle_vec[i] << " deg)";
			os.precision(n_prec);
			os.unsetf(ios_base::floatfield); // reset all flags to default
			os << ";  run failed" << endl;
		}
	}
	// Print frozen parameter information
	const Parameters &frz_ctl_pars = frozen_derivative_pars;
	if (frz_ctl_pars.size() > 0)
	{
		vector<string> keys = frz_ctl_pars.get_keys();
		std::sort(keys.begin(), keys.end());
		os << endl;
		os << "    Parameters frozen during best upgrade:" << endl;
		for (auto &ikey : keys)
		{
			auto iter = frz_ctl_pars.find(ikey);
			if (iter != frz_ctl_pars.end())
			{
				os << "      " << iter->first << " frozen at " << iter->second << endl;
			}
		}
	}

	// clean up run_manager memory
	run_manager.free_memory();

	// reload best parameters and set flag to switch to central derivatives next iteration
	if(cur_solution.get_phi() != 0 && !phiredswh_flag &&
		(cur_solution.get_phi()-best_upgrade_run.get_phi())/cur_solution.get_phi() < ctl_info->phiredswh)
	{
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
}



SVDASolver::~SVDASolver(void)
{
}
