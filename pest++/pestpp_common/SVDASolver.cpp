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

SVDASolver::SVDASolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_base_parameter_group_info_ptr, 
	const ParameterInfo *_ctl_par_info_ptr, const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
	const ParamTransformSeq &_par_transform, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, const Regularization *_regul_scheme, int _max_freeze_iter,
	OutputFileWriter &_output_file_writer, RestartController &_restart_controller, SVDSolver::MAT_INV _mat_inv)
	: SVDSolver(_ctl_info, _svd_info, _base_parameter_group_info_ptr, _ctl_par_info_ptr, _obs_info,
		_file_manager, _observations, _obj_func, _par_transform, _prior_info_ptr, _jacobian, 
		_regul_scheme, _max_freeze_iter, _output_file_writer, _restart_controller, _mat_inv, "super parameter solution")
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
	double limit_fac = 1.0;
	for(auto &ipar : init_numeric_pars)
	{
		par_limit = pair<bool, double>(false, 0.0);
		name = &(ipar.first);  // parameter name
		val_init = ipar.second; // inital parameter value
		auto iter_pu = upgrade_numeric_pars.find(*name);
		assert (iter_pu != upgrade_numeric_pars.end());
		val_upgrade = (*iter_pu).second;  // upgrade parameter value
		// only use relative change limits
		limit = abs((val_upgrade - val_init) / (val_init * ctl_info->relparmax));
		if (limit < 1)
		{
			limit_fac = min(limit_fac, 1.0 / limit);
		}
	}
	for(auto &ipar : upgrade_numeric_pars)
	{
		ipar.second *= limit_fac;
	}
			
	//convert parameters to their deriviative form and check any any that have exceeded their bounds
	Parameters derivative_update_pars = par_transform.numeric2derivative_cp(upgrade_numeric_pars);
	Parameters freeze_derivative_par;
	for (auto &ipar :  derivative_update_pars)
	{
		const string &name = ipar.first;
		const ParameterRec *p_info = ctl_par_info_ptr->get_parameter_rec_ptr(name);
		if(ipar.second > p_info->ubnd)
		{
			freeze_derivative_par[name] = p_info->ubnd;
		}
		// Check parameter lower bound
		else if(ipar.second < p_info->lbnd)
		{
			freeze_derivative_par[name] = p_info->lbnd;
		}
	}
	return freeze_derivative_par;
}

void SVDASolver::calc_upgrade_vec(double i_lambda, vector<Parameters> &frozen_derivative_par_vec, QSqrtMatrix &Q_sqrt, VectorXd &residuals_vec,
	vector<string> &numeric_par_names_vec, vector<string> &obs_names_vec, const Parameters &base_run_numeric_pars,
	Upgrade &ml_upgrade, Parameters &new_pars, MarquardtMatrix marquardt_type)
{
	////dew_clean_up
	//ml_upgrade = calc_lambda_upgrade_vec(jacobian, Q_sqrt, residuals_vec, numeric_par_names_vec, obs_names_vec,
	//	base_run_numeric_pars, Parameters(), i_lambda, marquardt_type);
	////Start out with new_pars as numeric parameters
	//new_pars = apply_upgrade(base_run_numeric_pars, ml_upgrade, 1.0);
	//Parameters new_frozen_pars = limit_parameters_ip(base_run_numeric_pars, new_pars);
	//frozen_derivative_par_vec.push_back(new_frozen_pars);
	////transform new_par to derivative parameters
	//par_transform.numeric2derivative_ip(new_pars);
	//// impose frozen parameters
	//for (auto &i : frozen_derivative_par_vec.back())
	//{
	//	new_pars[i.first] = i.second;
	//}

}
void SVDASolver::iteration(RunManagerAbstract &run_manager, TerminationController &termination_ctl, bool calc_init_obs)
{
	ostream &fout_restart = file_manager.get_ofstream("rst");
	fout_restart << "super_par_iteration" << endl;
	ostream &os = file_manager.rec_ofstream();
	ModelRun base_run(cur_solution);
	vector<string> obs_names_vec = base_run.get_obs_template().get_keys(); 
	vector<string> numeric_par_names_vec;


	// Calculate Jacobian
	cout << "  calculating jacobian... ";

	Parameters base_ctl_pars = base_run.get_ctl_pars();
	// make sure these are all in bounds
	for (auto &ipar :  base_ctl_pars)
	{
		const string &name = ipar.first;
		const ParameterRec *p_info = ctl_par_info_ptr->get_parameter_rec_ptr(name);
		if(ipar.second > p_info->ubnd)
		{
			ipar.second = p_info->ubnd;
		}
		// Check parameter lower bound
		else if(ipar.second < p_info->lbnd)
		{
			ipar.second = p_info->lbnd;
		}
	}

	while (true)
	{
		// fix frozen parameters in SVDA transformation
		par_transform.get_svda_ptr()->update_add_frozen_pars(base_run.get_frozen_ctl_pars());
		par_transform.get_svda_fixed_ptr()->reset(par_transform.get_svda_ptr()->get_frozen_derivative_pars());
		// need to reset parameters and the numeric parameters changed when the SVDA transformation was changed above
		//DEW_change
		base_run.set_ctl_parameters(base_run.get_ctl_pars());
		Parameters numeric_pars = par_transform.ctl2numeric_cp(base_run.get_ctl_pars());
		numeric_par_names_vec = numeric_pars.get_keys();
		Parameters derivative_pars = par_transform.ctl2derivative_cp(base_run.get_ctl_pars());
		set<string> out_of_bound_pars;
		if (!base_run.obs_valid() || calc_init_obs == true) {
		 calc_init_obs = true;
		}
		super_parameter_group_info = par_transform.get_svda_ptr()->build_par_group_info(*par_group_info_ptr);
		bool success = jacobian.build_runs(base_run, numeric_par_names_vec, par_transform,
			super_parameter_group_info, *ctl_par_info_ptr, run_manager, out_of_bound_pars,
			phiredswh_flag, calc_init_obs);

		jacobian.make_runs(run_manager);

		bool success2 = jacobian.process_runs(numeric_par_names_vec, par_transform,
			super_parameter_group_info, *ctl_par_info_ptr, run_manager, *prior_info_ptr, out_of_bound_pars,
			phiredswh_flag, calc_init_obs);

		success = success && success2;

		if (success)
		{
			break;
		}
		else
		{
			//add newly frozen parameters to list
			Parameters new_frz_derivative_pars;
			for (auto &ipar : out_of_bound_pars)
			{
				const auto iter = base_ctl_pars.find(ipar);
				assert(iter != base_ctl_pars.end());
				new_frz_derivative_pars.insert(iter->first, iter->second);
			}
			if (new_frz_derivative_pars.size() > 0)
			{
				base_run.add_frozen_ctl_parameters(new_frz_derivative_pars);
			}
		}
	}
	cout << endl;
	cout << "  computing upgrade vectors... " << endl;

	// update regularization weight factor
	double tikhonov_weight = regul_scheme_ptr->get_weight(base_run);
	// write out report for starting phi
	obj_func->phi_report(os, base_run.get_obs(), base_run.get_ctl_pars(), tikhonov_weight);
	// populate vectors with sorted observations (standard and prior info) and parameters
	{
		vector<string> prior_info_names = prior_info_ptr->get_keys();
		obs_names_vec.insert(obs_names_vec.end(), prior_info_names.begin(), prior_info_names.end());
	}

	//build residuals vector
	VectorXd residuals_vec = -1.0 * stlvec_2_egienvec(base_run.get_residuals_vec(obs_names_vec));

	//Build model runs
	run_manager.reinitialize(file_manager.build_filename("rnu"));
	vector<double> magnitude_vec;

	//Marquardt Lambda Update Vector
	Upgrade ml_upgrade;
	const Parameters &base_run_derivative_pars =  par_transform.ctl2derivative_cp(base_run.get_ctl_pars());
	Parameters base_numeric_pars = par_transform.ctl2numeric_cp(base_run.get_ctl_pars());
	vector<Parameters> frozen_derivative_par_vec;
	Parameters *new_frozen_par_ptr = 0;

	// build weights matrix sqrt(Q)
	QSqrtMatrix Q_sqrt(obs_info_ptr, prior_info_ptr, tikhonov_weight);

	double tmp_lambda[] = {1.0e-5, 0.1, 1.0, 10.0, 100.0, 1000.0};
	vector<double> lambda_vec(tmp_lambda, tmp_lambda+sizeof(tmp_lambda)/sizeof(double));
	lambda_vec.push_back(best_lambda);
	lambda_vec.push_back(best_lambda / 2.0);
	lambda_vec.push_back(best_lambda * 2.0);
	std::sort(lambda_vec.begin(), lambda_vec.end());
	auto iter = std::unique(lambda_vec.begin(), lambda_vec.end());
	lambda_vec.resize(std::distance(lambda_vec.begin(), iter));
	stringstream message;
	int i_update_vec = 0;
	for (double i_lambda : lambda_vec)
	{
		std::cout << string(message.str().size(), '\b');
		message.str("");
		message << "  computing upgrade vector (lambda = " << i_lambda << ")  " << ++i_update_vec << " / " << lambda_vec.size() << "             ";
		std::cout << message.str();

		Parameters new_pars;

		calc_upgrade_vec(i_lambda, frozen_derivative_par_vec, Q_sqrt, residuals_vec,
			numeric_par_names_vec, obs_names_vec, base_numeric_pars,
			ml_upgrade, new_pars, MarquardtMatrix::IDENT);

		//transform new_pars to model parameters
		par_transform.derivative2model_ip(new_pars);
		run_manager.add_run(new_pars, "IDEN", i_lambda);
		magnitude_vec.push_back(ml_upgrade.norm);
	}

	for (double i_lambda : lambda_vec)
	{
		std::cout << string(message.str().size(), '\b');
		message.str("");
		message << "  computing upgrade vector (lambda = " << i_lambda << ")  " << ++i_update_vec << " / " << lambda_vec.size() << "             ";
		std::cout << message.str();

		Parameters new_pars;

		calc_upgrade_vec(i_lambda, frozen_derivative_par_vec, Q_sqrt, residuals_vec,
			numeric_par_names_vec, obs_names_vec, base_numeric_pars,
			ml_upgrade, new_pars, MarquardtMatrix::JTQJ);

		//transform new_pars to model parameters
		par_transform.derivative2model_ip(new_pars);
		run_manager.add_run(new_pars, "DIAG", i_lambda);
		magnitude_vec.push_back(ml_upgrade.norm);
	}

	cout << endl;

	os << endl;

	// process model runs
	cout << "  testing upgrade vectors... ";
	fout_restart << "upgrade_model_runs_built " << run_manager.get_cur_groupid() << endl;
	run_manager.run();
	cout << endl;
	bool best_run_updated_flag = false;
	ModelRun best_upgrade_run(base_run);

	long jac_num_nonzero = jacobian.get_nonzero();
	long jac_num_total = jacobian.get_size();
	long jac_num_zero = jac_num_total - jac_num_nonzero;
	streamsize n_prec = os.precision(2);
	os << "    Number of terms in the jacobian equal to zero: " << jac_num_zero << " / " << jac_num_total
		<< " (" << double(jac_num_zero) / double(jac_num_total) * 100 << "%)" << endl << endl;
	os.precision(n_prec);

	os << "    Summary of upgrade runs:" << endl;
	for(int i=0; i<run_manager.get_nruns(); ++i) {
		ModelRun upgrade_run(base_run);
		Parameters tmp_pars;
		Observations tmp_obs;
		string lambda_type;
		double i_lambda;
		bool success = run_manager.get_run(i, tmp_pars, tmp_obs, lambda_type, i_lambda);
		if (success)
		{
			par_transform.model2ctl_ip(tmp_pars);
			upgrade_run.update_ctl(tmp_pars, tmp_obs);
			streamsize n_prec = os.precision(2);
			os << "      Lambda = ";
			os << setiosflags(ios::fixed) << setw(8) << i_lambda;
			os << "; Type: " << setw(4) << lambda_type;
			os << ";   length = " << magnitude_vec[i];
			os.precision(n_prec);
			os.unsetf(ios_base::floatfield); // reset all flags to default
			os << ";   phi = " << upgrade_run.get_phi(tikhonov_weight); 
			os.precision(2);
			os << setiosflags(ios::fixed);
			os << " ("  << upgrade_run.get_phi(tikhonov_weight)/base_run.get_phi(tikhonov_weight)*100 << "%)" << endl;
			os.precision(n_prec);
			os.unsetf(ios_base::floatfield); // reset all flags to default
			if ( upgrade_run.obs_valid() &&  !best_run_updated_flag || upgrade_run.get_phi() <  best_upgrade_run.get_phi()) {
				best_run_updated_flag = true;
				best_upgrade_run = upgrade_run;
				new_frozen_par_ptr = &frozen_derivative_par_vec[i];
				best_lambda = i_lambda;
			}
		}
		else
		{
			streamsize n_prec = os.precision(2);
			os << "     Marquardt Lambda = ";
			os << setiosflags(ios::fixed) << setw(8) << i_lambda;
			os << ";   length = " << magnitude_vec[i];
			os.precision(n_prec);
			os.unsetf(ios_base::floatfield); // reset all flags to default
			os << ";    run failed" << endl;
		}
	}
	// Print frozen parameter information for parameters frozen in SVD transformation
	const Parameters &frz_ctl_pars_svd = best_upgrade_run.get_frozen_ctl_pars();
	if (frz_ctl_pars_svd.size() > 0)
	{
		vector<string> keys = frz_ctl_pars_svd.get_keys();
		std::sort(keys.begin(), keys.end());
		os << endl;
		os << "    Parameters previously frozen in SVD transformation:" << endl;
		for (auto &ikey : keys)
		{
			auto iter = frz_ctl_pars_svd.find(ikey);
			if (iter != frz_ctl_pars_svd.end())
			{
				os << "      " << iter->first << " frozen at " << iter->second << endl;
			}
		}
	}

	if (new_frozen_par_ptr!= 0 && new_frozen_par_ptr->size() > 0)
	{
		vector<string> keys = new_frozen_par_ptr->get_keys();
		std::sort(keys.begin(), keys.end());
		os << endl;
		os << "    Parameters frozen during upgrades:" << endl;
		for (auto &ikey : keys)
		{
			auto iter = new_frozen_par_ptr->find(ikey);
			if (iter != new_frozen_par_ptr->end())
			{
				os << "      " << iter->first << " frozen at " << iter->second << endl;
			}
		}
	}
	if (new_frozen_par_ptr!=0) best_upgrade_run.add_frozen_ctl_parameters(*new_frozen_par_ptr);
	// clean up run_manager memory
	run_manager.free_memory();

	// reload best parameters and set flag to switch to central derivatives next iteration
	if(base_run.get_phi() != 0 && !phiredswh_flag &&
		(base_run.get_phi()-best_upgrade_run.get_phi())/base_run.get_phi() < ctl_info->phiredswh)
	{
		phiredswh_flag = true;
		os << endl << "      Switching to central derivatives:" << endl;
	}

	cout << "  Starting phi = " << base_run.get_phi() << ";  ending phi = " << best_upgrade_run.get_phi() <<
		"  ("  << best_upgrade_run.get_phi()/base_run.get_phi()*100 << "%)" << endl;
	cout << endl;
	os << endl;
	iteration_update_and_report(os, best_upgrade_run, termination_ctl);
	prev_phi_percent =  best_upgrade_run.get_phi()/base_run.get_phi()*100;
	cur_solution = best_upgrade_run;
}



SVDASolver::~SVDASolver(void)
{
}
