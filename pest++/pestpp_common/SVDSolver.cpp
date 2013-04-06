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
#include <algorithm>
#include "SVDSolver.h"
#include "RunManagerAbstract.h"
#include "QSqrtMatrix.h"
#include "eigen_tools.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "PriorInformation.h"
#include "Regularization.h"
#include "SVD_PROPACK.h"
#include "OutputFileWriter.h"
#include <sstream>

using namespace std;
using namespace pest_utils;
using namespace Eigen;


void SVDSolver::Upgrade::print(ostream &os) 
{
	//os << "Upgrade: " << n_sing_val_used << "out of " <<  tot_sing_val << "singular vales used" << endl;
	//int vec_size = par_name_vec.size();
	//for (int i=0; i< vec_size; ++i) {
	//	os << setw(21) << par_name_vec[i] << setw(20) << svd_uvec(i) << setw(20) <<grad_uvec(i) << endl;
	//}
}

SVDSolver::SVDSolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_par_group_info_ptr, const ParameterInfo *_ctl_par_info_ptr,
		const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, 
		const Regularization *_regul_scheme_ptr, int _n_rotation_fac, const string &_description)
		: ctl_info(_ctl_info), svd_info(_svd_info), par_group_info_ptr(_par_group_info_ptr), ctl_par_info_ptr(_ctl_par_info_ptr), obs_info_ptr(_obs_info), obj_func(_obj_func),
		  file_manager(_file_manager), observations_ptr(_observations), par_transform(_par_transform),
		  cur_solution(_obj_func, _par_transform, *_observations), phiredswh_flag(false), save_next_jacobian(true), prior_info_ptr(_prior_info_ptr), jacobian(_jacobian), prev_phi_percent(0.0),
		  num_no_descent(0), regul_scheme_ptr(_regul_scheme_ptr), n_rotation_fac(_n_rotation_fac), description(_description)
{
	svd_package = new SVD_EIGEN();
}

void SVDSolver::set_svd_package(PestppOptions::SVD_PACK _svd_pack)
{
	if(_svd_pack == PestppOptions::PROPACK){
		delete svd_package;
		svd_package = new SVD_PROPACK;
	}
	else {
		delete svd_package;
		svd_package = new SVD_EIGEN;
	}
	svd_package->set_max_sing(svd_info.maxsing);
	svd_package->set_eign_thres(svd_info.eigthresh);
}


SVDSolver::~SVDSolver(void)
{
	delete svd_package;
}

double SVDSolver::calc_angle_deg(const Upgrade &upgrade1, const Upgrade &upgrade2)
{
	const double PI = 3.141592;
	double angle = acos(min(1.0, (upgrade1.uvec).dot(upgrade2.uvec))) * 180.0 / PI;
	return angle;
}

Parameters SVDSolver::apply_upgrade(const Parameters &init_numeric_pars,const Upgrade &upgrade, double scale)
{
	Parameters upgrade_pars = init_numeric_pars;

	// Add upgrade to init parameters
	// This does not check whether parameters remain in bounds
	for(int ip = 0; ip<upgrade.uvec.size(); ++ip)
	{
		const string &p_name = upgrade.par_name_vec[ip];
		auto it = upgrade_pars.find(upgrade.par_name_vec[ip]);
		assert(it != upgrade_pars.end());
		it->second += upgrade.uvec[ip]*upgrade.norm*scale;
	}
	//Impose previously frozen parameters
	for (auto &ipar : upgrade.frozen_numeric_pars)
	{
		upgrade_pars[ipar.first] = ipar.second;
	}
	return upgrade_pars;
}

void SVDSolver::update_upgrade(Upgrade &upgrade, const Parameters &base_pars, const Parameters &new_pars, 
							   const Parameters &frozen_pars)
{

	//build map of parameter names to index in original par_name_vec vector
	unordered_map<string, int> par_vec_name_to_idx;
	for(int i=0; i<upgrade.par_name_vec.size(); ++i)
	{
		par_vec_name_to_idx[upgrade.par_name_vec[i]] = i;
	}

	//compute new upgrade
	for(int i=0; i<upgrade.par_name_vec.size(); ++i)
	{
		auto &it_base = base_pars.find(upgrade.par_name_vec[i]);
		assert(it_base != base_pars.end());
		auto &it_new = new_pars.find(upgrade.par_name_vec[i]);
		assert(it_new != new_pars.end());
		upgrade.uvec[i] = it_new->second - it_base->second;
	}
	//tranfere previously frozen componets of the ugrade vector to upgrade.uvec
	upgrade.frozen_numeric_pars = frozen_pars;
	for(auto &ipar : frozen_pars)
	{
		const string &par_name = ipar.first;
		auto &it = par_vec_name_to_idx.find(par_name);
		assert(it != par_vec_name_to_idx.end());
		auto &it_base = base_pars.find(ipar.first);
		upgrade.uvec(it->second) = ipar.second - it_base->second;
	}
	upgrade.norm = upgrade.uvec.norm();
	if (upgrade.norm != 0) upgrade.uvec *= 1.0 /upgrade.norm;

}


ModelRun& SVDSolver::solve(RunManagerAbstract &run_manager, TerminationController &termination_ctl, int max_iter, 
	const Parameters &ctl_pars, ModelRun &optimum_run)
{
	ostream &os = file_manager.rec_ofstream();
	cur_solution.set_ctl_parameters(ctl_pars);
	// Start Solution iterations
	bool save_nextjac = false;
	for (int iter_num=1; iter_num<=max_iter;++iter_num) {
		int global_iter_num = termination_ctl.get_iteration_number()+1;
		cout << "OPTIMISATION ITERATION NUMBER: " << global_iter_num << endl;
		os   << "OPTIMISATION ITERATION NUMBER: " << global_iter_num << endl;
		cout << "  Iteration type: " << get_description() << endl;
		os   << "  Iteration type: " << get_description() << endl;
		cout << "  SVD Package: " << svd_package->description << endl;
		os   << "  SVD Package: " << svd_package->description << endl;
		os   << "    Model calls so far : " << run_manager.get_total_runs() << endl;
		iteration(run_manager, termination_ctl, false);
		// write files that get wrtten at the end of each iteration
		stringstream filename;
		string complete_filename;
		// rei file for this iteration
		filename << "rei" << global_iter_num;
		OutputFileWriter::write_rei(file_manager.open_ofile_ext(filename.str()), global_iter_num, 
			*(cur_solution.get_obj_func_ptr()->get_obs_ptr()), 
			cur_solution.get_obs(), *(cur_solution.get_obj_func_ptr()),
			cur_solution.get_ctl_pars());
		file_manager.close_file(filename.str());
		// par file for this iteration
		filename.str(""); // reset the stringstream
		filename << "par" << global_iter_num;
		OutputFileWriter::write_par(file_manager.open_ofile_ext(filename.str()), cur_solution.get_ctl_pars(), *(cur_solution.get_par_tran().get_offset_ptr()), 
				*(cur_solution.get_par_tran().get_scale_ptr()));
		file_manager.close_file(filename.str());
		// sen file for this iteration
		OutputFileWriter::append_sen(file_manager.sen_ofstream(), global_iter_num, jacobian, *(cur_solution.get_obj_func_ptr()), *par_group_info_ptr);
		if (save_nextjac) {
			jacobian.save();
		}
		if (!optimum_run.obs_valid() || cur_solution.get_phi() < optimum_run.get_phi())
		{
			optimum_run = cur_solution;
			// save new optimum parameters to .par file
			OutputFileWriter::write_par(file_manager.open_ofile_ext("par"), optimum_run.get_ctl_pars(), *(optimum_run.get_par_tran().get_offset_ptr()), 
				*(optimum_run.get_par_tran().get_scale_ptr()));
			file_manager.close_file("par");
			// save new optimum residuals to .rei file
			OutputFileWriter::write_rei(file_manager.open_ofile_ext("rei"), global_iter_num, 
			*(optimum_run.get_obj_func_ptr()->get_obs_ptr()), 
			optimum_run.get_obs(), *(optimum_run.get_obj_func_ptr()),
			optimum_run.get_ctl_pars());
			file_manager.close_file("rei");
			jacobian.save();
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


SVDSolver::Upgrade SVDSolver::calc_svd_upgrade_vec(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
	const Eigen::VectorXd &Residuals, const vector<string> &par_name_vec, const vector<string> &obs_name_vec,
	const Parameters &base_numeric_pars, const Parameters &freeze_numeric_pars, int &tot_sing_val)
{
	Upgrade upgrade;
	upgrade.frozen_numeric_pars = freeze_numeric_pars;
	upgrade.par_name_vec = par_name_vec;
	//create a vector which only contains the names of the unfrozen parameters
	vector<string> p_name_nf_vec(par_name_vec.size());
	{
		auto it = std::copy_if(par_name_vec.begin(), par_name_vec.end(), p_name_nf_vec.begin(),
			[&freeze_numeric_pars](const string &s) ->bool{return (freeze_numeric_pars.find(s)==freeze_numeric_pars.end());});
		p_name_nf_vec.resize(std::distance(p_name_nf_vec.begin(), it));
	}

	//Compute effect of frozen parameters on the residuals vector
	Parameters delta_freeze_pars = freeze_numeric_pars;
	delta_freeze_pars -= base_numeric_pars;
	VectorXd del_residuals;
	{
		vector<string>frz_par_name_vec = freeze_numeric_pars.get_keys();
		VectorXd frz_del_par_vec = delta_freeze_pars.get_data_eigen_vec(frz_par_name_vec);
		MatrixXd jac_frz = jacobian.get_matrix(frz_par_name_vec, obs_name_vec);
		del_residuals = (jac_frz) *  frz_del_par_vec;
	}

	MatrixXd jac = jacobian.get_matrix(p_name_nf_vec, obs_name_vec);
	MatrixXd SqrtQ_J = Q_sqrt * jac;
	VectorXd Sigma;
	MatrixXd U;
	MatrixXd Vt;
	svd_package->solve_ip(SqrtQ_J, Sigma, U, Vt);
	//calculate the number of singluar values above the threshold
	int max_sing = (p_name_nf_vec.size() < svd_info.maxsing) ? par_name_vec.size() : svd_info.maxsing;
	MatrixXd SqrtQ_J_inv =SVD_inv(U, Sigma, Vt, max_sing, svd_info.eigthresh, max_sing);
	VectorXd tmp_svd_uvec =(SqrtQ_J_inv * Q_sqrt) * (Residuals + del_residuals);

	//build map of parameter names to index in original par_name_vec vector
	unordered_map<string, int> par_vec_name_to_idx;
	for(int i=0; i<par_name_vec.size(); ++i)
	{
		par_vec_name_to_idx[par_name_vec[i]] = i;
	}

	upgrade.uvec.resize(par_name_vec.size());
	//tranfere newly computed componets of the ugrade vector to upgrade.svd_uvec
	for(int i=0; i<p_name_nf_vec.size(); ++i)
	{
		auto &it = par_vec_name_to_idx.find(p_name_nf_vec[i]);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = tmp_svd_uvec(i);
	}
	//tranfere previously frozen componets of the ugrade vector to upgrade.svd_uvec
	for(auto &ipar : delta_freeze_pars)
	{
		auto &it = par_vec_name_to_idx.find(ipar.first);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = ipar.second;
	}
	tot_sing_val = Sigma.size();
	// Calculate svd unit vector
	upgrade.norm = upgrade.uvec.norm();
	if (upgrade.norm != 0) upgrade.uvec *= 1.0 /upgrade.norm;
	return upgrade;
}


SVDSolver::Upgrade SVDSolver::calc_lambda_upgrade_vec(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
	const Eigen::VectorXd &Residuals, const vector<string> &par_name_vec, const vector<string> &obs_name_vec,
	const Parameters &base_numeric_pars, const Parameters &freeze_numeric_pars, int &tot_sing_val, double lambda)
{
	Upgrade upgrade;
	upgrade.frozen_numeric_pars = freeze_numeric_pars;
	upgrade.par_name_vec = par_name_vec;
	//create a vector which only contains the names of the unfrozen parameters
	vector<string> p_name_nf_vec(par_name_vec.size());
	{
		auto it = std::copy_if(par_name_vec.begin(), par_name_vec.end(), p_name_nf_vec.begin(),
			[&freeze_numeric_pars](const string &s) ->bool{return (freeze_numeric_pars.find(s)==freeze_numeric_pars.end());});
		p_name_nf_vec.resize(std::distance(p_name_nf_vec.begin(), it));
	}

	//Compute effect of frozen parameters on the residuals vector
	Parameters delta_freeze_pars = freeze_numeric_pars;
	delta_freeze_pars -= base_numeric_pars;
	VectorXd del_residuals;
	{
		vector<string>frz_par_name_vec = freeze_numeric_pars.get_keys();
		VectorXd frz_del_par_vec = delta_freeze_pars.get_data_eigen_vec(frz_par_name_vec);
		MatrixXd jac_frz = jacobian.get_matrix(frz_par_name_vec, obs_name_vec);
		del_residuals = (jac_frz) *  frz_del_par_vec;
	}

	MatrixXd jac = jacobian.get_matrix(p_name_nf_vec, obs_name_vec);
	MatrixXd ml_mat = jac.transpose() * Q_sqrt * Q_sqrt * jac + lambda * MatrixXd::Identity(jac.cols(), jac.cols());
	VectorXd Sigma;
	MatrixXd U;
	MatrixXd Vt;
	svd_package->solve_ip(ml_mat, Sigma, U, Vt);
	//calculate the number of singluar values above the threshold
	int max_sing = (p_name_nf_vec.size() < svd_info.maxsing) ? par_name_vec.size() : svd_info.maxsing;
	MatrixXd ml_mat_inv =SVD_inv(U, Sigma, Vt, max_sing, svd_info.eigthresh, max_sing);
	VectorXd tmp_svd_uvec =(ml_mat_inv * jac.transpose() * Q_sqrt) * (Residuals + del_residuals);

	//build map of parameter names to index in original par_name_vec vector
	unordered_map<string, int> par_vec_name_to_idx;
	for(int i=0; i<par_name_vec.size(); ++i)
	{
		par_vec_name_to_idx[par_name_vec[i]] = i;
	}

	upgrade.uvec.resize(par_name_vec.size());
	//tranfere newly computed componets of the ugrade vector to upgrade.svd_uvec
	for(int i=0; i<p_name_nf_vec.size(); ++i)
	{
		auto &it = par_vec_name_to_idx.find(p_name_nf_vec[i]);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = tmp_svd_uvec(i);
	}
	//tranfere previously frozen componets of the ugrade vector to upgrade.svd_uvec
	for(auto &ipar : delta_freeze_pars)
	{
		auto &it = par_vec_name_to_idx.find(ipar.first);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = ipar.second;
	}
	tot_sing_val = Sigma.size();
	// Calculate svd unit vector
	upgrade.norm = upgrade.uvec.norm();
	if (upgrade.norm != 0) upgrade.uvec *= 1.0 /upgrade.norm;
	return upgrade;
}

SVDSolver::Upgrade SVDSolver::calc_grad_upgrade_vec(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
	const Eigen::VectorXd &Residuals, const vector<string> &par_name_vec, const vector<string> &obs_name_vec,
	const Parameters &base_numeric_pars, const Parameters &freeze_numeric_pars, double l2_norm)
{
	Upgrade upgrade;
	upgrade.frozen_numeric_pars = freeze_numeric_pars;
	upgrade.par_name_vec = par_name_vec;

	//Compute effect of frozen parameters on the residuals vector
	Parameters delta_freeze_pars = freeze_numeric_pars;
	delta_freeze_pars -= base_numeric_pars;
	VectorXd del_residuals;
	{
		vector<string>frz_par_name_vec = freeze_numeric_pars.get_keys();
		VectorXd frz_del_par_vec = delta_freeze_pars.get_data_eigen_vec(frz_par_name_vec);
		MatrixXd jac_frz = jacobian.get_matrix(frz_par_name_vec, obs_name_vec);
		del_residuals = (jac_frz) *  frz_del_par_vec;
	}

	//create a vector which only contains the names of the unfrozen parameters
	vector<string> p_name_nf_vec(par_name_vec.size());
	{
		auto it = std::copy_if(par_name_vec.begin(), par_name_vec.end(), p_name_nf_vec.begin(),
			[&freeze_numeric_pars](const string &s) ->bool{return freeze_numeric_pars.find(s)==freeze_numeric_pars.end();});
		p_name_nf_vec.resize(std::distance(p_name_nf_vec.begin(), it));
	}
	MatrixXd jac = jacobian.get_matrix(p_name_nf_vec, obs_name_vec);
	MatrixXd SqrtQ_J = Q_sqrt * jac;
	VectorXd tmp_grad_vec = jac.transpose() * Q_sqrt * Q_sqrt * (Residuals+del_residuals);

	//scale componets of new upgeade vector so that the entire vector have a norm of l2_norm
	double nf_norm_sqr = pow(tmp_grad_vec.norm(), 2.0);
	double freeze_norm_sqr = pow(delta_freeze_pars.l2_norm(), 2.0);
	if (nf_norm_sqr != 0) {
		tmp_grad_vec *= sqrt((pow(l2_norm, 2.0) - freeze_norm_sqr)/nf_norm_sqr);
	}

	//build map of parameter names to index in original par_name_vec vector
	unordered_map<string, int> par_vec_name_to_idx;
	for(int i=0; i<par_name_vec.size(); ++i)
	{
		par_vec_name_to_idx[par_name_vec[i]] = i;
	}

	upgrade.uvec.resize(par_name_vec.size());
	//tranfere newly computed componets of the ugrade vector to upgrade.svd_uvec
	for(int i=0; i<p_name_nf_vec.size(); ++i)
	{
		auto &it = par_vec_name_to_idx.find(p_name_nf_vec[i]);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = tmp_grad_vec(i);
	}
	//tranfere previously frozen componets of the ugrade vector to upgrade.uvec
	for(auto &ipar : delta_freeze_pars)
	{
		auto &it = par_vec_name_to_idx.find(ipar.first);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = ipar.second;
	}
	upgrade.norm = upgrade.uvec.norm();
	if (upgrade.norm != 0) upgrade.uvec *= 1.0 /upgrade.norm;
	return upgrade;
}

SVDSolver::Upgrade SVDSolver::get_rotated_upgrade(const Upgrade &upgrade_svd, const Upgrade &upgrade_grad,
		const Parameters &base_numeric_pars, double rot_fac,
		const Parameters &freeze_numeric_pars, double l2_norm)
{
	Upgrade upgrade;
	const double PI = 3.141592;
	assert(upgrade_svd.par_name_vec == upgrade_grad.par_name_vec);
	upgrade.par_name_vec = upgrade_svd.par_name_vec;
	upgrade.uvec = ((1.0 - rot_fac) * upgrade_svd.uvec + rot_fac * upgrade_grad.uvec);

	//Compute Frozen component of the upgrade vector
	Parameters delta_freeze_pars = freeze_numeric_pars;
	delta_freeze_pars -= base_numeric_pars;

	//Compute unfrozen component of the upgrade vector
	Parameters delta_unfrozen_pars;
	for(int i=0; i<upgrade.par_name_vec.size(); ++i)
	{
		auto iter = freeze_numeric_pars.find(upgrade.par_name_vec[i]);
		if (iter == freeze_numeric_pars.end()) {
			delta_unfrozen_pars[upgrade.par_name_vec[i]] = upgrade.uvec(i);
		}
	}

	//scale componets of new upgeade vector so that the entire vector have a norm of l2_norm
	double nf_norm_sqr = pow(delta_unfrozen_pars.l2_norm(), 2.0);
	double freeze_norm_sqr = pow(delta_freeze_pars.l2_norm(), 2.0);
	if (nf_norm_sqr != 0) {
		delta_unfrozen_pars *= sqrt((pow(l2_norm, 2.0) - freeze_norm_sqr)/nf_norm_sqr);
	}

	//build map of parameter names to index in original par_name_vec vector
	unordered_map<string, int> par_vec_name_to_idx;
	for(int i=0; i<upgrade.par_name_vec.size(); ++i)
	{
		par_vec_name_to_idx[upgrade.par_name_vec[i]] = i;
	}

	//tranfere newly computed componets of the ugrade vector
	for(auto &ipar : delta_unfrozen_pars)
	{
		auto &it = par_vec_name_to_idx.find(ipar.first);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = ipar.second;
	}
	//tranfere previously frozen componets of the ugrade vector to upgrade.uvec
	for(auto &ipar : delta_freeze_pars)
	{
		auto &it = par_vec_name_to_idx.find(ipar.first);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = ipar.second;
	}
	upgrade.norm = upgrade.uvec.norm();
	if (upgrade.norm != 0) upgrade.uvec *= 1.0 /upgrade.norm;
	return upgrade;
}

void SVDSolver::iteration(RunManagerAbstract &run_manager, TerminationController &termination_ctl, bool calc_init_obs)
{
	ostream &os = file_manager.rec_ofstream();
	const double PI = 3.141592;
	VectorXd residuals_vec;
	ModelRun base_run(cur_solution);
	vector<string> obs_names_vec;


	// Calculate Jacobian
	if (!cur_solution.obs_valid() || calc_init_obs == true) {
		 calc_init_obs = true;
	 }
	cout << "  calculating jacobian... ";
	jacobian.calculate(cur_solution, cur_solution.get_numeric_pars().get_keys(),  cur_solution.get_obs_template().get_keys(),
			*par_group_info_ptr, *ctl_par_info_ptr,run_manager,  *prior_info_ptr,
			phiredswh_flag, calc_init_obs);
	cout << endl;
	cout << "  computing upgrade vectors... " << endl;

	//Freeze Parameter for which the jacobian could not be calculated
	auto &failed_jac_pars = jacobian.get_failed_parameter_names();
	//auto  freeze_pars = cur_solution.get_numeric_pars().get_subset(failed_jac_pars.begin(), failed_jac_pars.end());
	//cur_solution.freeze_parameters(freeze_pars);

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

	//Build model runs
	run_manager.reinitialize();
	vector<double> rot_fac_vec;
	vector<double> rot_angle_vec;
	vector<double> magnitude_vec;


	//// compute svd update vector
	//double target_upgrade_norm = 0;
	//Upgrade svd_upgrade;
	//Parameters new_numeric_pars;
	//const Parameters &base_run_numeric_pars = base_run.get_numeric_pars();
	//vector<string> par_name_vec = base_run_numeric_pars.get_keys();
	//Parameters frozen_pars;
	//int tot_sing_val;
	//LimitType limit_type = LimitType::NONE;
	//while (true)
	//{
	//	svd_upgrade = calc_svd_upgrade_vec(jacobian, Q_sqrt, residuals_vec, par_name_vec, obs_names_vec,
	//		base_run.get_numeric_pars(), frozen_pars, tot_sing_val);
	//	new_numeric_pars = apply_upgrade(base_run_numeric_pars, svd_upgrade, 1.0);
	//	Parameters new_frozen_pars = limit_parameters_ip(base_run_numeric_pars, new_numeric_pars, limit_type, frozen_pars);
	//	if (new_frozen_pars.size() > 0 && (limit_type == LimitType::UBND || limit_type == LimitType::LBND) )
	//	{
	//		frozen_pars.insert(new_frozen_pars.begin(), new_frozen_pars.end());
	//	}
	//	update_upgrade(svd_upgrade, base_run_numeric_pars, new_numeric_pars, frozen_pars);
	//	if (limit_type != LimitType::UBND && limit_type != LimitType::LBND) break;
	//	if (frozen_pars.size() == new_numeric_pars.size()) break; // everything is frozen
	//}
	////add run svd run to run manager
	//run_manager.add_run(base_run.get_par_tran().numeric2model_cp(new_numeric_pars));
	//rot_fac_vec.push_back(0);
	//rot_angle_vec.push_back(0);
	//magnitude_vec.push_back(svd_upgrade.norm);
	//target_upgrade_norm = svd_upgrade.norm;


	//Marquardt Lambda Update Vector
	Parameters frozen_pars;
	Parameters new_numeric_pars;
	Upgrade ml_upgrade;
	int tot_sing_val;
	LimitType limit_type = LimitType::NONE;
	const Parameters &base_run_numeric_pars = base_run.get_numeric_pars();
	vector<string> par_name_vec = base_run_numeric_pars.get_keys();

	double tmp_lambda[] = {1.0, 10.0, 100.0, 1000.0};
	vector<double> lambda_vec(tmp_lambda, tmp_lambda+sizeof(tmp_lambda)/sizeof(double));
	for (double i_lambda : lambda_vec)
	{
		frozen_pars.clear();
		while (true)
		{
			ml_upgrade = calc_lambda_upgrade_vec(jacobian, Q_sqrt, residuals_vec, par_name_vec, obs_names_vec,
				base_run.get_numeric_pars(), frozen_pars, tot_sing_val, i_lambda);
			new_numeric_pars = apply_upgrade(base_run_numeric_pars, ml_upgrade, 1.0);
			Parameters new_frozen_pars = limit_parameters_ip(base_run_numeric_pars, new_numeric_pars, limit_type, frozen_pars);
			if (new_frozen_pars.size() > 0 && (limit_type == LimitType::UBND || limit_type == LimitType::LBND) )
			{
				frozen_pars.insert(new_frozen_pars.begin(), new_frozen_pars.end());
			}
			update_upgrade(ml_upgrade, base_run_numeric_pars, new_numeric_pars, frozen_pars);
			if (limit_type != LimitType::UBND && limit_type != LimitType::LBND) break;
			if (frozen_pars.size() == new_numeric_pars.size()) break; // everything is frozen
		}
		//add run svd run to run manager
		run_manager.add_run(base_run.get_par_tran().numeric2model_cp(new_numeric_pars));
		rot_fac_vec.push_back(0);
		rot_angle_vec.push_back(0);
		magnitude_vec.push_back(ml_upgrade.norm);
	}



	////check for another run svd run with limited parameters frozen
	//if (limit_type != LimitType::NONE)
	//{
	//	Upgrade svd_frz_limited_upgrade;
	//	while (true)
	//	{
	//		svd_frz_limited_upgrade = calc_svd_upgrade_vec(jacobian, Q_sqrt, residuals_vec, par_name_vec, obs_names_vec,
	//			base_run.get_numeric_pars(), frozen_pars, tot_sing_val);
	//		new_numeric_pars = apply_upgrade(base_run_numeric_pars, svd_frz_limited_upgrade, 1.0);
	//		Parameters new_frozen_pars = limit_parameters_ip(base_run_numeric_pars, new_numeric_pars, limit_type, frozen_pars);
	//		if (new_frozen_pars.size() > 0)
	//		{
	//			frozen_pars.insert(new_frozen_pars.begin(), new_frozen_pars.end());
	//		}
	//		update_upgrade(svd_frz_limited_upgrade, base_run_numeric_pars, new_numeric_pars, frozen_pars);
	//		if (limit_type == LimitType::NONE) break;
	//		if (frozen_pars.size() == new_numeric_pars.size()) break; // everything is frozen
	//	}
	//	//add second svd run with limited parameters frozen to run manager
	//	run_manager.add_run(base_run.get_par_tran().numeric2model_cp(new_numeric_pars));
	//	rot_fac_vec.push_back(0);
	//	rot_angle_vec.push_back(0);
	//	magnitude_vec.push_back(svd_frz_limited_upgrade.norm);
	//	target_upgrade_norm = svd_frz_limited_upgrade.norm;
	//}


	////compute gradient descent upgrade vector
	//frozen_pars.clear();
	//Upgrade grad_upgrade;
	//Parameters grad_numeric_pars;
	//limit_type = LimitType::NONE;
	//double tmp_alpha = alpha;
	//if (alpha <= 0)
	//{
	//	alpha = .5;
	//}
	//else if (precent_grad_phi > 100)
	//{
	//	alpha /= 2.0;
	//}
	//else if (alpha < alpha_prev && precent_grad_phi < precent_grad_phi_prev)
	//{
	//	alpha *= 0.80;
	//}
	//else
	//{
	//	alpha *= 1.2;
	//}
	//alpha_prev = tmp_alpha;
	//precent_grad_phi_prev = precent_grad_phi;
	//while (true)
	//{
	//	grad_upgrade = calc_grad_upgrade_vec(jacobian, Q_sqrt, residuals_vec, par_name_vec, obs_names_vec,
	//		base_run.get_numeric_pars(), frozen_pars, svd_upgrade.norm);
	//	grad_numeric_pars = apply_upgrade(base_run_numeric_pars, grad_upgrade, alpha);
	//	Parameters new_frozen_pars = limit_parameters_ip(base_run_numeric_pars, grad_numeric_pars, limit_type, frozen_pars);
	//	frozen_pars.insert(new_frozen_pars.begin(), new_frozen_pars.end());
	//	update_upgrade(grad_upgrade, base_run_numeric_pars, grad_numeric_pars, frozen_pars);
	//	if (limit_type != LimitType::UBND && limit_type != LimitType::LBND) break;
	//	if (frozen_pars.size() == new_numeric_pars.size()) break; // everything is frozen
	//}

	os << endl;
	os << "      SVD information:" << endl;
	//os << "        number of singular values used: " << upgrade.n_sing_val_used << "/" << upgrade.tot_sing_val << endl;
	//os << "        upgrade vector magnitude (without limits or bounds) = " << svd_upgrade.norm << endl;
	//os << "        angle to direction of greatest descent: ";
	//os <<  calc_angle_deg(svd_upgrade, grad_upgrade) << " deg" << endl;
	os << endl;

	//add rotated runs to run manager
	//for(int i=1; i<n_rotation_fac-1; ++i) 
	//{
	//	double rot_fac = (i == 0) ? 0 : i / double(n_rotation_fac-1);
	//	Upgrade tmp_upgrade;
	//	limit_type = LimitType::NONE;
	//	frozen_pars.clear();
	//	while (true)
	//	{
	//		tmp_upgrade = get_rotated_upgrade(svd_upgrade, grad_upgrade, base_run.get_numeric_pars(),
	//				rot_fac, frozen_pars, svd_upgrade.norm);
	//		new_numeric_pars = apply_upgrade(base_run_numeric_pars, tmp_upgrade, 1.0);
	//		Parameters new_frozen_pars = limit_parameters_ip(base_run_numeric_pars, new_numeric_pars, limit_type, frozen_pars);
	//		frozen_pars.insert(new_frozen_pars.begin(), new_frozen_pars.end());
	//		update_upgrade(tmp_upgrade, base_run_numeric_pars, new_numeric_pars, frozen_pars);
	//		if (limit_type != LimitType::UBND && limit_type != LimitType::LBND) break;
	//		if (frozen_pars.size() == new_numeric_pars.size()) break; // everything is frozen
	//	}
	//	run_manager.add_run(base_run.get_par_tran().numeric2model_cp(new_numeric_pars));
	//	rot_fac_vec.push_back(rot_fac);
	//	rot_angle_vec.push_back(calc_angle_deg(svd_upgrade, tmp_upgrade));
	//	magnitude_vec.push_back(tmp_upgrade.norm);
	//}

	////add gradient descent run to run manager
	//run_manager.add_run(base_run.get_par_tran().numeric2model_cp(grad_numeric_pars));
	//rot_fac_vec.push_back(1.0);
	//rot_angle_vec.push_back(calc_angle_deg(svd_upgrade, grad_upgrade));
	//magnitude_vec.push_back(grad_upgrade.norm);

	// process model runs
	cout << "  testing upgrade vectors... ";
	run_manager.run();
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

	//// Print frozen parameter information
	//if (svd_upgrade.frozen_numeric_pars.size() > 0)
	//{
	//	Parameters ctl_frz_par = svd_upgrade.frozen_numeric_pars;
	//	vector<string> keys = ctl_frz_par.get_keys();
	//	std::sort(keys.begin(), keys.end());
	//	base_run.get_par_tran().numeric2ctl_ip(ctl_frz_par);
	//	os << endl;
	//	os << "  Parameters frozen during SVD upgrade:" << endl;
	//	// ctl_frz_par also contains fixed parameters so we need to loop through the keys 
	//	// of the original numeric parameters
	//	for (auto &ikey : keys)
	//	{
	//		auto iter = ctl_frz_par.find(ikey);
	//		if (iter != ctl_frz_par.end())
	//		{
	//			os << "    " << iter->first << " frozen at " << iter->second << endl;
	//		}
	//	}
	//}
	//if (grad_upgrade.frozen_numeric_pars.size() > 0)
	//{
	//	Parameters ctl_frz_par = grad_upgrade.frozen_numeric_pars;
	//	vector<string> keys = ctl_frz_par.get_keys();
	//	std::sort(keys.begin(), keys.end());
	//	base_run.get_par_tran().numeric2ctl_ip(ctl_frz_par);
	//	os << endl;
	//	os << "  Parameters frozen during Gradient Descent upgrade:" << endl;
	//	// ctl_frz_par also contains fixed parameters so we need to loop through the keys 
	//	// of the original numeric parameters
	//	for (auto &ikey : keys)
	//	{
	//		auto iter = ctl_frz_par.find(ikey);
	//		if (iter != ctl_frz_par.end())
	//		{
	//			os << "    " << iter->first << " frozen at " << iter->second << endl;
	//		}
	//	}
	//}

	// clean up run_manager memory
	run_manager.free_memory();

	// set flag to switch to central derivatives next iteration
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
	// Clear Frozen Parameters
	cur_solution.thaw_parameters();
}

Parameters SVDSolver::limit_parameters_ip(const Parameters &init_numeric_pars, Parameters &upgrade_numeric_pars, 
										  LimitType &limit_type, const Parameters &frozen_pars)
{
	map<string, LimitType> limit_type_map;
	const string *name;
	const ParameterRec *p_info;
	double p_init;
	double p_upgrade;
	double p_limit;
	double b_facorg_lim;
	Parameters bound_parameters;
	pair<bool, double> par_limit;
	Parameters ctl_parameters_at_limit;
	Parameters upgrade_ctl_pars  = par_transform.numeric2ctl_cp(upgrade_numeric_pars);

	Parameters init_ctl_pars = init_numeric_pars;
	//remove frozen parameters from numeric parameters
	auto end_iter = init_ctl_pars.end();
	for(auto iter = init_ctl_pars.begin(); iter != end_iter;) {
		if (frozen_pars.find((*iter).first) != frozen_pars.end())
		{
		   init_ctl_pars.erase(iter++);
		}
		else 
		{
			++iter;
	   }
	}
	par_transform.numeric2ctl_ip(init_ctl_pars);

	for(auto &ipar : init_ctl_pars)
	{
		par_limit = pair<bool, double>(false, 0.0);
		name = &(ipar.first);  // parameter name
		p_init = ipar.second; // inital parameter value
		auto pu_iter = upgrade_ctl_pars.find(*name);
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
			limit_type_map[*name] = LimitType::REL;
		}

		// Check Factor Change Limit
		else if(p_info->chglim == "FACTOR") {
			if (b_facorg_lim > 0 && p_upgrade < b_facorg_lim/ctl_info->facparmax ) 
			{
				par_limit.first = true;
				par_limit.second = b_facorg_lim / ctl_info->facparmax;
				limit_type_map[*name] = LimitType::FACT;
			}
			else if (b_facorg_lim > 0 && p_upgrade > b_facorg_lim*ctl_info->facparmax ) 
			{
				par_limit.first = true;
				par_limit.second = b_facorg_lim * ctl_info->facparmax;
				limit_type_map[*name] = LimitType::FACT;
			}
			else if (b_facorg_lim < 0 && p_upgrade < b_facorg_lim*ctl_info->facparmax)
			{
				par_limit.first = true;
				par_limit.second = b_facorg_lim * ctl_info->facparmax;
				limit_type_map[*name] = LimitType::FACT;
			}
			else if (b_facorg_lim < 0 && p_upgrade > b_facorg_lim/ctl_info->facparmax)
			{
				par_limit.first = true;
				par_limit.second = b_facorg_lim / ctl_info->facparmax;
				limit_type_map[*name] = LimitType::FACT;
			}
		}
		// Check parameter upper bound
		if((!par_limit.first && p_upgrade > p_info->ubnd) || 
			(par_limit.first && par_limit.second > p_info->ubnd)) {
				par_limit.first = true;
				par_limit.second = p_info->ubnd;
				bound_parameters[*name] = p_info->ubnd;
				limit_type_map[*name] = LimitType::UBND;
		}
		// Check parameter lower bound
		else if((!par_limit.first && p_upgrade < p_info->lbnd) || 
			(par_limit.first && par_limit.second < p_info->lbnd)) {
				par_limit.first = true;
				par_limit.second = p_info->lbnd;
				bound_parameters[*name] = p_info->lbnd;
				limit_type_map[*name] = LimitType::LBND;
		}
		// Add any limited parameters to model_parameters_at_limit
		if (par_limit.first) {
			ctl_parameters_at_limit.insert(*name, par_limit.second);
		}
	}
	// Calculate most stringent limit factor on a PEST parameter
	double limit_factor= 1.0;
	double tmp_limit;
	string limit_parameter_name = "";
	Parameters numeric_parameters_at_limit = par_transform.model2numeric_cp(ctl_parameters_at_limit);
	for(auto &ipar : numeric_parameters_at_limit)
	{
			name = &(ipar.first);
			p_limit = ipar.second;
			p_init = init_numeric_pars.get_rec(*name);
			p_upgrade = upgrade_numeric_pars.get_rec(*name);
			tmp_limit = (p_limit - p_init) / (p_upgrade - p_init);
			if (tmp_limit < limit_factor)  {
				limit_factor = (p_limit - p_init) / (p_upgrade - p_init);
				limit_parameter_name = *name;
			}
	}
	// Apply limit factor to PEST upgrade parameters
	if (limit_factor != 1.0)
	{
		for(auto &ipar : upgrade_numeric_pars)
		{
			name = &(ipar.first);
			p_init = init_numeric_pars.get_rec(*name);
			ipar.second = p_init + (ipar.second - p_init) *  limit_factor;
		}
	}
	// Impose frozen Parameters
	for (auto &ipar : frozen_pars)
	{
		auto iter = upgrade_numeric_pars.find(ipar.first);
		if (iter != upgrade_numeric_pars.end())
		{
			upgrade_numeric_pars[ipar.first] = ipar.second;
		}
	}
	limit_type = LimitType::NONE;
	Parameters::iterator iter;
	iter = numeric_parameters_at_limit.find(limit_parameter_name);
	if (iter != numeric_parameters_at_limit.end()) {
		bound_parameters.clear();
		bound_parameters[limit_parameter_name] = numeric_parameters_at_limit[limit_parameter_name];
		auto iter = limit_type_map.find(limit_parameter_name);
		assert(iter != limit_type_map.end());
		limit_type = iter->second;
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