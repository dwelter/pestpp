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
		const Regularization *_regul_scheme_ptr, int _max_freeze_iter, OutputFileWriter &_output_file_writer, Pest::RestartCommand &_restart_command, const string &_description)
		: ctl_info(_ctl_info), svd_info(_svd_info), par_group_info_ptr(_par_group_info_ptr), ctl_par_info_ptr(_ctl_par_info_ptr), obs_info_ptr(_obs_info), obj_func(_obj_func),
		  file_manager(_file_manager), observations_ptr(_observations), par_transform(_par_transform),
		  cur_solution(_obj_func, *_observations), phiredswh_flag(false), save_next_jacobian(true), prior_info_ptr(_prior_info_ptr), jacobian(_jacobian), prev_phi_percent(0.0),
		  num_no_descent(0), regul_scheme_ptr(_regul_scheme_ptr), max_freeze_iter(_max_freeze_iter), output_file_writer(_output_file_writer), description(_description), best_lambda(20.0), restart_command(_restart_command)
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
		const auto &it_base = base_pars.find(upgrade.par_name_vec[i]);
		assert(it_base != base_pars.end());
		const auto &it_new = new_pars.find(upgrade.par_name_vec[i]);
		assert(it_new != new_pars.end());
		upgrade.uvec[i] = it_new->second - it_base->second;
	}
	//tranfere previously frozen componets of the ugrade vector to upgrade.uvec
	upgrade.frozen_numeric_pars = frozen_pars;
	for(auto &ipar : frozen_pars)
	{
		const string &par_name = ipar.first;
		const auto &it = par_vec_name_to_idx.find(par_name);
		assert(it != par_vec_name_to_idx.end());
		const auto &it_base = base_pars.find(ipar.first);
		upgrade.uvec(it->second) = ipar.second - it_base->second;
	}
	upgrade.norm = upgrade.uvec.norm();
	if (upgrade.norm != 0) upgrade.uvec *= 1.0 /upgrade.norm;

}


ModelRun& SVDSolver::solve(RunManagerAbstract &run_manager, TerminationController &termination_ctl, int max_iter, 
	ModelRun &cur_run, ModelRun &optimum_run)
{
	ostream &os = file_manager.rec_ofstream();
	ostream &fout_restart = file_manager.get_ofstream("rst");
	cur_solution = cur_run;
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
		fout_restart << "start_iteration " << global_iter_num << endl;
		termination_ctl.save_state(fout_restart);
		// write head for SVD file
		output_file_writer.write_svd_iteration(global_iter_num);
		iteration(run_manager, termination_ctl, false);

		// write files that get wrtten at the end of each iteration
		stringstream filename;
		string complete_filename;

		// rei file for this iteration
		filename << "rei" << global_iter_num;
		output_file_writer.write_rei(file_manager.open_ofile_ext(filename.str()), global_iter_num, 
			*(cur_solution.get_obj_func_ptr()->get_obs_ptr()), 
			cur_solution.get_obs(), *(cur_solution.get_obj_func_ptr()),
			cur_solution.get_ctl_pars());
		file_manager.close_file(filename.str());
		// par file for this iteration
		filename.str(""); // reset the stringstream
		filename << "par" << global_iter_num;
		output_file_writer.write_par(file_manager.open_ofile_ext(filename.str()), cur_solution.get_ctl_pars(), *(par_transform.get_offset_ptr()), 
				*(par_transform.get_scale_ptr()));
		file_manager.close_file(filename.str());
		// sen file for this iteration
		output_file_writer.append_sen(file_manager.sen_ofstream(), global_iter_num, jacobian, *(cur_solution.get_obj_func_ptr()), *par_group_info_ptr);
		if (save_nextjac) {
			jacobian.save();
		}
		if (!optimum_run.obs_valid() || cur_solution.get_phi() < optimum_run.get_phi())
		{
			optimum_run.set_ctl_parameters(cur_solution.get_ctl_pars());
			optimum_run.set_observations(cur_solution.get_obs());
			// save new optimum parameters to .par file
			output_file_writer.write_par(file_manager.open_ofile_ext("par"), optimum_run.get_ctl_pars(), *(par_transform.get_offset_ptr()), 
				*(par_transform.get_scale_ptr()));
			file_manager.close_file("par");
			// save new optimum residuals to .rei file
			output_file_writer.write_rei(file_manager.open_ofile_ext("rei"), global_iter_num, 
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

VectorXd SVDSolver::calc_residual_corrections(const Jacobian &jacobian, const Parameters &del_numeric_pars, 
							   const vector<string> obs_name_vec)
{
	VectorXd del_residuals;
	vector<string>frz_par_name_vec = del_numeric_pars.get_keys();
	//remove the parameters for which the jaocbian could not be computed
	const set<string> &failed_jac_par_names = jacobian.get_failed_parameter_names();
	auto end_iter = remove_if(frz_par_name_vec.begin(), frz_par_name_vec.end(), 	
		[&failed_jac_par_names](string &str)->bool{return failed_jac_par_names.find(str)!=failed_jac_par_names.end();});
	frz_par_name_vec.resize(std::distance(frz_par_name_vec.begin(), end_iter));


	VectorXd frz_del_par_vec = del_numeric_pars.get_data_eigen_vec(frz_par_name_vec);

	MatrixXd jac_frz = jacobian.get_matrix(obs_name_vec, frz_par_name_vec);
	del_residuals = (jac_frz) *  frz_del_par_vec;
	return del_residuals;
}

SVDSolver::Upgrade SVDSolver::calc_lambda_upgrade_vec(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
	const Eigen::VectorXd &Residuals, const vector<string> &par_name_vec, const vector<string> &obs_name_vec,
	const Parameters &base_numeric_pars, const Parameters &freeze_numeric_pars, int &tot_sing_val, 
	double lambda, MarquardtMatrix marquardt_type)
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
	VectorXd Sigma;
	Eigen::SparseMatrix<double> U;
	Eigen::SparseMatrix<double> Vt;
	Parameters delta_freeze_pars = freeze_numeric_pars;
	delta_freeze_pars -= base_numeric_pars;
	VectorXd del_residuals = calc_residual_corrections(jacobian, delta_freeze_pars, obs_name_vec);
	Eigen::SparseMatrix<double> q_sqrt = Q_sqrt.get_sparse_matrix(obs_name_vec);
	{ 
		Eigen::SparseMatrix<double> jac = jacobian.get_matrix(obs_name_vec, p_name_nf_vec);
		Eigen::SparseMatrix<double> SqrtQ_J = q_sqrt * jac;
		svd_package->solve_ip(SqrtQ_J, Sigma, U, Vt);
	}

	//calculate the number of singluar values above the threshold
	int num_sing_used = 0;
	int max_sing = (p_name_nf_vec.size() < svd_info.maxsing) ? par_name_vec.size() : svd_info.maxsing;
	for (int i = 0; i < max_sing; ++i) {
		if (Sigma(i) != 0 && i<max_sing && Sigma(i) / Sigma(0) > svd_info.eigthresh) {
			++num_sing_used;
		}
	}

	//Trim the Matricies based on the number of singular values to be used
	Sigma = VectorXd(Sigma.head(num_sing_used));

	Vt = Eigen::SparseMatrix<double>(Vt.topRows(num_sing_used));
	U = Eigen::SparseMatrix<double>(U.leftCols(num_sing_used));

	//Only add lambda to singular values above the threshhold 
	if (marquardt_type == MarquardtMatrix::IDENT)
	{
		Sigma = Sigma.array() + lambda;
	}
	else
	{
		//this needs checking 
		Sigma = Sigma.array() + (Sigma.cwiseProduct(Sigma).array() + lambda * lambda).sqrt();
	}
	VectorXd Sigma_inv = Sigma.array().inverse();

	Eigen::SparseVector<double> tmp_residual = (Residuals + del_residuals).sparseView();
	Eigen::VectorXd tmp_svd_uvec;
	tmp_svd_uvec = Vt.transpose() * Sigma_inv.asDiagonal() * U.transpose() * q_sqrt  * (Residuals + del_residuals);
	output_file_writer.write_svd(Sigma, Vt, num_sing_used, lambda, freeze_numeric_pars);

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
		const auto &it = par_vec_name_to_idx.find(p_name_nf_vec[i]);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = tmp_svd_uvec(i);
	}
	//tranfere previously frozen componets of the ugrade vector to upgrade.svd_uvec
	for(auto &ipar : delta_freeze_pars)
	{
		const auto &it = par_vec_name_to_idx.find(ipar.first);
		assert(it != par_vec_name_to_idx.end());
		upgrade.uvec(it->second) = ipar.second;
	}

	// Calculate svd unit vector
	upgrade.norm = upgrade.uvec.norm();
	if (upgrade.norm != 0) upgrade.uvec *= 1.0 /upgrade.norm;
	return upgrade;
}

void SVDSolver::iteration(RunManagerAbstract &run_manager, TerminationController &termination_ctl, bool calc_init_obs)
{
	ostream &fout_restart = file_manager.get_ofstream("rst");
	fout_restart << "native_par_iteration" << endl;
	ostream &os = file_manager.rec_ofstream();
	const double PI = 3.141592;
	VectorXd residuals_vec;
	set<string> out_ofbound_pars;

	vector<string> obs_names_vec = cur_solution.get_obs_template().get_keys();
	vector<string> numeric_parname_vec = par_transform.ctl2numeric_cp(cur_solution.get_ctl_pars()).get_keys();

	if (restart_command == Pest::RestartCommand::REUSE_JACOBIAN)
	{
		restart_command = Pest::RestartCommand::NONE;
		cout << "  reading previosuly computed jacobian... ";
		{
			jacobian.read(file_manager.build_filename("jco"), *prior_info_ptr);
		}

		cout << endl << endl;
		cout << "  running the model once with the current parameters... ";
		run_manager.reinitialize(file_manager.build_filename("rnu"));
		int run_id = run_manager.add_run( par_transform.ctl2model_cp(cur_solution.get_ctl_pars()));
		run_manager.run();
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(run_id, tmp_pars, tmp_obs);
		if (success)
		{
			par_transform.model2ctl_ip(tmp_pars);
			cur_solution.update_ctl(tmp_pars, tmp_obs);
			goto restart_reuse_jacoboian;
		}
		else
		{
			throw(PestError("Error: Base parameter run failed.  Can not continue."));
		}
	}

	// Calculate Jacobian
	if (!cur_solution.obs_valid() || calc_init_obs == true) {
		 calc_init_obs = true;
	 }
	cout << "  calculating jacobian... ";
	jacobian.build_runs(cur_solution, numeric_parname_vec, par_transform,
			*par_group_info_ptr, *ctl_par_info_ptr, run_manager, out_ofbound_pars,
			phiredswh_flag, calc_init_obs);
	jacobian.make_runs(run_manager);
	jacobian.process_runs(cur_solution, numeric_parname_vec, obs_names_vec, par_transform,
			*par_group_info_ptr, *ctl_par_info_ptr, run_manager, *prior_info_ptr, out_ofbound_pars,
			phiredswh_flag, calc_init_obs);

	restart_reuse_jacoboian:
	cout << endl;


	//Freeze Parameter for which the jacobian could not be calculated
	auto &failed_jac_pars_names = jacobian.get_failed_parameter_names();
	auto  failed_jac_pars = cur_solution.get_ctl_pars().get_subset(failed_jac_pars_names.begin(), failed_jac_pars_names.end());

	cout << endl;
	cout << "  computing upgrade vectors... " << endl;
	// update regularization weight factor
	double tikhonov_weight = regul_scheme_ptr->get_weight(cur_solution);
	// write out report for starting phi
	obj_func->phi_report(os, cur_solution.get_obs(), cur_solution.get_ctl_pars(), tikhonov_weight);
	// write failed jacobian parameters out
	if (failed_jac_pars.size() > 0)
	{
		os << endl;
		os << "  the following parameters have been frozen as the runs to compute their derivatives failed: " << endl;
		for (auto &ipar : failed_jac_pars)
		{
			os << "    " << ipar.first << " frozen at " << ipar.second << endl; 
		}
	}
	os << endl;
	// populate vectors with sorted observations (standard and prior info) and parameters
	{
		vector<string> prior_info_names = prior_info_ptr->get_keys();
		obs_names_vec.insert(obs_names_vec.end(), prior_info_names.begin(), prior_info_names.end());
	}
	// build weights matrix sqrt(Q)
	QSqrtMatrix Q_sqrt(obs_info_ptr, prior_info_ptr, tikhonov_weight);
	//build residuals vector
	residuals_vec = -1.0 * stlvec_2_egienvec(cur_solution.get_residuals_vec(obs_names_vec));

	//Build model runs
	run_manager.reinitialize(file_manager.build_filename("rnu"));
	vector<double> rot_fac_vec;
	vector<double> rot_angle_vec;
	vector<double> magnitude_vec;
	vector<Parameters> frozen_par_vec;

	//Marquardt Lambda Update Vector
	Parameters frozen_derivative_pars;
	Upgrade ml_upgrade;
	int tot_sing_val;
	LimitType limit_type = LimitType::NONE;
	const Parameters &base_run_numeric_pars =  par_transform.ctl2numeric_cp(cur_solution.get_ctl_pars());
	vector<string> par_name_vec = base_run_numeric_pars.get_keys();

	double tmp_lambda[] = {0.1, 1.0, 10.0, 100.0, 1000.0};
	vector<double> lambda_vec(tmp_lambda, tmp_lambda+sizeof(tmp_lambda)/sizeof(double));
	lambda_vec.push_back(best_lambda);
	lambda_vec.push_back(best_lambda / 2.0);
	lambda_vec.push_back(best_lambda * 2.0);
	std::sort(lambda_vec.begin(), lambda_vec.end());
	auto iter = std::unique(lambda_vec.begin(), lambda_vec.end());
	lambda_vec.resize(std::distance(lambda_vec.begin(), iter));
	int max_freeze_iter = 1;
	int i_update_vec = 0;
	stringstream message;
	for (double i_lambda : lambda_vec)
	{
		std::cout << string(message.str().size(), '\b');
		message.str("");
		message << "  computing vector (lambda = " << i_lambda << ")  " << ++i_update_vec << " / " << lambda_vec.size() << "             ";

		std::cout << message.str();
		
		Parameters new_numeric_pars;
		frozen_derivative_pars = failed_jac_pars;
		for (int i_freeze=1; ; ++i_freeze)
		{
			Parameters new_frozen_derivative_pars;
			Parameters frozen_numeric_pars = par_transform.derivative2numeric_cp(frozen_derivative_pars);
			// need to remove parameters frozen due to failed jacobian runs when calling calc_lambda_upgrade_vec
			ml_upgrade = calc_lambda_upgrade_vec(jacobian, Q_sqrt, residuals_vec, par_name_vec, obs_names_vec,
				base_run_numeric_pars, frozen_numeric_pars, tot_sing_val, i_lambda);
			new_numeric_pars = apply_upgrade(base_run_numeric_pars, ml_upgrade, 1.0);
			if (i_freeze > max_freeze_iter)
			{
				new_frozen_derivative_pars = limit_parameters_ip(base_run_numeric_pars, new_numeric_pars, limit_type, frozen_numeric_pars);
			}
			else
			{
				new_frozen_derivative_pars = limit_parameters_freeze_all_ip(base_run_numeric_pars, new_numeric_pars, frozen_numeric_pars);
			}

			if (new_frozen_derivative_pars.size() > 0 && (limit_type != LimitType::FACT || limit_type != LimitType::REL ) )
			{
				frozen_derivative_pars.insert(new_frozen_derivative_pars.begin(), new_frozen_derivative_pars.end());
			}
			else
			{
				break;
			}
			update_upgrade(ml_upgrade, base_run_numeric_pars, new_numeric_pars, frozen_numeric_pars);
			if (limit_type != LimitType::UBND && limit_type != LimitType::LBND) break;
			if (frozen_numeric_pars.size() == new_numeric_pars.size()) break; // everything is frozen
		}
		//add run to run manager
		Parameters new_model_pars = par_transform.numeric2ctl_cp(new_numeric_pars);
		// make sure these are all in bounds
		for (auto &ipar :  new_model_pars)
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
		par_transform.ctl2model_ip(new_model_pars);
		run_manager.add_run(new_model_pars);
		rot_fac_vec.push_back(0);
		rot_angle_vec.push_back(0);
		magnitude_vec.push_back(ml_upgrade.norm);
		frozen_par_vec.push_back(frozen_derivative_pars);
	}
	cout << endl;
	fout_restart << "upgrade_model_runs_built " << run_manager.get_cur_groupid() << endl;
	cout << "  performing upgrade vector runs... ";
	run_manager.run();

	// process model runs
	cout << endl;
	cout << "  testing upgrade vectors... ";
	cout << endl;
	bool best_run_updated_flag = false;
	ModelRun best_upgrade_run(cur_solution);

	for(int i=0; i<run_manager.get_nruns(); ++i) {
		double rot_fac = rot_fac_vec[i];
		ModelRun upgrade_run(cur_solution);
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(i, tmp_pars, tmp_obs);
		if (success)
		{
			par_transform.model2ctl_ip(tmp_pars);
			upgrade_run.update_ctl(tmp_pars, tmp_obs);
			upgrade_run.set_frozen_ctl_parameters(frozen_par_vec[i]);
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
	const Parameters &frz_ctl_pars = best_upgrade_run.get_frozen_ctl_pars();
		
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


void SVDSolver::check_limits(const Parameters &init_derivative_pars, const Parameters &upgrade_derivative_pars,
						map<string, LimitType> &limit_type_map, Parameters &derivative_parameters_at_limit)
{
	const string *name;
	double p_init;
	double p_upgrade;
	double b_facorg_lim;
	pair<bool, double> par_limit;
	const ParameterRec *p_info;

	for(auto &ipar : init_derivative_pars)
	{
		par_limit = pair<bool, double>(false, 0.0);
		name = &(ipar.first);  // parameter name
		p_init = ipar.second; // inital parameter value
		auto pu_iter = upgrade_derivative_pars.find(*name);
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
				limit_type_map[*name] = LimitType::UBND;
		}
		// Check parameter lower bound
		else if((!par_limit.first && p_upgrade < p_info->lbnd) || 
			(par_limit.first && par_limit.second < p_info->lbnd)) {
				par_limit.first = true;
				par_limit.second = p_info->lbnd;
				limit_type_map[*name] = LimitType::LBND;
		}
		// Add any limited parameters to model_parameters_at_limit
		if (par_limit.first) {
			derivative_parameters_at_limit.insert(*name, par_limit.second);
		}
	}
}


Parameters SVDSolver::limit_parameters_ip(const Parameters &init_numeric_pars, Parameters &upgrade_numeric_pars, 
										  LimitType &limit_type, const Parameters &frozen_pars)
{
	map<string, LimitType> limit_type_map;
	const string *name;
	double p_init;
	double p_upgrade;
	double p_limit;
	Parameters bound_derivative_parameters;
	pair<bool, double> par_limit;
	Parameters derivative_parameters_at_limit;
	Parameters upgrade_derivative_pars  = par_transform.numeric2derivative_cp(upgrade_numeric_pars);

	Parameters init_derivative_pars = init_numeric_pars;
	//remove frozen parameters from numeric parameters
	init_derivative_pars.erase(frozen_pars);
	par_transform.numeric2derivative_ip(init_derivative_pars);

	check_limits(init_derivative_pars, upgrade_derivative_pars, limit_type_map, derivative_parameters_at_limit);

	// Calculate most stringent limit factor on a PEST parameter
	double limit_factor= 1.0;
	double tmp_limit;
	string limit_parameter_name = "";
	Parameters numeric_parameters_at_limit = par_transform.derivative2numeric_cp(derivative_parameters_at_limit);
	for(auto &ipar : numeric_parameters_at_limit)
	{
		name = &(ipar.first);
		p_limit = ipar.second;
		p_init = init_numeric_pars.get_rec(*name);
		p_upgrade = upgrade_numeric_pars.get_rec(*name);
		tmp_limit = (p_limit - p_init) / (p_upgrade - p_init);
		if (tmp_limit < limit_factor)  {
			limit_factor = tmp_limit;
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
	bound_derivative_parameters.clear();
	iter = numeric_parameters_at_limit.find(limit_parameter_name);
	if (iter != numeric_parameters_at_limit.end()) {
		bound_derivative_parameters[limit_parameter_name] = numeric_parameters_at_limit[limit_parameter_name];
		auto iter = limit_type_map.find(limit_parameter_name);
		assert(iter != limit_type_map.end());
		limit_type = iter->second;
	}
	par_transform.numeric2derivative_ip(bound_derivative_parameters);
	return bound_derivative_parameters;
}

Parameters SVDSolver::limit_parameters_freeze_all_ip(const Parameters &init_numeric_pars,
					Parameters &upgrade_numeric_pars, const Parameters &frozen_pars)
{
	map<string, LimitType> limit_type_map;
	const string *name;
	double p_init;
	double p_upgrade;
	double p_limit;
	Parameters frozen_derivative_parameters;
	pair<bool, double> par_limit;
	Parameters derivative_parameters_at_limit;
	Parameters upgrade_derivative_pars  = par_transform.numeric2derivative_cp(upgrade_numeric_pars);

	Parameters init_derivative_pars = init_numeric_pars;
	//remove frozen parameters from numeric parameters
	init_derivative_pars.erase(frozen_pars);
	par_transform.numeric2derivative_ip(init_derivative_pars);

	check_limits(init_derivative_pars, upgrade_derivative_pars, limit_type_map, derivative_parameters_at_limit);
	// Remove parameters at there upper and lower bound limits as these will be frozen
	vector<string> pars_at_bnds;
	for (auto ipar : limit_type_map)
	{
		if (ipar.second == LimitType::LBND || ipar.second == LimitType::UBND)
		{
			pars_at_bnds.push_back(ipar.first);
		}
	}
	derivative_parameters_at_limit.erase(pars_at_bnds);

	// Calculate most stringent limit factor on a PEST parameter
	double limit_factor = 1.0;
	double tmp_limit;
	string limit_parameter_name = "";
	Parameters numeric_parameters_at_limit = par_transform.derivative2numeric_cp(derivative_parameters_at_limit);
	for(auto &ipar : numeric_parameters_at_limit)
	{
		name = &(ipar.first);
		p_limit = ipar.second;
		p_init = init_numeric_pars.get_rec(*name);
		p_upgrade = upgrade_numeric_pars.get_rec(*name);
		tmp_limit = (p_limit - p_init) / (p_upgrade - p_init);
		if (tmp_limit < limit_factor)  {
			limit_factor = tmp_limit;
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

	//Transform parameters back their derivitive state and freeze any that violate their bounds
	upgrade_derivative_pars = par_transform.numeric2derivative_cp(upgrade_numeric_pars);
	check_limits(init_derivative_pars, upgrade_derivative_pars, limit_type_map, derivative_parameters_at_limit);
	for (auto &ipar : derivative_parameters_at_limit)
	{
		if(limit_type_map[ipar.first] == LimitType::UBND || limit_type_map[ipar.first] == LimitType::LBND)
		{
			upgrade_derivative_pars[ipar.first] = ipar.second;
			frozen_derivative_parameters.insert(ipar);
		}
	}
	//Transform parameters back their numeric state
	upgrade_numeric_pars = par_transform.derivative2numeric_cp(upgrade_derivative_pars);

	// Impose frozen Parameters
	for (auto &ipar : frozen_pars)
	{
		auto iter = upgrade_numeric_pars.find(ipar.first);
		if (iter != upgrade_numeric_pars.end())
		{
			upgrade_numeric_pars[ipar.first] = ipar.second;
		}
	}
	return frozen_derivative_parameters;
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
	double max_fac_change = 0;
	double max_rel_change = 0;
	const string *max_fac_par = 0;
	const string *max_rel_par = 0;
	const Parameters &old_ctl_pars = cur_solution.get_ctl_pars();
	const Parameters &new_ctl_pars = upgrade.get_ctl_pars();

	os << "    Parameter Upgrades (Control File Parameters)" << endl;
	os << "     Parameter      Current       Previous       Factor       Relative" << endl;
	os << "        Name         Value         Value         Change        Change" << endl;
	os << "    ------------  ------------  ------------  ------------  ------------" << endl;

	for( const auto &ipar : new_ctl_pars)
	{
		p_name = &(ipar.first);
		p_new =  ipar.second;
		p_old = old_ctl_pars.get_rec(*p_name);
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
		//os << left;
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
	os << endl;
	os << "   Maximum changes in transformed numeric parameters:" << endl;
	os << "     Maximum relative change = " << max_rel_change << "   [" << *max_rel_par << "]" << endl;
	os << "     Maximum factor change = " << max_fac_change << "   [" << *max_fac_par << "]" << endl;

	max_fac_change = 0;
	max_rel_change = 0;
	max_fac_par = 0;
	max_rel_par = 0;
	os << "    Parameter Upgrades (Transformed Numeric Parameters)" << endl;
	os << "     Parameter      Current       Previous       Factor       Relative" << endl;
	os << "        Name         Value         Value         Change        Change" << endl;
	os << "    ------------  ------------  ------------  ------------  ------------" << endl;

	const Parameters old_numeric_pars = par_transform.ctl2numeric_cp(cur_solution.get_ctl_pars());
	const Parameters new_numeric_pars = par_transform.ctl2numeric_cp(upgrade.get_ctl_pars());
	for( const auto &ipar : new_numeric_pars)
	{
		p_name = &(ipar.first);
		p_new =  ipar.second;
		p_old = old_numeric_pars.get_rec(*p_name);
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
