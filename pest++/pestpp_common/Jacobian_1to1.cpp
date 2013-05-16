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

#include <cstdlib>
#include <vector>
#include <set>
#include <fstream>
#include <algorithm>
#include "Jacobian_1to1.h"
#include "Transformable.h"
#include "ParamTransformSeq.h"
#include "pest_error.h"
#include "pest_data_structs.h"
#include "ModelRunPP.h"
#include "RunManagerAbstract.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "PriorInformation.h"

using namespace std;
using namespace pest_utils;

Jacobian_1to1::Jacobian_1to1(FileManager &_file_manager) : Jacobian(_file_manager) 
{
}

Jacobian_1to1::~Jacobian_1to1() {
}

bool Jacobian_1to1::calculate(ModelRun &init_model_run, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, set<string> &out_of_bound_par, bool phiredswh_flag, bool calc_init_obs)
{
	bool success;
	success = calculate(init_model_run, init_model_run.get_numeric_pars().get_keys(),  init_model_run.get_obs_template().get_keys(),
		group_info, ctl_par_info, run_manager, prior_info, out_of_bound_par, phiredswh_flag, calc_init_obs);
	return success;
}


bool Jacobian_1to1::calculate(ModelRun &init_model_run, vector<string> numeric_par_names, vector<string> obs_names,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, set<string> &out_of_bound_par, bool phiredswh_flag, bool calc_init_obs)
{
	int i_run;
	string *par_name;
	Parameters model_parameters(init_model_run.get_model_pars());
	Observations observations(init_model_run.get_obs_template());
	base_numeric_parameters = init_model_run.get_numeric_pars();
	run_manager.reinitialize(file_manager.build_filename("rnj"));
	const vector<string> &model_par_name_vec = run_manager.get_par_name_vec();
	const vector<string> &obs_name_vec = run_manager.get_obs_name_vec();
	size_t n_par = model_par_name_vec.size();

	failed_parameter_names.clear();
	// compute runs for to jacobain calculation as it is influenced by derivative type( forward or central)
	vector<JacobianRun> del_numeric_par_vec;
	if (calc_init_obs) {
		del_numeric_par_vec.push_back(JacobianRun("", 0));
		run_manager.add_run(model_parameters);
	}
	Parameters new_derivative_pars;
	bool success;
	Parameters base_derivative_parameters = init_model_run.get_par_tran().numeric2derivative_cp(base_numeric_parameters);
	//Loop through derivative parameters and build the parameter sets necessary for computing the jacobian
	for (auto &i_name : numeric_par_names)
	{
		assert(base_derivative_parameters.find(i_name) != base_derivative_parameters.end());
		vector<JacobianRun> tmp_del_numeric_par_vec;
		double derivative_par_value = base_derivative_parameters.get_rec(i_name);
		success = get_derivative_parameters(i_name, derivative_par_value, init_model_run.get_par_tran(), group_info, ctl_par_info,
			tmp_del_numeric_par_vec, phiredswh_flag);
		if (success && !tmp_del_numeric_par_vec.empty())
		{
			del_numeric_par_vec.insert(del_numeric_par_vec.end(), tmp_del_numeric_par_vec.begin(), tmp_del_numeric_par_vec.end());
			// update changed model parameters in model_parameters
			for (const auto &par : tmp_del_numeric_par_vec)
			{
				Parameters org_pars;
				Parameters new_pars;
				org_pars.insert(make_pair(par.par_name, model_parameters.get_rec(par.par_name)));
				new_pars.insert(make_pair(par.par_name, par.numeric_value));
				init_model_run.get_par_tran().derivative2model_ip(new_pars);
				for (auto &ipar : new_pars)
				{
					model_parameters[ipar.first] = ipar.second;
				}
				run_manager.add_run(model_parameters);
				for (const auto &ipar : org_pars)
				{
					model_parameters[ipar.first] = ipar.second;
				}
			}
		}
		else 
		{
			failed_parameter_names.insert(i_name);
		}
	}
	ofstream &fout_restart = file_manager.get_ofstream("rst");
	fout_restart << "jacobian_model_runs_begin_group_id " << run_manager.get_cur_groupid() << endl;
	// make model runs
	run_manager.run();
	fout_restart << "jacobian_model_runs_end_group_id " << run_manager.get_cur_groupid() << endl;

	// calculate jacobian
	base_sim_obs_names = obs_names;

	if(matrix.rows() != base_sim_obs_names.size() || matrix.cols() !=numeric_par_names.size())
	{
		matrix.resize(base_sim_obs_names.size(), numeric_par_names.size());
	}
	// initialize prior information
	prior_info_sen.clear();

	unordered_map<string, int> par2col_map = get_par2col_map();
	unordered_map<string, int>::iterator found;
	unordered_map<string, int>::iterator not_found = par2col_map.end();

	i_run = 0;
	// if initial run was performed, get the newly calculated values

	if (calc_init_obs) {
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(i_run, tmp_pars, tmp_obs);
		if (!success)
		{
			throw(PestError("Error: Base parameter run failed.  Can not compute the Jacobian"));
		}
		init_model_run.update(tmp_pars, tmp_obs);
		++i_run;
	}

	// process the parameter pertubation runs
	int nruns = del_numeric_par_vec.size();
	list<ModelRun> run_list;
	base_numeric_par_names.clear();
	ModelRun tmp_model_run = init_model_run;
	int icol = 0;
	for(; i_run<nruns; ++i_run)
	{
		par_name = &del_numeric_par_vec[i_run].par_name;
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(i_run, tmp_pars, tmp_obs);
		if (success)
		{
			tmp_model_run.update(tmp_pars, tmp_obs);
			run_list.push_back(tmp_model_run);
		}

		if(i_run+1>=nruns || *par_name !=  del_numeric_par_vec[i_run+1].par_name)
		{
			if (!run_list.empty())
			{
				base_numeric_par_names.push_back(*par_name);
				run_list.push_front(init_model_run);
				calc_derivative(*par_name, icol, run_list, group_info, ctl_par_info, prior_info);
				icol++;
			}
			else
			{
				failed_parameter_names.insert(*par_name);
			}
		}
	}
	matrix.conservativeResize(base_sim_obs_names.size(), base_numeric_par_names.size());
	// clean up
	fout_restart << "jacobian_saved" << endl;
	run_manager.free_memory();
	return true;
}

bool Jacobian_1to1::get_derivative_parameters(const string &par_name, double par_value, const ParamTransformSeq &par_trans, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		vector<JacobianRun> &del_numeric_par_vec, bool phiredswh_flag)
{
	bool success = false;

	double par_0 =  par_value;
	double par_1;
	const ParameterGroupRec *g_rec = group_info.get_group_rec_ptr(par_name);

	if (g_rec->forcen == "ALWAYS_3" || phiredswh_flag == true) {
		// Central Difference
		vector<double> new_par_vec;
		vector<Parameters> dir_numeric_pars_vec;
		success = central_diff(par_name, par_value, group_info, ctl_par_info, par_trans, new_par_vec, dir_numeric_pars_vec);
		if (success)
		{
			for (auto & ipar : new_par_vec)
			{
				del_numeric_par_vec.push_back(JacobianRun(par_name, ipar));
			}
		}
	}
	if (!success) {
		// Forward Difference
		success = forward_diff(par_name, par_value, group_info, ctl_par_info, par_trans, par_1);
		if(success)
		{
			del_numeric_par_vec.push_back(JacobianRun(par_name, par_1));
		}
	}
	return success;
}

bool Jacobian_1to1::forward_diff(const string &par_name, double base_derivative_val, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, double &new_par_val)
{
	const ParameterRec *par_info_ptr = ctl_par_info.get_parameter_rec_ptr(par_name);
	Parameters new_par;
	bool out_of_bound_forward;
	bool out_of_bound_backward;
	vector<string> out_of_bound__forward_par_vec;
	vector<string> out_of_bound__backard_par_vec;
	string tmp_name;

	// perturb derivative parameters
	double incr = derivative_inc(par_name, group_info, ctl_par_info, base_derivative_val, false);
	new_par_val = new_par[par_name] = base_derivative_val + incr;
	// try forward derivative
	out_of_bound_forward = out_of_bounds(par_trans.derivative2ctl_cp(new_par), group_info, ctl_par_info);
	if (!out_of_bound_forward) {
		return true;
	}
	// try backward derivative if forward derivative didn't work
	new_par.clear();
	new_par_val = new_par[par_name] = base_derivative_val - incr;
	out_of_bound_backward = out_of_bounds(par_trans.derivative2ctl_cp(new_par), group_info, ctl_par_info);
	if (!out_of_bound_backward)
	{	
		return true;
	}
	return false;
}

bool Jacobian_1to1::central_diff(const string &par_name, double base_derivative_val, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, vector<double> &new_par_vec, 
		vector<Parameters>  &perturb_derivative_par_vec)
{
	double new_par;
	bool out_of_bnds_forward, out_of_bnds_back, out_of_bnds;
	Parameters perturb_derivative_pars;
	string tmp_name;

	double incr = derivative_inc(par_name, group_info, ctl_par_info, base_derivative_val, true);
	// try backward difference
	new_par = perturb_derivative_pars[par_name] = base_derivative_val - incr;
	out_of_bnds_back = out_of_bounds(par_trans.derivative2ctl_cp(perturb_derivative_pars), group_info, ctl_par_info);

	if (!out_of_bnds_back) {
		new_par_vec.push_back(new_par);
		perturb_derivative_par_vec.push_back(perturb_derivative_pars);
	}
	// try forward derivative
	new_par = perturb_derivative_pars[par_name] = base_derivative_val + incr;
	out_of_bnds_forward = out_of_bounds(par_trans.derivative2ctl_cp(perturb_derivative_pars), group_info, ctl_par_info);
	if (!out_of_bnds_forward) {
		new_par_vec.push_back(new_par);
		perturb_derivative_par_vec.push_back(perturb_derivative_pars);
	}
	// if backward difference was out of bounds do a second forward derivative
	if (out_of_bnds_back) {
		new_par = perturb_derivative_pars[par_name] = base_derivative_val + 2.0 * incr;
		out_of_bnds = out_of_bounds(par_trans.numeric2ctl_cp(perturb_derivative_pars), group_info, ctl_par_info);
		if (!out_of_bnds) {
			new_par_vec.push_back(new_par);
			perturb_derivative_par_vec.push_back(perturb_derivative_pars);
		}
		else
		{
			return false;  // can't do central difference without going out of bounds
		}
	}
	// if forward difference was out of bounds do a second backward derivative
	if (out_of_bnds_forward) {
		new_par = perturb_derivative_pars[par_name] = base_derivative_val - 2.0 * incr;
		out_of_bnds = out_of_bounds(par_trans.numeric2ctl_cp(perturb_derivative_pars), group_info, ctl_par_info);
		if (!out_of_bnds) {
			new_par_vec.insert(new_par_vec.begin(), new_par);
			perturb_derivative_par_vec.push_back(perturb_derivative_pars);
		}
		else
		{
			return false;  // can't do central difference without going out of bounds
		}
	}
	return true;
}

bool Jacobian_1to1::out_of_bounds(const Parameters &ctl_parameters, const ParameterGroupInfo &group_info,
	const ParameterInfo &ctl_par_info) const
{
	const string *par_name;
	double min, max;
	const ParameterRec *par_info_ptr;
	bool out_of_bounds=false;

	for (auto &p : ctl_parameters)
	{
			par_name = &(p.first);
			par_info_ptr = ctl_par_info.get_parameter_rec_ptr(*par_name);
			max = par_info_ptr->ubnd;
			min = par_info_ptr->lbnd;
		if (p.second > max || p.second < min) {
			out_of_bounds = true;
		}
	}
	return out_of_bounds;
}