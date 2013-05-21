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
#include <fstream>
#include "Jacobian.h"
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
using namespace Eigen;

Jacobian::Jacobian(FileManager &_file_manager) : file_manager(_file_manager) 
{
}

Jacobian::~Jacobian() {
}


vector<string> Jacobian::obs_and_reg_list() const
{
	vector<string> all_obs = base_sim_obs_names;
	vector<string> prior_info_obs = get_map_keys(prior_info_sen);
	all_obs.insert(all_obs.end(), prior_info_obs.begin(), prior_info_obs.end());
	return all_obs;
}

unordered_map<string, int> Jacobian::get_par2col_map() const
{
	unordered_map<string, int> par2col_map;
	int icol_old = 0;
	for (vector<string>::const_iterator b=base_numeric_par_names.begin(), e=base_numeric_par_names.end();
		b!=e; ++b, ++icol_old)
	{
		par2col_map[(*b)] = icol_old;
	}
	return par2col_map;
}
	
unordered_map<string, int> Jacobian::get_obs2row_map() const
{
	unordered_map<string, int> obs2row_map;
	int irow_old = 0;
	for (vector<string>::const_iterator b=base_sim_obs_names.begin(), e=base_sim_obs_names.end();
		b!=e; ++b, ++irow_old)
	{
		obs2row_map[(*b)] = irow_old;
	}
	return obs2row_map;
}


MatrixXd Jacobian::get_matrix(const vector<string> &obs_names, const vector<string> & par_names) const
{
	int n_rows = obs_names.size();
	int n_cols = par_names.size();
	int irow_old;
	int icol_old;
	int irow_new;
	int icol_new;

	MatrixXd new_matrix(n_rows, n_cols);
	unordered_map<string, int> par2col_map = get_par2col_map();
	unordered_map<string, int> obs2row_map = get_obs2row_map();
	unordered_map<string, int> par2col_new_map;

	// Build mapping of parameter names column number in new matrix to be returned
	icol_new = 0;
	for (vector<string>::const_iterator b=par_names.begin(), e=par_names.end();
		b!=e; ++b, ++icol_new) {
		par2col_new_map[(*b)] = icol_new;
	}

	//build jacobian
	new_matrix.setConstant(0.0);  // initialize all entries to 0
	unordered_map<string, int>::const_iterator found;
	unordered_map<string, int>::const_iterator not_found_obs2row_map = obs2row_map.end();
	unordered_map<string, int>::const_iterator not_found_par2col_map = par2col_map.end();
	map<string, map<string, double>>::const_iterator found_prior_info;
	map<string, map<string, double>>::const_iterator not_found_prior_info = prior_info_sen.end();

	irow_new=0;
	for (vector<string>::const_iterator b_obs=obs_names.begin(), e_obs=obs_names.end();
		b_obs!=e_obs; ++b_obs, ++irow_new) {
		// check to see if this is a standard observation (not prior information)
		found = obs2row_map.find(*b_obs);
		found_prior_info =  prior_info_sen.find(*b_obs);
		if (found !=  not_found_obs2row_map) {
			irow_old = (*found).second;
			icol_new=0;
			for (vector<string>::const_iterator b_par=par_names.begin(), e_par=par_names.end();
				b_par!=e_par; ++b_par, ++icol_new)
			{
				found = par2col_map.find(*b_par);
				assert(found != not_found_par2col_map);  // parameter in par_names not found in current jacobian
				icol_old = (*found).second;
				new_matrix(irow_new, icol_new) = matrix(irow_old, icol_old);
			}
		}
		// check for observation in prior information
		else if (found_prior_info != not_found_prior_info)
		{
			const string *par_name;
			unordered_map<string, int>::const_iterator found_par;
			unordered_map<string, int>::const_iterator not_found_par=par2col_new_map.end();
			for(map<string, double>::const_iterator b_pi_rec=(*found_prior_info).second.begin(), e_pi_rec=(*found_prior_info).second.end();
				b_pi_rec != e_pi_rec; ++b_pi_rec) 
			{
				par_name = &(*b_pi_rec).first;
				found_par = par2col_new_map.find(*par_name);
				if (found_par != not_found_par) {
					icol_new = (*found_par).second;
					new_matrix(irow_new, icol_new) = (*b_pi_rec).second;
				}
			}
		}
		else {
			assert(true); //observation name in obs_names not found in current jacobian
		}
	}
	return new_matrix;
}


bool Jacobian::calculate(ModelRun &init_model_run, ParamTransformSeq &par_transform, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, set<string> &out_of_bound_par, bool phiredswh_flag, bool calc_init_obs)
{
	bool success;
	vector<string> numeric_par_names = par_transform.ctl2numeric_cp(init_model_run.get_ctl_pars()).get_keys();
	success = calculate(init_model_run, numeric_par_names,  init_model_run.get_obs_template().get_keys(), par_transform,
		group_info, ctl_par_info, run_manager, prior_info, out_of_bound_par, phiredswh_flag, calc_init_obs);
	return success;

}

bool Jacobian::calculate(ModelRun &init_model_run, vector<string> numeric_par_names, vector<string> obs_names, ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, set<string> &out_of_bound_par, bool phiredswh_flag, bool calc_init_obs)
{
	int i_run;
	string *par_name;
	Observations observations(init_model_run.get_obs_template());

	run_manager.reinitialize();
	const vector<string> &par_name_vec = run_manager.get_par_name_vec();
	const vector<string> &obs_name_vec = run_manager.get_obs_name_vec();

	// compute runs for to jacobain calculation as it is influenced by derivative type( forward or central)
	vector<JacobianRun> del_numeric_par_vec;
	if (calc_init_obs) {
		del_numeric_par_vec.push_back(JacobianRun("", 0));
		Parameters model_pars = par_transform.ctl2model_cp(init_model_run.get_ctl_pars());
		run_manager.add_run(model_pars);
	}
	out_of_bound_par.clear();
	bool all_par_in_bound = true;
	Parameters numeric_pars = par_transform.ctl2numeric_cp(init_model_run.get_ctl_pars());
	for(const auto &ipar_name : numeric_par_names)
	{
		// need to optimize already computing model pars in get_derivative_parameters.  should not need to compute them again
		bool success = get_derivative_parameters(ipar_name, numeric_pars, par_transform, group_info, ctl_par_info, del_numeric_par_vec, phiredswh_flag, out_of_bound_par);
		if (success)
		{
			Parameters numeric_parameters = par_transform.ctl2numeric_cp(init_model_run.get_ctl_pars());
			numeric_parameters.update_rec(ipar_name, del_numeric_par_vec.back().numeric_value);
			Parameters model_parameters = par_transform.numeric2model_cp(numeric_parameters);
			run_manager.add_run(model_parameters);
		}
		else
		{
			all_par_in_bound = false;
		}
	}
	if (!all_par_in_bound)
	{
		return false;
	}
	
	// make model runs
	run_manager.run();

	// calculate jacobian
	base_sim_obs_names = obs_names;
	matrix.resize(base_sim_obs_names.size(), numeric_par_names.size());
	// initialize prior information
	prior_info_sen.clear();

	i_run = 0;
	// if initial run was was get the newly calculated values

	if (calc_init_obs) {
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(i_run, tmp_pars, tmp_obs);
		if (!success)
		{
			throw(PestError("Error: Base parameter run failed.  Cannot compute the Jacobian"));
		}
		par_transform.model2ctl_ip(tmp_pars);
		init_model_run.update_ctl(tmp_pars, tmp_obs);
		++i_run;
	}
	base_numeric_parameters =  par_transform.ctl2numeric_cp(init_model_run.get_ctl_pars());

	// process the parameter pertubation runs
	int nruns = del_numeric_par_vec.size();
	int icol;
	list<ModelRun> run_list;
	base_numeric_par_names.clear();
	for(; i_run<nruns; ++i_run)
	{
		par_name = &del_numeric_par_vec[i_run].par_name;
		ModelRun tmp_model_run = init_model_run;
		Parameters numeric_parameters = base_numeric_parameters;
		numeric_parameters.update_rec(*par_name, del_numeric_par_vec[i_run].numeric_value);
		tmp_model_run.set_ctl_parameters(par_transform.numeric2ctl_cp(numeric_parameters));
		Parameters tmp_pars;
		Observations tmp_obs;
		int status=run_manager.get_run(i_run, tmp_pars, tmp_obs);
		if (status <1)
		{
			throw(PestError("Error: run failed.  Cannot compute the Jacobian"));
		}
		par_transform.model2ctl_ip(tmp_pars);
		tmp_model_run.update_ctl(tmp_pars, tmp_obs);
		run_list.push_back(tmp_model_run);

		if(i_run+1>=nruns || *par_name !=  del_numeric_par_vec[i_run+1].par_name)
		{
			base_numeric_par_names.push_back(*par_name);
			run_list.push_front(init_model_run);
			icol = base_numeric_par_names.size() - 1;
			calc_derivative(*par_name, icol, run_list, par_transform, group_info, ctl_par_info, prior_info);
			run_list.clear();
		}
	}
	matrix.conservativeResize(base_sim_obs_names.size(), base_numeric_par_names.size());
	// clean up
	run_manager.free_memory();
	return true;
}

bool Jacobian::get_derivative_parameters(const string &par_name, Parameters &numeric_pars, ParamTransformSeq &par_transform, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		vector<JacobianRun> &del_numeric_par_vec, bool phiredswh_flag, set<string> &out_of_bound_par)
{
	bool success = false;
	double par_0;
	double par_1;
	const ParameterGroupRec *g_rec;

	par_0 =  numeric_pars.find(par_name)->second;
	g_rec = group_info.get_group_rec_ptr(par_name);

	if (g_rec->forcen == "ALWAYS_3" || phiredswh_flag == true) {
		// Central Difference
		vector<double> new_par_vec;
		vector<Parameters> dir_numeric_pars_vec;
		success = central_diff(par_name, numeric_pars, group_info, ctl_par_info, par_transform, new_par_vec, dir_numeric_pars_vec, out_of_bound_par);
		if (success)
		{
			del_numeric_par_vec.push_back(JacobianRun(par_name, new_par_vec[0]));
		}
	}
	if (!success) {
		// Forward Difference
		success = forward_diff(par_name, numeric_pars, group_info, ctl_par_info, par_transform, par_1, out_of_bound_par);
		if (success)
		{
			del_numeric_par_vec.push_back(JacobianRun(par_name, par_1));
		}
	}
	return success;
}


void Jacobian::calc_derivative(const string &numeric_par_name, int jcol, list<ModelRun> &run_list,  const ParamTransformSeq &par_trans, 
							   const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const PriorInformation &prior_info)
{
	const ParameterGroupRec *g_rec;
	double del_par;
	double del_obs;
	double der;
	int irow;
	Parameters::const_iterator par_iter;
	list<pair<ModelRun &, Parameters>> numeric_run_list;;
	for(auto &irun : run_list)
	{
		numeric_run_list.push_back(pair<ModelRun &, Parameters>(irun, par_trans.ctl2numeric_cp(irun.get_ctl_pars())));
	}
	// sort run_list the parameter numeric_par_name;
	auto compare = [&numeric_par_name] (const pair<ModelRun &, Parameters> &a, const pair<ModelRun &, Parameters> &b)
		{return a.second.get_rec(numeric_par_name) > b.second.get_rec(numeric_par_name);};
	numeric_run_list.sort(compare);
	auto &run_first = numeric_run_list.front();
	auto  &run_last = numeric_run_list.back();

	//p_rec = group_info.get_parameter_rec_ptr(*par_name);
	g_rec = group_info.get_group_rec_ptr(numeric_par_name);

	if (numeric_run_list.size() ==3 && g_rec->dermthd == "PARABOLIC")
	{
		// Central Difference Parabola
		// Solve Ac = o for c to get the equation for a parabola where:
		//        | p0**2  p0  1 |               | c0 |            | o0 |
		//   A =  | p1**2  p1  1 |          c =  | c1 |        y = | o1 |
		//        | p2**2  p2  1 |               | c2 |            | o2 |
		// then compute the derivative as:
		//   dy/dx = 2 * c0 * p1 + c1
		auto &run2 = (*(++numeric_run_list.begin()));
		MatrixXd a_mat(3,3);
		VectorXd c(3), y(3);
		irow = 0;
		for (auto &iobs_name : base_sim_obs_names)
		{
			// assemble A matrix
			int i=0;
			for (auto &irun_pair : numeric_run_list)
				//list<ModelRun>::iterator b=run_list.begin(), e=run_list.end(); b!=e; ++b, ++i) 
			{
				// assemble A matrix
				double par_value = irun_pair.second.get_rec(numeric_par_name);
				a_mat(i,0) = par_value*par_value; a_mat(i,1) = par_value; a_mat(i,2) = 1;
				// assemble y vector
				y(i) =  irun_pair.first.get_obs().get_rec(iobs_name);
			}
			c = a_mat.colPivHouseholderQr().solve(y);
			// assume for now that derivative is to be calculated around the second parameter location
			der = 2.0 * c(0) *  run2.second.get_rec(numeric_par_name) + c(1);
			matrix(irow++,jcol) = der;
		}
	}
	else 
	{
		// Forward Difference and Central Difference Outer
		irow=0;
		del_par = run_last.second.get_rec(numeric_par_name) - run_first.second.get_rec(numeric_par_name);
		for (auto &iobs_name : base_sim_obs_names)
		{
			del_obs = run_last.first.get_obs().get_rec(iobs_name) - run_first.first.get_obs().get_rec(iobs_name);
			matrix(irow++,jcol) = del_obs / del_par;
		}
	}
	// Prior Information allways calculated using outer model runs even for central difference
	del_par = run_last.second.get_rec(numeric_par_name) - run_first.second.get_rec(numeric_par_name);
	calc_prior_info_sen(numeric_par_name, run_first.first.get_ctl_pars(), run_last.first.get_ctl_pars(), del_par,  prior_info);
}


bool Jacobian::forward_diff(const string &par_name, const Parameters &numeric_parameters, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		const ParamTransformSeq &par_trans, double &new_par, set<string> &out_of_bound_par)
{
	// if the transformation is one to one, call the simpler and more effiecient version of this routine
	// designed specifically for this case
	bool out_of_bound_forward;
	bool out_of_bound_backward;
	string tmp_name;

	double incr = derivative_inc(par_name, group_info, ctl_par_info, numeric_parameters.get_rec(par_name), false);
	// try forward derivative
	// perturb numeric paramateres
	Parameters numeric_derivative_pars(numeric_parameters);
	new_par = numeric_derivative_pars[par_name] += incr;
	out_of_bound_forward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_derivative_pars), group_info, ctl_par_info, out_of_bound_par);
	if (!out_of_bound_forward) {
		return true;
	}
	// try backward derivative if forward derivative didn't work
	numeric_derivative_pars = numeric_parameters;
	new_par = numeric_derivative_pars[par_name] -= incr;
	out_of_bound_backward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_derivative_pars), group_info, ctl_par_info, out_of_bound_par);
	if (!out_of_bound_backward)
	{	
		return true;
	}
	return false;
}


bool Jacobian::central_diff(const string &par_name, const Parameters &pest_parameters, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, vector<double> &new_par_vec, 
		vector<Parameters>  &numeric_dir_par_vec, set<string> &out_of_bound_par)
{
	double new_par;
	bool out_of_bnds_forward, out_of_bnds_back, out_of_bnds;
	Parameters numeric_dir_pars;
	string tmp_name;

	double incr = derivative_inc(par_name, group_info, ctl_par_info, pest_parameters.get_rec(par_name), true);
	// try backward difference
	numeric_dir_pars = pest_parameters;
	new_par = numeric_dir_pars[par_name] -= incr;
	out_of_bnds_back = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), group_info, ctl_par_info, out_of_bound_par);

	if (!out_of_bnds_back) {
		new_par_vec.push_back(new_par);
		numeric_dir_par_vec.push_back(numeric_dir_pars);
	}
	// try forward derivative
	numeric_dir_pars = pest_parameters;
	new_par = numeric_dir_pars[par_name] += incr;
	out_of_bnds_forward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), group_info, ctl_par_info, out_of_bound_par);
	if (!out_of_bnds_forward) {
		new_par_vec.push_back(new_par);
		numeric_dir_par_vec.push_back(numeric_dir_pars);
	}
	// if backward difference was out of bounds do a second forward derivative
	if (out_of_bnds_back) {
		set<string> tmp_out_of_bound_par;
		numeric_dir_pars = pest_parameters;
		new_par = numeric_dir_pars[par_name] += 2.0 * incr;
		out_of_bnds = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), group_info, ctl_par_info, tmp_out_of_bound_par);
		if (!out_of_bnds) {
			new_par_vec.push_back(new_par);
			numeric_dir_par_vec.push_back(numeric_dir_pars);
		}
		else
		{
			return false;  // can't do central difference without going out of bounds
		}
	}
	// if forward difference was out of bounds do a second backward derivative
	if (out_of_bnds_forward) {
		set<string> tmp_out_of_bound_par;
		numeric_dir_pars = pest_parameters;
		new_par = numeric_dir_pars[par_name] -= 2.0 * incr;
		out_of_bnds = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), group_info, ctl_par_info, tmp_out_of_bound_par);
		if (!out_of_bnds) {
			new_par_vec.insert(new_par_vec.begin(), new_par);
			numeric_dir_par_vec.insert(numeric_dir_par_vec.begin(), numeric_dir_pars);
		}
		else
		{
			return false;  // can't do central difference without going out of bounds
		}
	}
	return true;
}


void Jacobian::calc_prior_info_sen(const string &par_name, const Parameters &ctl_pars_1, const Parameters &ctl_pars_2, double del_numeric_par,  const PriorInformation &prior_info)
{
	double del_prior_info;
	const string *prior_info_name;
	const PriorInformationRec *pi_rec;

	for (auto &i_pinfo : prior_info)
	{
		prior_info_name = &(i_pinfo.first);
		pi_rec = &(i_pinfo.second);
		del_prior_info = i_pinfo.second.calc_phi(ctl_pars_2) - i_pinfo.second.calc_phi(ctl_pars_1);
		if (del_prior_info != 0) {
			prior_info_sen[*prior_info_name][par_name] =
			(pi_rec->calc_residual(ctl_pars_2) - pi_rec->calc_residual(ctl_pars_1)) / del_numeric_par;
		}
	}
}

bool Jacobian::out_of_bounds(const Parameters &ctl_parameters, const ParameterGroupInfo &group_info,
	const ParameterInfo &ctl_par_info, set<string> &out_of_bound_par) const
{
	const string *par_name;
	double min, max;
	const ParameterRec *par_info_ptr;
	bool out_of_bounds=false;

	for(Parameters::const_iterator b=ctl_parameters.begin(), e=ctl_parameters.end();
		b!=e; ++b) {
			par_name = &(*b).first;
			par_info_ptr = ctl_par_info.get_parameter_rec_ptr(*par_name);
			max = par_info_ptr->ubnd;
			min = par_info_ptr->lbnd;
		if ((*b).second > max || (*b).second < min) {
			out_of_bounds = true;
			out_of_bound_par.insert(*par_name);
		}
	}
	return out_of_bounds;
}

double Jacobian::derivative_inc(const string &name, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, double cur_par_value, bool central)
{
	const ParameterGroupRec *g_rec;
	double incr = 0.0;	

	//// to do add error checking
	g_rec = group_info.get_group_rec_ptr(name);
	if (g_rec->inctyp == "ABSOLUTE") {
		incr = g_rec->derinc;}
	else if (g_rec->inctyp == "RELATIVE") {
		incr =  g_rec->derinc * abs(cur_par_value);
	}
	// apply derincmul for central derivatives
	if (central) {
		incr *= g_rec->derincmul;
	}
	// apply lower bound
	if (incr < g_rec->derinclb) {
		incr = g_rec->derinclb;
	}
	return incr;
}

int Jacobian::size_prior_info_sen() const
{
	int n_sen = 0;
	for(map<string, map<string, double>>::const_iterator b=prior_info_sen.begin(), e=prior_info_sen.end();
		b!=e; ++b)
	{
		n_sen += (*b).second.size();
	}
	return n_sen;
}
const set<string>& Jacobian::get_failed_parameter_names() const
{
	return failed_parameter_names;
}

void Jacobian::save(const string &ext) const
{
	ofstream &jout = file_manager.open_ofile_ext(ext, ios::out |ios::binary);

	int n_par = base_numeric_par_names.size();
	int n_standard_obs = base_sim_obs_names.size();
	int n_obs_and_pi = n_standard_obs + prior_info_sen.size();
	int i,j,n;
	int tmp;
	double data;
	const string *par_name_ptr;
	char par_name[12];
	char obs_name[20];

	// write header
	tmp  = -n_par;
	jout.write((char*) &tmp, sizeof(tmp));
	tmp = -n_obs_and_pi;
	jout.write((char*) &tmp, sizeof(tmp));

	////write number nonzero elements in jacobian 
	n = 0;
	for(int i=0, n_sobs=matrix.rows(), n_spars=matrix.cols(); i<n_sobs; ++i){
		for(j=0;j<n_spars; ++j) {
			if(matrix(i,j)!=0) ++n;
		}
	}
	//get number of nonzero derivatives in prior information
	n += size_prior_info_sen();  
	jout.write((char*)&n, sizeof(n));

	////write matrix
	n = 0;
	map<string, double>::const_iterator found_pi_par;
	map<string, double>::const_iterator not_found_pi_par;
	for(j=0;j<n_par; ++j) {
		par_name_ptr = &base_numeric_par_names[j];
		// standard observations
		for(i=0; i<n_standard_obs; ++i){
			data = matrix(i,j);
			++n;
			if (data != 0.0) {
				jout.write((char*) &(n), sizeof(n));
				jout.write((char*) &(data), sizeof(data));
			}
		}
		//prior information

		for(map<string, map<string, double>>::const_iterator b=prior_info_sen.begin(), e=prior_info_sen.end();
		b!=e; ++b)
		{
			++n;
			not_found_pi_par = (*b).second.end();
			found_pi_par = (*b).second.find(*par_name_ptr);
			if (found_pi_par != not_found_pi_par) {
				data = (*found_pi_par).second;
				jout.write((char*) &(n), sizeof(n));
				jout.write((char*) &(data), sizeof(data));
			}
		}
	}

	//save parameter names
	for(vector<string>::const_iterator b=base_numeric_par_names.begin(), e=base_numeric_par_names.end();
		b!=e; ++b) {
		string_to_fortran_char(*b, par_name, 12);
		jout.write(par_name, 12);
	}

	//save observation names (part 1 standard observations)
	for(vector<string>::const_iterator b=base_sim_obs_names.begin(), e=base_sim_obs_names.end();
		b!=e; ++b) {
		string_to_fortran_char(*b, obs_name, 20);
		jout.write(obs_name, 20);
	}
	//save observation names (part 2 prior information)
	for(map<string, map<string, double>>::const_iterator b=prior_info_sen.begin(), e=prior_info_sen.end();
		b!=e; ++b)
	{
		string_to_fortran_char((*b).first, obs_name, 20);
		jout.write(obs_name, 20);
	}
	file_manager.close_file(ext);
}


//void Jacobian::read(const string &filename)
//{
//	ifstream fin;
//	fin.open(filename.c_str(), ifstream::binary);
//
//	int n_par = base_numeric_par_names.size();
//	int n_standard_obs = base_sim_obs_names.size();
//	int n_nonzero;
//	int n_obs_and_pi;
//	int i,j,n;
//	double data;
//	char par_name[12];
//	char obs_name[20];
//
//	// read header
//	fin.read((char*) &n_par, sizeof(n_par));
//	fin.read((char*) &n_obs_and_pi, sizeof(n_obs_and_pi));
//	n_par = -n_par;
//	n_obs_and_pi = -n_obs_and_pi;
//
//	if (n_obs_and_pi == base_sim_obs_names.size() + prior_info_sen.size())
//	{
//		cerr << "Error Reading Jacobian: Prior number of observations and prior information records in current problem and file are inconsistent" << endl;
//		throw(PestError("Error Reading Jacobian: Prior number of observations and prior information records in current problem and file are inconsistent"));
//	}
//
//	////read number nonzero elements in jacobian (observations + prior information)
//	fin.read((char*)&n_nonzero, sizeof(n_nonzero));
//	// read matrix
//	LaGenMatDouble tmp_matrix = LaGenMatDouble::zeros(n_obs_and_pi, n_par);
//	for (int i_rec=0; i_rec<n_nonzero; ++ i_rec)
//	{
//		fin.read((char*) &(n), sizeof(n));
//		fin.read((char*) &(data), sizeof(data));
//		j = int(n/n_par);
//		i = n % n_par;
//		tmp_matrix(i,j) = data;
//	}
//	//read parameter names
//	vector<string> tmp_par_names;
//	for (int i_rec=0; i_rec<n_par; ++i_rec)
//	{
//		fin.read(par_name, 12);
//		string temp_par = string(par_name, 12);
//trouble here		strip_ip(temp_par);
//		tmp_par_names.push_back(temp_par);
//	}
//	
//	//read observation and Prior info names
//	vector<string> tmp_obs_pi_names;
//	for (int i_rec=0; i_rec<n_obs_and_pi; ++i_rec)
//	{
//		fin.read(obs_name, 20);
//		tmp_obs_pi_names.push_back(strip_cp(string(obs_name, 20)));
//	}
//		//prior information
//
//	//	prior_info_sen.clear();
//	//	for (int i_rec=n_standard_obs; i_rec<n_obs_and_pi; ++i_rec)
//	//	{
//
//	//	}
//
//	//	for(map<string, map<string, double>>::const_iterator b=prior_info_sen.begin(), e=prior_info_sen.end();
//	//	b!=e; ++b)
//	//	{
//	//		++n;
//	//		not_found_pi_par = (*b).second.end();
//	//		found_pi_par = (*b).second.find(*par_name_ptr);
//	//		if (found_pi_par != not_found_pi_par) {
//	//			data = (*found_pi_par).second;
//	//			fout.write((char*) &(n), sizeof(n));
//	//			fout.write((char*) &(data), sizeof(data));
//	//		}
//	//	}
////	}
//	fin.close();
//}

