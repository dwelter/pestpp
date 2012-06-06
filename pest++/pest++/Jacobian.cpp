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
#include <lapackpp.h>
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

Jacobian::Jacobian(FileManager &_file_manager) : file_manager(_file_manager) 
{
}

Jacobian::~Jacobian() {
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


LaGenMatDouble Jacobian::get_matrix(const vector<string> & par_names, const vector<string> &obs_names) const
{
	int n_rows = obs_names.size();
	int n_cols = par_names.size();
	int irow_old;
	int icol_old;
	int irow_new;
	int icol_new;

	LaGenMatDouble new_matrix(n_rows, n_cols);
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
	new_matrix = 0.0;  // initialize all entries to 0
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


void Jacobian::calculate(ModelRun &init_model_run, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, bool phiredswh_flag, bool calc_init_obs)
{
	calculate(init_model_run, init_model_run.get_numeric_pars().get_keys(),  init_model_run.get_obs_template().get_keys(),
		group_info, ctl_par_info, run_manager, prior_info, phiredswh_flag, calc_init_obs);

}

void Jacobian::calculate(ModelRun &init_model_run, vector<string> numeric_par_names, vector<string> obs_names,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, bool phiredswh_flag, bool calc_init_obs)
{
	int i_run;
	string *par_name;
	Parameters model_parameters = init_model_run.get_model_pars();
	Parameters numeric_parameters = init_model_run.get_numeric_pars();
	Observations observations = init_model_run.get_obs_template();


	// compute runs for to jacobain calculation as it is influenced by derivative type( forward or central)
	vector<JacobianRun> del_numeric_par_vec;
	if (calc_init_obs) {
		del_numeric_par_vec.push_back(JacobianRun("", 0));
	}
	for(int b=0, e=numeric_par_names.size(); b!=e; ++b)
	{
		get_derivative_parameters(numeric_par_names[b], init_model_run, group_info, ctl_par_info, del_numeric_par_vec, phiredswh_flag);
	}


	// calculate array sizes and number of model runs
	int nrun = del_numeric_par_vec.size();
	if (calc_init_obs) ++nrun;

	// allocate arrays
	run_manager.allocate_memory(model_parameters, observations, nrun);
	

	// fill arrays
	for(int i=0, e=del_numeric_par_vec.size(); i!=e; ++i)
	{
		par_name = &del_numeric_par_vec[i].par_name;
		if (*par_name != "") {
			numeric_parameters = init_model_run.get_numeric_pars();
			numeric_parameters.update_rec(*par_name, del_numeric_par_vec[i].numeric_value);
			model_parameters = init_model_run.get_par_tran().numeric2model_cp(numeric_parameters);
		}
		else {
			model_parameters = init_model_run.get_model_pars();
		}
		run_manager.add_run(model_parameters);
	}
	
	// make model runs
	run_manager.run();

	// calculate jacobian
	base_numeric_par_names = numeric_par_names;
	base_sim_obs_names = obs_names;
	if(matrix.size(0) != base_sim_obs_names.size() && matrix.size(1) !=base_numeric_par_names.size())
	{
		matrix.resize(base_sim_obs_names.size(), base_numeric_par_names.size());
	}
	// initialize prior information
	prior_info_sen.clear();

	unordered_map<string, int> par2col_map = get_par2col_map();
	unordered_map<string, int>::iterator found;
	unordered_map<string, int>::iterator not_found = par2col_map.end();

	i_run = 0;
	// if initial run was was get the newly calculated values

	if (calc_init_obs) {
		run_manager.get_run(init_model_run, i_run);
		++i_run;
	}

	// process the parameter pertubation runs
	int nruns = del_numeric_par_vec.size();
	int icol;
	list<ModelRun> run_list;
	for(; i_run<nruns; ++i_run)
	{
		par_name = &del_numeric_par_vec[i_run].par_name;
		ModelRun tmp_model_run = init_model_run;
		numeric_parameters = init_model_run.get_numeric_pars();
		numeric_parameters.update_rec(*par_name, del_numeric_par_vec[i_run].numeric_value);
		tmp_model_run.set_numeric_parameters(numeric_parameters);
		run_manager.get_run(tmp_model_run, i_run);
		run_list.push_back(tmp_model_run);

		if(i_run+1>=nruns || *par_name !=  del_numeric_par_vec[i_run+1].par_name)
		{
			run_list.push_front(init_model_run);
			found = par2col_map.find(*par_name);
			assert(found != not_found);
			icol = (*found).second;
			calc_derivative(*par_name, icol, run_list, group_info, ctl_par_info, prior_info);
		}
	}
	// clean up
	run_manager.free_memory();
}

bool Jacobian::get_derivative_parameters(const string &par_name, ModelRun &init_model_run, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		vector<JacobianRun> &del_numeric_par_vec, bool phiredswh_flag)
{
	bool success = false;
	double par_0;
	double par_1;
	const ParameterGroupRec *g_rec;

	par_0 =  init_model_run.get_numeric_pars().find(par_name)->second;
	g_rec = group_info.get_group_rec_ptr(par_name);

	if (g_rec->forcen == "ALWAYS_3" || phiredswh_flag == true) {
		// Central Difference
		vector<double> new_par_vec;
		vector<Parameters> dir_numeric_pars_vec;
		success = central_diff(par_name, init_model_run.get_numeric_pars(), group_info, ctl_par_info, init_model_run.get_par_tran(), new_par_vec, dir_numeric_pars_vec);
		if (success)
		{
			del_numeric_par_vec.push_back(JacobianRun(par_name, new_par_vec[0]));
		}
	}
	if (!success) {
		// Forward Difference
		Parameters numeric_dir_pars;
		success = forward_diff(par_name, init_model_run.get_numeric_pars(), group_info, ctl_par_info, init_model_run.get_par_tran(), par_1, numeric_dir_pars);
		del_numeric_par_vec.push_back(JacobianRun(par_name, par_1));
	}
	return success;
}


void Jacobian::calc_derivative(const string &numeric_par_name, int jcol, list<ModelRun> &run_list,  const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const PriorInformation &prior_info)
{
	const string *obs_name;
	const ParameterGroupRec *g_rec;
	double del_par;
	double del_obs;
	double der;
	int irow;
	Parameters::const_iterator par_iter;
	LaGenMatDouble a_mat(3,3);
	LaVectorDouble c(3), y(3);

	// sort run_list the parameter numeric_par_name;
	ModelRunAbstractBase::Compare compare(numeric_par_name, ModelRunAbstractBase::Compare::NUMERIC_PAR);
	run_list.sort(compare);
	//p_rec = group_info.get_parameter_rec_ptr(*par_name);
	g_rec = group_info.get_group_rec_ptr(numeric_par_name);


	ModelRun &run_first = run_list.front();
	ModelRun &run_last = run_list.back();

	if (run_list.size() ==3 && g_rec->dermthd == "PARABOLIC")
	{
		// Central Difference Parabola
		// Solve Ac = o for c to get the equation for a parabola where:
		//        | p0**2  p0  1 |               | c0 |            | o0 |
		//   A =  | p1**2  p1  1 |          c =  | c1 |        y = | o1 |
		//        | p2**2  p2  1 |               | c2 |            | o2 |
		// then compute the derivative as:
		//   dy/dx = 2 * c0 * p1 + c1
		ModelRun &run2 = (*(++run_list.begin()));
		LaGenMatDouble a_mat(3,3);
		LaVectorDouble c(3), y(3);
		irow = 0;
		for (vector<string>::const_iterator b_o =base_sim_obs_names.begin(), e_o=base_sim_obs_names.end();
			b_o!=e_o; ++b_o)
		{
			obs_name = &(*b_o);
			// assemble A matrix
			int i=0;
			for (list<ModelRun>::iterator b=run_list.begin(), e=run_list.end(); b!=e; ++b, ++i) 
			{
				// assemble A matrix
				double par_value = (*b).get_numeric_pars().get_rec(numeric_par_name);
				a_mat(i,0) = par_value*par_value; a_mat(i,1) = par_value; a_mat(i,2) = 1;
				// assemble y vector
				y(i) =  (*b).get_obs().get_rec(*obs_name);
			}
			LaLinearSolve(a_mat, c, y);
			// assume for now that derivative is to be calculated around the second parameter location
			der = 2.0 * c(0) *  run2.get_numeric_pars().get_rec(numeric_par_name) + c(1);
			matrix(irow++,jcol) = der;
		}
	}
	else 
	{
		// Forward Difference and Central Difference Outer
		ModelRun &run_first = run_list.front();
		ModelRun &run_last = run_list.back();
		irow=0;
		del_par = run_last.get_numeric_pars().get_rec(numeric_par_name) - run_first.get_numeric_pars().get_rec(numeric_par_name);
		for (vector<string>::const_iterator b_o =base_sim_obs_names.begin(), e_o=base_sim_obs_names.end();
			b_o!=e_o; ++b_o)
		{
			obs_name = &(*b_o);
			del_obs = run_last.get_obs().get_rec(*obs_name) - run_first.get_obs().get_rec(*obs_name);
			matrix(irow++,jcol) = del_obs / del_par;
		}
	}
	// Prior Information allways calculated using outer model runs even for central difference
	calc_prior_info_sen(numeric_par_name, run_first, run_last, prior_info);
}

bool Jacobian::forward_diff(const string &par_name, const Parameters &numeric_parameters, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, double &new_par, 
		Parameters &numeric_derivative_pars)
{
	bool out_of_bound_forward;
	bool out_of_bound_backward;
	vector<string> out_of_bound__forward_par_vec;
	vector<string> out_of_bound__backard_par_vec;
	string tmp_name;

	double incr = derivative_inc(par_name, group_info, ctl_par_info, numeric_parameters, false);
	// try forward derivative
	// copy numeric parameters to Model_parameters
	numeric_derivative_pars = numeric_parameters;
	// perturb numeric paramateres
	new_par = numeric_derivative_pars[par_name] += incr;
	out_of_bound_forward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_derivative_pars), group_info, ctl_par_info, out_of_bound__forward_par_vec);
	if (!out_of_bound_forward) {
		return true;
	}
	// try backward derivative if forward derivative didn't work
	numeric_derivative_pars = numeric_parameters;
	new_par = numeric_derivative_pars[par_name] -= incr;
	out_of_bound_backward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_derivative_pars), group_info, ctl_par_info, out_of_bound__backard_par_vec);
	if (!out_of_bound_backward)
	{	
		return true;
	}
	return false;
}


bool Jacobian::central_diff(const string &par_name, const Parameters &pest_parameters, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, vector<double> &new_par_vec, 
		vector<Parameters>  &numeric_dir_par_vec)
{
	vector<string> out_of_bound_par_vec;
	double new_par;
	bool out_of_bnds_forward, out_of_bnds_back, out_of_bnds;
	Parameters numeric_dir_pars;
	string tmp_name;

	double incr = derivative_inc(par_name, group_info, ctl_par_info, pest_parameters, true);
	// try backward difference
	numeric_dir_pars = pest_parameters;
	new_par = numeric_dir_pars[par_name] -= incr;
	out_of_bnds_back = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), group_info, ctl_par_info, out_of_bound_par_vec);

	if (!out_of_bnds_back) {
		new_par_vec.push_back(new_par);
		numeric_dir_par_vec.push_back(numeric_dir_pars);
	}
	// try forward derivative
	numeric_dir_pars = pest_parameters;
	new_par = numeric_dir_pars[par_name] += incr;
	out_of_bnds_forward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), group_info, ctl_par_info, out_of_bound_par_vec);
	if (!out_of_bnds_forward) {
		new_par_vec.push_back(new_par);
		numeric_dir_par_vec.push_back(numeric_dir_pars);
	}
	// if backward difference was out of bounds do a second forward derivative
	if (out_of_bnds_back) {
		numeric_dir_pars = pest_parameters;
		new_par = numeric_dir_pars[par_name] += 2.0 * incr;
		out_of_bnds = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), group_info, ctl_par_info, out_of_bound_par_vec);
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
		numeric_dir_pars = pest_parameters;
		new_par = numeric_dir_pars[par_name] -= 2.0 * incr;
		out_of_bnds = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), group_info, ctl_par_info, out_of_bound_par_vec);
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


void Jacobian::calc_prior_info_sen(const string &par_name, ModelRun &run1, ModelRun &run2, const PriorInformation &prior_info)
{
	double del_par = *run2.get_numeric_pars().get_rec_ptr(par_name) - *run1.get_numeric_pars().get_rec_ptr(par_name);
	double del_prior_info;
	const string *prior_info_name;
	Parameters run2_ctl_pars = run2.get_ctl_pars();
	Parameters run1_ctl_pars = run1.get_ctl_pars();

	for (PriorInformation::const_iterator b=prior_info.begin(), e=prior_info.end(); b!=e; ++b) {
		prior_info_name = &((*b).first);
		del_prior_info = (*b).second.calc_phi(run2_ctl_pars) - (*b).second.calc_phi(run1_ctl_pars);
		if (del_prior_info != 0) {
			prior_info_sen[*prior_info_name][par_name] = del_prior_info / del_par;
		}
	}
}

bool Jacobian::out_of_bounds(const Parameters &ctl_parameters, const ParameterGroupInfo &group_info,
	const ParameterInfo &ctl_par_info, vector<string> &out_of_bound_par_vec) const
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
			out_of_bound_par_vec.push_back(*par_name);
		}
	}
	return out_of_bounds;
}

double Jacobian::derivative_inc(const string &name, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,  const Parameters &parameters, bool central)
{
	const ParameterGroupRec *g_rec;
	double incr = 0.0;	

	//// to do add error checking
	g_rec = group_info.get_group_rec_ptr(name);
	if (g_rec->inctyp == "ABSOLUTE") {
		incr = g_rec->derinc;}
	else if (g_rec->inctyp == "RELATIVE") {
		incr =  g_rec->derinc * abs((parameters.find(name)->second));
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


void Jacobian::save(const string &filename) const
{
	ofstream fout;
	fout.open(filename.c_str(), ofstream::binary);

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
	fout.write((char*) &tmp, sizeof(tmp));
	tmp = -n_obs_and_pi;
	fout.write((char*) &tmp, sizeof(tmp));

	////write number nonzero elements in jacobian 
	n = 0;
	for(int i=0, n_sobs=matrix.rows(), n_spars=matrix.cols(); i<n_sobs; ++i){
		for(j=0;j<n_spars; ++j) {
			if(matrix(i,j)!=0) ++n;
		}
	}
	//get number of nonzero derivatives in prior information
	n += size_prior_info_sen();  
	fout.write((char*)&n, sizeof(n));

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
				fout.write((char*) &(n), sizeof(n));
				fout.write((char*) &(data), sizeof(data));
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
				fout.write((char*) &(n), sizeof(n));
				fout.write((char*) &(data), sizeof(data));
			}
		}
	}

	//save parameter names
	for(vector<string>::const_iterator b=base_numeric_par_names.begin(), e=base_numeric_par_names.end();
		b!=e; ++b) {
		string_to_fortran_char(*b, par_name, 12);
		fout.write(par_name, 12);
	}

	//save observation names (part 1 standard observations)
	for(vector<string>::const_iterator b=base_sim_obs_names.begin(), e=base_sim_obs_names.end();
		b!=e; ++b) {
		string_to_fortran_char(*b, obs_name, 20);
		fout.write(obs_name, 20);
	}
	//save observation names (part 2 prior information)
	for(map<string, map<string, double>>::const_iterator b=prior_info_sen.begin(), e=prior_info_sen.end();
		b!=e; ++b)
	{
		string_to_fortran_char((*b).first, obs_name, 20);
		fout.write(obs_name, 20);
	}
	fout.close();
}

//JacobianAnalytic::JacobianAnalytic(FileManager &_file_manager, const vector<string>  *_ctl_ordered_par_names,
//	const vector<string> *_ctl_ordered_obs_names, int _jacfile) : Jacobian(_file_manager), jacfile(_jacfile)
//{
//	ctl_ordered_par_names = _ctl_ordered_par_names;
//	ctl_ordered_obs_names = _ctl_ordered_obs_names;
//}
//
//void JacobianAnalytic::calculate(ModelRun &init_model_run,  vector<string> numeric_par_names, vector<string> obs_names, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
//	RunManager &run_manager, const PriorInformation &prior_info, bool phiredswh_flag, bool calc_init_obs)
//{
//	const string *par_name;
//	int j;
//	const ParameterRec *p_rec;
//	vector<ModelRun> model_run_vec;
//	LaGenMatDouble a_mat(3,3);
//	LaVectorDouble c(3), y(3);
//
//	// First Caluclate Numeric Derivatives for Parameters where dercom != 0
//	if (calc_init_obs) 
//	{
//		run_manager.run_model(init_model_run);
//		// save base parameters and observations locally in jacobian class
//		base_numeric_par_names = numeric_par_names;
//		base_numeric_parameters = Parameters(init_model_run.get_numeric_pars(), numeric_par_names);
//		base_sim_obs_names = obs_names;
//		base_sim_observations = Observations(init_model_run.get_obs(), obs_names);
//	}
//
//	// define temp_model_run so we don't overwrite the base model run
//	ModelRun temp_model_run = init_model_run;
//	if(matrix.size(0) != temp_model_run.get_obs().size() && matrix.size(1) !=temp_model_run.get_numeric_pars().size())
//	{
//		matrix.resize(temp_model_run.get_obs().size(), temp_model_run.get_numeric_pars().size());
//	}
//	j = 0;
//
//	// First Caluclate Numeric Derivatives for Parameters where dercom != 0
//	for(Parameters::const_iterator b=temp_model_run.get_numeric_pars().begin(), 
//		e=temp_model_run.get_numeric_pars().end(); b!=e; ++b, j++)
//	{
//		par_name = &((*b).first);
//		p_rec = ctl_par_info.get_parameter_rec_ptr(*par_name);
//		if (p_rec->dercom != 0) {
//			calc_derivative(*par_name, j, temp_model_run, group_info, ctl_par_info, run_manager, prior_info, phiredswh_flag);
//		}
//	}
//	run_manager.run_model(temp_model_run, jacfile);
//	read_analytic_file(file_manager.get_analytic_derivative_filename(), ctl_par_info);
//	init_model_run.get_par_tran().analytic_derivative_to_numeric(*this, init_model_run.get_model_pars());
//}
//
//void JacobianAnalytic::read_analytic_file(const string &filename, const ParameterInfo &ctl_par_info)
//{
//	ifstream fin;
//	string line;
//	vector<string> tokens;
//	int n_par, n_obs;
//	int lnum = 1;
//	int ipar, iobs;
//	int i;
//	double data;
//	vector<int> obs_row;
//	vector<int> par_col;
//
//	// maps parameter name to row in jacobian matrix
//	i = 0;
//	{
//		map<const string, int> par_map;
//		map<const string, int>::const_iterator it;
//		for(vector<string>::const_iterator b=base_numeric_par_names.begin(), e=base_numeric_par_names.end();
//			b!=e; ++b, ++i) {
//				par_map[(*b)] = i;
//		}
//		for(vector<string>::const_iterator b=ctl_ordered_par_names->begin(), e=ctl_ordered_par_names->end();
//			b!=e; ++b) {
//				it = par_map.find(*b);
//				par_col.push_back((*it).second);
//		}
//	}
//
//	{
//		i = 0;
//		map<const string, int> obs_map;
//		map<const string, int>::const_iterator it;
//		for(vector<string>::const_iterator b=base_sim_obs_names.begin(), e=base_sim_obs_names.end();
//		b!=e; ++b, ++i) {
//			obs_map[(*b)] = i;
//		}
//		for(vector<string>::const_iterator b=ctl_ordered_obs_names->begin(), e=ctl_ordered_obs_names->end();
//			b!=e; ++b) {
//				it = obs_map.find(*b);
//				obs_row.push_back((*it).second);
//		}
//	}
//
//
//	fin.open(filename.c_str());
//
//	if (! fin) {throw(PestFileError(filename));}
//	try {
//		getline(fin, line);
//		strip_ip(line);
//		tokens.clear();
//		tokenize(line, tokens);
//		convert_ip(tokens[0], n_par);
//		convert_ip(tokens[1], n_obs);
//		++lnum;
//		for(iobs = 0, ipar = 0; getline(fin, line); ) {
//			strip_ip(line);
//			tokens.clear();
//			tokenize(line, tokens);
//			for(vector<string>::const_iterator b=tokens.begin(), e=tokens.end();
//				b!=e; ++b) 
//			{
//				convert_ip(*b, data);
//				if (data != -1.11E33) {
//					matrix(obs_row[iobs], par_col[ipar]) = data;
//				}
//				if (++ipar >= n_par) {
//					ipar=0;
//					++iobs;
//				}
//			}
//		}
//	}
//	catch (PestConversionError &e) {
//		std::stringstream out;
//		out << "Error parsing \"" << filename << "\" on line number " << lnum << endl;
//		e.add_front(out.str());
//		e.raise();
//	}
//}
//
