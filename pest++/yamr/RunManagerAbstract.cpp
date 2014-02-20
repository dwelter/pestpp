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
#include "RunManagerAbstract.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>
#include <cstring>
#include "Transformable.h"
#include "utilities.h"

RunManagerAbstract::RunManagerAbstract(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, int _max_n_failure)
	: total_runs(0), comline_vec(_comline_vec), tplfile_vec(_tplfile_vec),
	inpfile_vec(_inpfile_vec), insfile_vec(_insfile_vec), outfile_vec(_outfile_vec),
	file_stor(stor_filename), max_n_failure(_max_n_failure)
{
	cout << endl;
	cout << "             Generalized Run Manager Interface" << endl;
	cout << "                       developed by:" << endl << endl;
	cout << "                        Dave Welter" << endl;
	cout << "          Computational Water Resource Engineering" << endl;
	cout << endl << endl;
	cur_group_id = -1;
}

void RunManagerAbstract::initialize(const Parameters &model_pars, const Observations &obs, const string &_filename)
{
	file_stor.reset(model_pars.get_keys(), obs.get_keys(), _filename);
}

void RunManagerAbstract::initialize(const std::vector<std::string> &par_names, std::vector<std::string> &obs_names, const string &_filename)
{
	file_stor.reset(par_names, obs_names, _filename);
}

void RunManagerAbstract::reinitialize(const string &_filename)
{
	vector<string> par_names = get_par_name_vec();
	vector<string> obs_names = get_obs_name_vec();
	file_stor.reset(par_names, obs_names, _filename);
}

void RunManagerAbstract::initialize_restart(const std::string &_filename)
{

	file_stor.init_restart(_filename);
}

int RunManagerAbstract::add_run(const vector<double> &model_pars)
{
	int run_id = file_stor.add_run(model_pars);
	return run_id;
}

int RunManagerAbstract::add_run(const Parameters &model_pars)
{
	int run_id = file_stor.add_run(model_pars);
	return run_id;
}

int RunManagerAbstract::add_run(const Eigen::VectorXd &model_pars)
{
	int run_id = file_stor.add_run(model_pars);
	return run_id;
}

 const vector<string>& RunManagerAbstract::get_par_name_vec() const
 {
	return file_stor.get_par_name_vec();
 }

 const vector<string>& RunManagerAbstract::get_obs_name_vec() const
{
	return file_stor.get_obs_name_vec();
}

bool RunManagerAbstract::get_run(int run_id, Parameters &pars, Observations &obs)
{
	bool success = false;
	int status = file_stor.get_run(run_id, &pars, &obs);
	if (status > 0) success = true; 
	return success;
}

bool  RunManagerAbstract::get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs)
{
	bool success = false;
	int status = file_stor.get_run(run_id, pars, npars, obs, nobs);
	if (status > 0) success = true; 
	return success;
}


void  RunManagerAbstract::free_memory()
{
}

bool RunManagerAbstract::n_run_failures_exceeded(int id)
{
	bool ret_val;
	int istatus = file_stor.get_run_status(id);
	if (istatus<= -max_n_failure  &&  istatus > -100)
	 {
		 ret_val = true;
	 }
	else
	{
		ret_val = false;
	}
	return ret_val;
}

const std::set<int> RunManagerAbstract::get_failed_run_ids()
{
	std::set<int> failed_runs;
	int n_runs = file_stor.get_nruns();
	for (int id=0; id<n_runs; ++id)
	{
		if(n_run_failures_exceeded(id))
		 {
			 failed_runs.insert(id);
		 }
	}
	return failed_runs;
}

int RunManagerAbstract::get_num_failed_runs(void)
{
	int n_failed = 0;
	int n_runs = file_stor.get_nruns();
	for (int id=0; id<n_runs; ++id)
	 {
		if(n_run_failures_exceeded(id))
		{
			++n_failed;
		 }
	 }
	 return n_failed;
}

Parameters RunManagerAbstract::get_model_parameters(int run_id)
 {
	 return file_stor.get_parameters(run_id); 
 }

 Observations RunManagerAbstract::get_obs_template(double value) const
 {
	Observations ret_obs;
	const vector<string> &obs_name_vec = file_stor.get_obs_name_vec();
	int nobs = obs_name_vec.size();
	for(int i=0; i<nobs; ++i)
	{
		ret_obs[obs_name_vec[i]] = value;
	}
	return ret_obs;
 }
 int RunManagerAbstract::get_cur_groupid()
 {
	 return cur_group_id;
 }
 
 bool RunManagerAbstract::run_requried(int run_id)
 {
	 bool ret_val;
	 int istatus = file_stor.get_run_status(run_id);
	 if (istatus > -max_n_failure)
	 {
		 ret_val = true;
	 }
	 else 
	 {
		 ret_val = false;
	 }
	 return ret_val;
 }

 vector<int> RunManagerAbstract::get_outstanding_run_ids()
 {
	 vector<int> run_ids;
	 int n_runs = file_stor.get_nruns();
	 for (int id=0; id<n_runs; ++id)
	 {
		 if(run_requried(id))
		 {
			 run_ids.push_back(id);
		 }
	 }
	 return run_ids;
 }
