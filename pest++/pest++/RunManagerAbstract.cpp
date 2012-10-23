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
#include "ModelRunPP.h"
#include "Pest.h"
#include "utilities.h"


RunManagerAbstract::RunManagerAbstract(const ModelExecInfo &_model_exec_info, const string &stor_filename)
	 : total_runs(0), file_stor(stor_filename)
{
	comline_vec = _model_exec_info.comline_vec;
	tplfile_vec = _model_exec_info.tplfile_vec;
	inpfile_vec =_model_exec_info.inpfile_vec;
	insfile_vec = _model_exec_info.insfile_vec;
	outfile_vec =_model_exec_info.outfile_vec;
}

void RunManagerAbstract::allocate_memory(const Parameters &model_pars, const Observations &obs)
{
	file_stor.reset(model_pars, obs);
    par_name_vec = model_pars.get_keys();
    obs_name_vec = obs.get_keys();
}


int RunManagerAbstract::add_run(const Parameters &model_pars)
{
	int run_id = file_stor.add_run(model_pars);
	return run_id;
}

void RunManagerAbstract::get_run(ModelRun &model_run, int run_id, PAR_UPDATE update_type)
{
	Parameters pars;
	Observations obs;
	file_stor.get_run(run_id, &pars, &obs);

	//Must set parameters before observations
	if(update_type == FORCE_PAR_UPDATE || model_run.get_par_tran().is_one_to_one())
	{
		// transform to numeric parameters
		model_run.get_par_tran().model2numeric_ip(pars);
		model_run.set_numeric_parameters(pars);
	}

	// Process Observations
	model_run.set_observations(obs);
}


void  RunManagerAbstract::free_memory()
{
	file_stor.free_memory();
}


Parameters RunManagerAbstract::get_model_parameters(int run_id)
 {
	 return file_stor.get_parameters(run_id); 
 }

 Observations RunManagerAbstract::get_obs_template(double value) const
 {
	Observations ret_obs;

	for(int i=0, nobs=obs_name_vec.size(); i<nobs; ++i)
	{
		ret_obs[obs_name_vec[i]] = value;
	}
	return ret_obs;
 }