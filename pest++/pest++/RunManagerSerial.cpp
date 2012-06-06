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
#include "RunManagerSerial.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>
#include <cstring>
#include <map>
#include <direct.h>
#include "Transformable.h"
#include "ModelRunPP.h"
#include "Pest.h"
#include "utilities.h"

using namespace std;
using namespace pest_utils;


extern "C"
{

	void WRTTPL(int *,
		char *,
		char *,
		int *,
		char *,
		double *,
		int *);

	void READINS(int *,
		char *,
		char *,
		int *,
		char *,
		double *,
		int *);
}

string RunManagerSerial::tpl_err_msg(int i)
{
	string err_msg;
	switch (i)
	{
	case 0:
		err_msg = "Routine completed successfully";
		break;
	case 1:
		err_msg = "TPL file does not exist";
		break;
	case 2:
		err_msg = "Not all parameters listed in TPL file";
		break;
	case 3:
		err_msg = "Illegal header specified in TPL file";
		break;
	case 4:
		err_msg = "Error getting parameter name from template";
		break;
	case 5:
		err_msg = "Error getting parameter name from template";
		break;
	case 10:
		err_msg = "Error writing to model input file";
	}
	return err_msg;
}


string RunManagerSerial::ins_err_msg(int i)
{
	string err_msg;
	switch (i)
	{
	case 0:
		err_msg = "Routine completed successfully";
		break;
	default:
		err_msg = "";
	}
	return err_msg;
}

RunManagerSerial::RunManagerSerial(const ModelExecInfo &_model_exec_info, const string &_run_dir)
	: RunManagerAbstract(_model_exec_info), run_dir(_run_dir), max_runs(0)
{
}


void RunManagerSerial::allocate_memory(const Parameters &model_pars, const Observations &obs, int _max_runs)
{
	// For serial runs there is no need to allocate memory.  Just initialize parameter and observation names
	nruns = 0;
	max_runs = _max_runs;
	int npar = model_pars.size();
	int nobs = obs.size();
	pval.clear();
	pval.reserve(npar);
	oval.clear();
	// populate parameter names
	par_name_vec.clear();
	vector<string> names = model_pars.get_keys();
	for (int i=0; i<npar; ++i) {
		par_name_vec.push_back(names[i]);
	}
	// populate observation names
	obs_name_vec.clear();
	names = obs.get_keys();
	for (int i=0; i<nobs; ++i) {
		obs_name_vec.push_back(names[i]);
	}
}


void RunManagerSerial::add_run(const Parameters &model_pars)
{
	int i_run = nruns;
	nruns++;
	assert (model_pars.size() == par_name_vec.size());

	string *name;
	for(int i=0, npar=par_name_vec.size(); i<npar; ++i)
	{
		name = &(par_name_vec[i]);
		pval.push_back (model_pars.get_rec(*name));
	}
}


void RunManagerSerial::run()
{
	int ifail;
	int i_run;
	int success_runs;
	int npar = par_name_vec.size();
	int nobs = obs_name_vec.size();
	int ntpl = tplfile_vec.size();
    int nins = insfile_vec.size();
	stringstream message;

	success_runs = 0;
	failed_runs.clear();
	oval.resize(nruns*nobs,0.0);
	for (i_run=0; i_run<nruns; ++i_run)
	{
		try {
			std::cout << string(message.str().size(), '\b');
			message.str("");
			message << "(" << success_runs << "/" << nruns << " runs complete)";
			std::cout << message.str();
			_chdir(run_dir.c_str());
			WRTTPL(&ntpl, StringvecFortranCharArray(tplfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
				StringvecFortranCharArray(inpfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
				&npar, StringvecFortranCharArray(par_name_vec, 50, pest_utils::TO_LOWER).get_prt(),
				&pval[i_run*npar], &ifail);
			if(ifail != 0)
			{
				throw PestError("Error processing template file");
			}


			for (int i=0, n_exec=comline_vec.size(); i<n_exec; ++i)
			{
				system(comline_vec[i].c_str());
			}

			READINS(&nins, StringvecFortranCharArray(insfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
				StringvecFortranCharArray(outfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
				&nobs, StringvecFortranCharArray(obs_name_vec, 50, pest_utils::TO_LOWER).get_prt(),
				&oval[i_run*nobs], &ifail);
			if(ifail != 0)
			{
				throw PestError("Error processing template file");
			}
			success_runs +=1;
		}
		catch(...)
		{
			failed_runs.insert(i_run);
			cerr << "   Aborting model run" << endl;
		}
	}
	total_runs += success_runs;
	std::cout << string(message.str().size(), '\b');
	message.str("");
	message << "(" << success_runs << "/" << nruns << " runs complete)";
	std::cout << message.str();
	if (success_runs < i_run)
	{
		cout << endl << endl;
		cout << "WARNING: " << i_run-success_runs << " out of " <<i_run << " runs failed" << endl << endl;
	}
}



void RunManagerSerial::get_run(ModelRun &model_run, int run_num, PAR_UPDATE update_type) const
{
	const string *name;
	if (failed_runs.count(run_num) > 0)
	{
		throw PestError("model run failed");
	}
	//Must set parameters before observations
	if(update_type == FORCE_PAR_UPDATE || model_run.get_par_tran().is_one_to_one())
	{
		Parameters new_model_par = model_run.get_model_pars();
		assert(new_model_par.size() == par_name_vec.size());
		for(int i=0, npar=par_name_vec.size(); i<npar; ++i)
		{
			name = &(par_name_vec[i]);
			new_model_par.update_rec(*name, pval[run_num*npar+i]);
		}
		model_run.set_numeric_parameters(model_run.get_par_tran().model2numeric_cp(new_model_par));
	}

	// Process Observations
	Observations new_obs = model_run.get_obs_template();
	assert(new_obs.size() == obs_name_vec.size());
	for(int i=0, nobs=obs_name_vec.size(); i<nobs; ++i)
	{
		name = &(obs_name_vec[i]);
		new_obs.update_rec(*name, oval[run_num*nobs+i]);
	}
	model_run.set_observations(new_obs);

}

 Parameters RunManagerSerial::get_model_parameters(int run_num) const
 {
	Parameters ret_pars;

	for(int i=0, npar = par_name_vec.size(); i<npar; ++i)
	{
		ret_pars[par_name_vec[i]] = pval[run_num*npar+i];
	}
	return ret_pars;
 }

RunManagerSerial::~RunManagerSerial(void)
{
}
