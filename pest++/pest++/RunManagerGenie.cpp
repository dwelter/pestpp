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
#include "RunManagerGenie.h"
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

using namespace std;
using namespace pest_utils;

const int RunManagerGenie::LEN_PARAMETER_NAME = 12+1;
const int RunManagerGenie::LEN_OBSERVATION_NAME = 20+1;

extern "C" { 
	int GENIE_INTERFACE(int*, int*, char*,
		int*, int*, char*, char*, double*, double*,
		int*, int*, char*, char*, char*,
		char*, char*, char*, int*); 
	int GENIE_KILL_GMAN(char*,char*);
}



RunManagerGenie::RunManagerGenie(const ModelExecInfo &_model_exec_info, 
	const std::string &_host, const std::string &_id)
: RunManagerAbstract(_model_exec_info), pval(0), oval(0), host(_host), id(_id)
{
}

void RunManagerGenie::allocate_memory(const Parameters &model_pars, const Observations &obs, int _max_runs)
{
	max_runs = _max_runs;
	nruns = 0;
	int npar = model_pars.size();
	int nobs = obs.size();

	assert(max_runs > 0);
	try {
		pval = new double[npar*max_runs];
		oval = new double[nobs*max_runs];
	}
	catch(bad_alloc e)
	{
		cerr << "Error: can not allocate memeory required to run GENIE." << endl;
		cerr << "  program location: Jacobian::calculate_genie" << endl;
		throw e;
	}
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

void RunManagerGenie::add_run(const Parameters &model_pars)
{
	int i_run = nruns;
	nruns++;
	assert (nruns <=max_runs);
	assert (model_pars.size() == par_name_vec.size());

	string *name;
	for(int i=0, npar=par_name_vec.size(); i<npar; ++i)
	{
		name = &(par_name_vec[i]);
		pval[i_run*npar+i] = model_pars.get_rec(*name);
	}
}

void RunManagerGenie::run()
{
	int nexec = comline_vec.size();
	int npar = par_name_vec.size();
	int nobs = obs_name_vec.size();
	int ntpl = tplfile_vec.size();
    int nins = insfile_vec.size();
	int ikill = 1;
	stringstream apar;
	stringstream aobs;
	stringstream execnames;
	stringstream tplfle;
	stringstream infle;
	stringstream insfle;
	stringstream outfle;

	std::copy(par_name_vec.begin(), par_name_vec.end(),std::ostream_iterator<std::string>(apar,"\n"));
	apar.str(lower_cp(apar.str()));
	std::copy(obs_name_vec.begin(), obs_name_vec.end(),std::ostream_iterator<std::string>(aobs,"\n"));
	aobs.str(lower_cp(aobs.str()));
	std::copy(comline_vec.begin(), comline_vec.end(),std::ostream_iterator<std::string>(execnames,"\n"));
	//execnames.str(lower_cp(execnames.str()));
	std::copy(tplfile_vec.begin(), tplfile_vec.end(),std::ostream_iterator<std::string>(tplfle,"\n"));
	//tplfle.str(lower_cp(tplfle.str()));
	std::copy(inpfile_vec.begin(), inpfile_vec.end(),std::ostream_iterator<std::string>(infle,"\n"));
	//infle.str(lower_cp(infle.str()));
	std::copy(insfile_vec.begin(), insfile_vec.end(),std::ostream_iterator<std::string>(insfle,"\n"));
	//insfle.str(lower_cp(insfle.str()));
	std::copy(outfile_vec.begin(), outfile_vec.end(),std::ostream_iterator<std::string>(outfle,"\n"));
	//outfle.str(lower_cp(outfle.str()));

	failed_runs.clear();  //not implemented yet
	GENIE_INTERFACE(&nruns, &nexec, String2CharPtr(execnames.str()).get_char_ptr(), &npar, &nobs,
		String2CharPtr(apar.str()).get_char_ptr(), String2CharPtr(aobs.str()).get_char_ptr(),
                    pval, oval, &ntpl, &nins, String2CharPtr(tplfle.str()).get_char_ptr(),
                    String2CharPtr(infle.str()).get_char_ptr(), String2CharPtr(insfle.str()).get_char_ptr(),
					String2CharPtr(outfle.str()).get_char_ptr(),
                    String2CharPtr(host).get_char_ptr(), String2CharPtr(id).get_char_ptr(), &ikill);
	total_runs += nruns;
}


void RunManagerGenie::get_run(ModelRun &model_run, int run_num, PAR_UPDATE update_type) const
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


 Parameters RunManagerGenie::get_model_parameters(int run_num) const
 {
	Parameters ret_pars;

	for(int i=0, npar = par_name_vec.size(); i<npar; ++i)
	{
		ret_pars[par_name_vec[i]] = pval[run_num*npar+i];
	}
	return ret_pars;
 }


void RunManagerGenie::free_memory()
{
	if (pval !=0) {
		delete[] pval;
		pval = 0;
	}
	if (oval !=0) {
		delete[] oval;
		oval = 0;
	}
}


RunManagerGenie::~RunManagerGenie(void)
{
	//close GMAN run manager
	cout << "Closing GMAN run manager..." << endl;
	GENIE_KILL_GMAN(String2CharPtr(id).get_char_ptr(), String2CharPtr(host).get_char_ptr());

	free_memory();
}

