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
#include "system_variables.h"
#include "Transformable.h"
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

RunManagerSerial::RunManagerSerial(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, const string &_run_dir)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
		_insfile_vec, _outfile_vec, stor_filename),
		run_dir(_run_dir)
{

}

void RunManagerSerial::run()
{
	int ifail;
	int i_run;
	int success_runs;
    const vector<string> &par_name_vec = file_stor.get_par_name_vec();
    const vector<string> &obs_name_vec = file_stor.get_obs_name_vec();
	int npar = par_name_vec.size();
	int nobs = obs_name_vec.size();
	int ntpl = tplfile_vec.size();
    int nins = insfile_vec.size();
	stringstream message;

	success_runs = 0;
	failed_runs.clear();
    int nruns = file_stor.get_nruns();
	for (i_run=0; i_run<nruns; ++i_run)
	{
        Parameters pars = file_stor.get_parameters(i_run);
        Observations obs;
		vector<double> par_values;
		for(auto &i : par_name_vec)
		{
			par_values.push_back(pars.get_rec(i));
		}
		try {
			std::cout << string(message.str().size(), '\b');
			message.str("");
			message << "(" << success_runs << "/" << nruns << " runs complete)";
			std::cout << message.str();
			OperSys::chdir(run_dir.c_str());
			WRTTPL(&ntpl, StringvecFortranCharArray(tplfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
				StringvecFortranCharArray(inpfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
				&npar, StringvecFortranCharArray(par_name_vec, 50, pest_utils::TO_LOWER).get_prt(),
				&par_values[0], &ifail);
			if(ifail != 0)
			{
				throw PestError("Error processing template file");
			}


			for (int i=0, n_exec=comline_vec.size(); i<n_exec; ++i)
			{
				system(comline_vec[i].c_str());
			}

		    std::vector<double> obs_vec;
		    obs_vec.resize(nobs, -9999.00);
			READINS(&nins, StringvecFortranCharArray(insfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
				StringvecFortranCharArray(outfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
				&nobs, StringvecFortranCharArray(obs_name_vec, 50, pest_utils::TO_LOWER).get_prt(),
				&obs_vec[0], &ifail);
			if(ifail != 0)
			{
				throw PestError("Error processing template file");
			}
			success_runs +=1;
            pars.clear();
            pars.insert(par_name_vec, par_values);
            obs.clear();
            obs.insert(obs_name_vec, obs_vec);
            file_stor.update_run(i_run, pars, obs);
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


RunManagerSerial::~RunManagerSerial(void)
{
}
