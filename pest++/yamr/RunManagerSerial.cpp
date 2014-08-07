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
#include <algorithm>
#include "system_variables.h"
#include "Transformable.h"
#include "utilities.h"
#include "iopp.h"

using namespace std;
using namespace pest_utils;


extern "C"
{

	void wrttpl_(int *,
		char *,
		char *,
		int *,
		char *,
		double *,
		int *);

	void readins_(int *,
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
	const string &stor_filename, const string &_run_dir, int _max_run_fail,
	const bool _io_fortran)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
	_insfile_vec, _outfile_vec, stor_filename, _max_run_fail),
	run_dir(_run_dir), io_fortran(_io_fortran)
{
	cout << "              starting serial run manager ..." << endl << endl;
}

void RunManagerSerial::run()
{
	int ifail;
	int success_runs = 0;
	int prev_sucess_runs = 0;
	const vector<string> &par_name_vec = file_stor.get_par_name_vec();
	const vector<string> &obs_name_vec = file_stor.get_obs_name_vec();
	int npar = par_name_vec.size();
	int nobs = obs_name_vec.size();
	int ntpl = tplfile_vec.size();
	int nins = insfile_vec.size();
	stringstream message;
	bool isDouble = true;
	bool forceRadix = true;
	TemplateFiles tpl_files(isDouble,forceRadix,tplfile_vec,inpfile_vec,par_name_vec);
	//InstructionFiles ins_files(insfile_vec,outfile_vec,obs_name_vec);
	std::vector<double> obs_vec;
	// This is necessary to support restart as some run many already be complete
	vector<int> run_id_vec;
	int nruns = get_outstanding_run_ids().size();
	while (!(run_id_vec = get_outstanding_run_ids()).empty())
	{
		for (int i_run : run_id_vec)
		{
			Observations obs;
			vector<double> par_values;
			Parameters pars;
			file_stor.get_parameters(i_run, pars);

			for (auto &i : par_name_vec)
			{
				par_values.push_back(pars.get_rec(i));
			}
			try {
				std::cout << string(message.str().size(), '\b');
				message.str("");
				message << "(" << success_runs << "/" << nruns << " runs complete)";
				std::cout << message.str();
				OperSys::chdir(run_dir.c_str());
				if (std::any_of(par_values.begin(), par_values.end(), OperSys::double_is_invalid))
				{
					throw PestError("Error running model: invalid parameter value returned");
				}
				if (io_fortran)
				{
					wrttpl_(&ntpl, StringvecFortranCharArray(tplfile_vec, 50).get_prt(),
						StringvecFortranCharArray(inpfile_vec, 50).get_prt(),
						&npar, StringvecFortranCharArray(par_name_vec, 50, pest_utils::TO_LOWER).get_prt(),
						&par_values[0], &ifail);
					if (ifail != 0)
					{
						throw PestError("Error processing template file");
					}
				}
				else
				{					
					tpl_files.write(par_name_vec, par_values);

					//throw PestError("non-fortran IO not implemented for TPL files");
				}
								
				for (int i = 0, n_exec = comline_vec.size(); i < n_exec; ++i)
				{
					system(comline_vec[i].c_str());
				}
				obs_vec.resize(nobs, RunStorage::no_data);
				if (io_fortran)
				{
					readins_(&nins, StringvecFortranCharArray(insfile_vec, 50).get_prt(),
						StringvecFortranCharArray(outfile_vec, 50).get_prt(),
						&nobs, StringvecFortranCharArray(obs_name_vec, 50, pest_utils::TO_LOWER).get_prt(),
						&obs_vec[0], &ifail);
				}
				else
				{
					throw PestError("non-fortran IO not implemented for INS files");
				}
				if (ifail != 0)
				{
					throw PestError("Error processing instruction file");
				}
				//obs_vec = ins_files.readins();

				// check parameters and observations for inf and nan
				if (std::any_of(par_values.begin(), par_values.end(), OperSys::double_is_invalid))
				{
					throw PestError("Error running model: invalid parameter value returned");
				}
				if (std::any_of(obs_vec.begin(), obs_vec.end(), OperSys::double_is_invalid))
				{
					throw PestError("Error running model: invalid observation value returned");
				}
				success_runs += 1;
				pars.clear();
				pars.insert(par_name_vec, par_values);
				obs.clear();
				obs.insert(obs_name_vec, obs_vec);
				file_stor.update_run(i_run, pars, obs);
			}
			catch (const std::exception& ex)
			{
				file_stor.update_run_failed(i_run);
				cerr << endl;
				cerr << "  " << ex.what() << endl;
				cerr << "  Aborting model run" << endl << endl;
			}
			catch (...)
			{
				file_stor.update_run_failed(i_run);
				cerr << endl;
				cerr << "  Error running model" << endl;
				cerr << "  Aborting model run" << endl << endl;
			}
		}
	}
	total_runs += success_runs;
	std::cout << string(message.str().size(), '\b');
	message.str("");
	message << "(" << success_runs << "/" << nruns << " runs complete";
	if (prev_sucess_runs > 0)
	{
		message << " and " << prev_sucess_runs << " additional run completed previously";
	}
	message << ")";
	std::cout << message.str();
	if (success_runs < nruns)
	{			cout << endl << endl;
		cout << "WARNING: " << nruns - success_runs << " out of " << nruns << " runs failed" << endl << endl;
	}
}


RunManagerSerial::~RunManagerSerial(void)
{
}
