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
#ifndef RUNMANAGERGENIE_H
#define RUNMANAGERGENIE_H

#include<string>
#include <vector>
#include <set>
#include "RunManagerAbstract.h"

class Parameters;
class Observations;
class ModelRun;
class ModelExecInfo;

class RunManagerGenie : public RunManagerAbstract
{
public:
	static const int LEN_PARAMETER_NAME;
	static const int LEN_OBSERVATION_NAME;
	RunManagerGenie(const ModelExecInfo &_mode_exec_info, const std::string &_host, const std::string &_id="PPEST");
	virtual void allocate_memory(const Parameters &pars, const Observations &obs, int _nruns);
	virtual void free_memory();
	virtual void add_run(const Parameters &model_pars);
	virtual void run();
	virtual void get_run(ModelRun &model_run, int run_num, PAR_UPDATE update_type=DEFAULT_PAR_UPDATE) const;
	virtual Parameters get_model_parameters(int run_num) const;
	virtual ~RunManagerGenie(void);
protected:
	int max_runs;
	double *pval;
	double *oval;
	std::string id;
	std::string host; 
	std::set<int> failed_runs;  // not implemented yet
};

#endif /* RUNMANAGERGENIE_H */