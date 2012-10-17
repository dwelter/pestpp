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
#ifndef RUNMANAGERABSTRACT_H
#define RUNMANAGERABSTRACT_H

#include<string>
#include <vector>

class Parameters;
class Observations;
class ModelRun;
class ModelExecInfo;

class RunManagerAbstract
{
public:
	enum PAR_UPDATE{DEFAULT_PAR_UPDATE, FORCE_PAR_UPDATE};
	RunManagerAbstract::RunManagerAbstract(const ModelExecInfo &_mode_exec_info);
	virtual void allocate_memory(const Parameters &pars, const Observations &obs, int _nruns) = 0;
	virtual void free_memory() = 0;
	virtual int add_run(const Parameters &model_pars) = 0;
	virtual void run() = 0;
	virtual void get_run(ModelRun &model_run, int run_num, PAR_UPDATE update_type=DEFAULT_PAR_UPDATE)  = 0;
	virtual Parameters get_model_parameters(int run_num) const = 0;
	virtual Observations get_obs_template(double value = -9999.0) const;
	virtual int get_total_runs(void) const {return total_runs;}
	virtual int get_nruns() {return nruns;}
	virtual ~RunManagerAbstract(void) {}
protected:
	int total_runs;
	int nruns;
	std::vector<std::string> par_name_vec;
	std::vector<std::string> obs_name_vec;
	std::vector<std::string> comline_vec;
	std::vector<std::string> tplfile_vec;
	std::vector<std::string> inpfile_vec;
	std::vector<std::string> insfile_vec;
	std::vector<std::string> outfile_vec;
};

#endif /*  RUNMANAGERABSTRACT_H */