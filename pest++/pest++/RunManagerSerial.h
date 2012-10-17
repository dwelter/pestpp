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
#ifndef RUNMANAGERSERIAL_H
#define RUNMANAGERSERIAL_H

#include "runmanagergenie.h"
#include <string>

class RunManagerSerial : public RunManagerAbstract
{
public:
	RunManagerSerial(const ModelExecInfo &_mode_exec_info, const std::string &run_dir);
	virtual void allocate_memory(const Parameters &pars, const Observations &obs, int _nruns);
	virtual void free_memory() {pval.clear();}
	virtual int add_run(const Parameters &model_pars);
	virtual void run();
	virtual void get_run(ModelRun &model_run, int run_num, PAR_UPDATE update_type=DEFAULT_PAR_UPDATE);
	virtual Parameters get_model_parameters(int run_num) const;
	virtual int get_total_runs(void) const {return total_runs;}
	virtual int get_nruns() {return nruns;}
	~RunManagerSerial(void);
private:
	std::string run_dir;
	std::vector<std::string> par_name_vec;
	std::vector<double> pval;
	std::vector<double> oval;
	std::set<int> failed_runs;
	int max_runs;
	static std::string tpl_err_msg(int i);
	static std::string ins_err_msg(int i);
};

#endif /* RUNMANAGERSERIAL_H */