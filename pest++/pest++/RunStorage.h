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

#ifndef RUN_STORAGE_H_
#define RUN_STORAGE_H_

#include <string>
#include <fstream>
#include "Transformable.h"

class RunStorage {
public:
	RunStorage(const string &_filename);
	void reset(const Parameters &pars, const Observations &obs);
	int add_run(const Parameters &pars);
	void update_run(int run_id, const Parameters &pars, const Observations &obs);
	void update_run(int run_id, const vector<char> serial_data);
	int get_nruns();
	void get_run(int run_id, Parameters *pars, Observations *obs);
	Parameters get_parameters(int run_id);
	vector<char> get_serial_pars(int run_id);
	void free_memory();
	~RunStorage();
private:
	std::string filename;
	std::fstream buf_stream;
	unsigned long run_byte_size;
	unsigned long pars_byte_size;
	int n_runs;
	Parameters default_pars;
	Observations default_obs;
	void check_rec_size(const vector<char> &serial_data) const;
	void check_rec_id(int run_id) const;
};

#endif //RUN_STORAGE_H_