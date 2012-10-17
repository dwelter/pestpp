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


#include <sstream> 
#include <cstdio>
#include <cassert>
#include <iostream>
#include <fstream>
#include "RunStorage.h"
#include "Serialization.h"

using namespace std;

RunStorage::RunStorage(const string &_filename) :filename(_filename), n_runs(0), run_byte_size(0)
{
}

void RunStorage::reset(const Parameters &pars, const Observations &obs)
{
	n_runs = 0;
	default_pars = pars;
	default_obs = obs;
	for(auto &i : default_pars)
	{
		i.second = Transformable::no_data;
	}
	for(auto &i : default_obs)
	{
		i.second = Transformable::no_data;
	}
	free_memory();
	// a file needs to exist before it can beo pened it with read and write 
	// permission.   So open it with write permission to crteate it, close 
	// and then reopen it with read and write permisssion.
	buf_stream.open(filename.c_str(), ios_base::out | ios_base::binary);
	buf_stream.close();
	buf_stream.open(filename.c_str(), ios_base::out | ios_base::in | ios_base::binary);
	assert(buf_stream.good() == true);
	if (!buf_stream.good())
	{
		throw PestFileError(filename);
	}
	// calculate the number of bytes required to store a model run
	vector<char> serial_data  = Serialization::serialize(default_pars);
	pars_byte_size = serial_data.size();
	serial_data = Serialization::serialize(default_obs);
	run_byte_size = pars_byte_size + serial_data.size();
}

int RunStorage::add_run(const Parameters &pars)
{
	int run_id = 0;
	vector<char> serial_data = Serialization::serialize(pars, default_obs);
	n_runs++;
	run_id = n_runs - 1;
	// do error checking
	check_rec_size(serial_data);
	buf_stream.seekp(run_byte_size*run_id, ios_base::beg);
	buf_stream.write(serial_data.data(), serial_data.size());
	return run_id;
}

void RunStorage::update_run(int run_id, const Parameters &pars, const Observations &obs)
{
	vector<char> serial_data = Serialization::serialize(pars, obs);
	update_run(run_id, serial_data);
}

void RunStorage::update_run(int run_id, const vector<char> serial_data)
{
	check_rec_size(serial_data);
	check_rec_id(run_id);
	buf_stream.seekp(run_byte_size*run_id, ios_base::beg);
	buf_stream.write(serial_data.data(), serial_data.size());
}

void RunStorage::get_run(int run_id, Parameters *pars, Observations *obs)
{
	vector<char> serial_data;
	check_rec_id(run_id);
	serial_data.resize(run_byte_size);
	buf_stream.seekg(run_byte_size*run_id, ios_base::beg);
	buf_stream.read(serial_data.data(), serial_data.size());
	Serialization::unserialize(serial_data, *pars, *obs);
}

vector<char> RunStorage::get_serial_pars(int run_id)
{
	vector<char> serial_data;
	check_rec_id(run_id);
	serial_data.resize(pars_byte_size);
	buf_stream.seekg(run_byte_size*run_id, ios_base::beg);
	buf_stream.read(serial_data.data(), serial_data.size());
	return serial_data;
}

void RunStorage::free_memory()
{
	if (buf_stream.is_open()) {
		buf_stream.close();
		remove(filename.c_str());
	}
}

void RunStorage::check_rec_size(const vector<char> &serial_data) const
{
	if (serial_data.size() != run_byte_size)
	{ 
		throw PestError("Error in RunStorage routine.  Size of serial data is different from what is expected");
	}
}

void RunStorage::check_rec_id(int run_id) const
{
	if ( run_id + 1 > n_runs)
	{
		ostringstream msg;
		msg << "Error in RunStorage routine: run id = " << run_id << "is not valid.  Valid values are 0 to " << n_runs - 1;
		throw PestError(msg.str());
	}
}

RunStorage::~RunStorage()
{
	free_memory();
}