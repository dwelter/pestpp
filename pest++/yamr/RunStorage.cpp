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
#include "Transformable.h"

using namespace std;

RunStorage::RunStorage(const string &_filename) :filename(_filename), n_runs(0), run_byte_size(0)
{
}

void RunStorage::reset(const vector<string> &_par_names, const vector<string> &_obs_names, const string &_filename)
{
	n_runs = 0;
	par_names = _par_names;
	obs_names = _obs_names;
	// a file needs to exist before it can be opened it with read and write 
	// permission.   So open it with write permission to crteate it, close 
	// and then reopen it with read and write permisssion.
	if (_filename.size() > 0)
	{
		filename = _filename;
	}
	if (buf_stream.is_open())
	{
		buf_stream.close();
	}
	buf_stream.open(filename.c_str(), ios_base::out | ios_base::binary);
	buf_stream.close();
	buf_stream.open(filename.c_str(), ios_base::out | ios_base::in | ios_base::binary);
	assert(buf_stream.good() == true);
	if (!buf_stream.good())
	{
		throw PestFileError(filename);
	}
	// calculate the number of bytes required to store a model run
	run_par_byte_size = par_names.size() * sizeof(double);
	run_data_byte_size = run_par_byte_size + obs_names.size() * sizeof(double);
	run_byte_size =  sizeof(std::int8_t) + run_data_byte_size;
	vector<char> serial_pnames(Serialization::serialize(par_names));
	vector<char> serial_onames(Serialization::serialize(obs_names));
	beg_obs_name = serial_pnames.size();
	beg_run0 = beg_obs_name + serial_onames.size();
	buf_stream.seekp(0, ios_base::beg);
	buf_stream.write(serial_pnames.data(), serial_pnames.size());
	buf_stream.write(serial_onames.data(), serial_onames.size());
}

int RunStorage::get_nruns()
{
	return n_runs;
}

const std::vector<string>& RunStorage::get_par_name_vec()const
{
	return par_names;
}

const std::vector<string>& RunStorage::get_obs_name_vec()const
{
	return obs_names;
}

streamoff RunStorage::get_stream_pos(int run_id)
{
	streamoff pos = beg_run0 + run_byte_size*run_id;
	return pos;
}

 int RunStorage::add_run(const vector<double> &model_pars)
 {
	n_runs++;
	int run_id = n_runs - 1;
	std::int8_t r_status = 0;
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.write(reinterpret_cast<const char*>(&model_pars[0]), model_pars.size()*sizeof(double));
	return run_id;
 }

 int RunStorage::add_run(const Eigen::VectorXd &model_pars)
 {
	n_runs++;
	int run_id = n_runs - 1;
	std::int8_t r_status = 0;
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.write(reinterpret_cast<const char*>(&model_pars(0)), model_pars.size()*sizeof(model_pars(0)));
	return run_id;
 }


int RunStorage::add_run(const Parameters &pars)
{
	vector<double> data(pars.get_data_vec(par_names));
	int run_id = add_run(data);
	return run_id;
}

void RunStorage::update_run(int run_id, const Parameters &pars, const Observations &obs)
{
	check_rec_id(run_id);
	vector<double> par_data(pars.get_data_vec(par_names));
	vector<double> obs_data(obs.get_data_vec(obs_names));
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	//set run status flage to complete
	std::int8_t r_status = 1;
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	//write data
	buf_stream.write(reinterpret_cast<char*>(par_data.data()), par_data.size() * sizeof(double));
	buf_stream.write(reinterpret_cast<char*>(obs_data.data()), obs_data.size() * sizeof(double));
}

void RunStorage::update_run(int run_id, const vector<char> serial_data)
{
	check_rec_size(serial_data);
	check_rec_id(run_id);
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	//set run status flage to complete
	std::int8_t r_status = 1;
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	//write data
	buf_stream.write(serial_data.data(), serial_data.size());
}


void RunStorage::update_run_failed(int run_id)
{
	std::int8_t r_status = get_run_status_native(run_id);
	if (r_status < 1)
	{
		--r_status;
		check_rec_id(run_id);
		buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
		//update run status flag
		buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	}


}

std::int8_t RunStorage::get_run_status_native(int run_id)
{
	std::int8_t  r_status;
	check_rec_id(run_id);
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	return r_status;
}

int RunStorage::get_run_status(int run_id)
{
	int status = get_run_status_native(run_id);
	return status;
}

int RunStorage::get_run(int run_id, Parameters *pars, Observations *obs)
{

	std::int8_t  r_status;
	size_t n_par = par_names.size();
	size_t n_obs = obs_names.size();
	vector<double> par_data;
	vector<double> obs_data;
	par_data.resize(n_par);
	obs_data.resize(n_obs);

	check_rec_id(run_id);
	
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&par_data[0]), n_par * sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(&obs_data[0]), n_obs * sizeof(double));
	pars->update(par_names, par_data);
	obs->update(obs_names, obs_data);
	int status = r_status;
	return status;
}

int RunStorage::get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs)
{
	std::int8_t r_status;
	check_rec_id(run_id);

	size_t p_size = par_names.size();
	size_t o_size = obs_names.size();

	assert(npars == p_size);
	assert(nobs == o_size);

	if (npars != npars) {
		throw(PestIndexError("RunStorage::get_run: parameter dimension in incorrect"));
	}
	if (nobs != nobs) {
		throw(PestIndexError("RunStorage::get_run: observation dimension in incorrect"));
	}

	p_size = min(p_size, npars);
	o_size = min(o_size, nobs);
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(pars), p_size * sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(obs), o_size * sizeof(double));
	int status = r_status;
	return status;
}


vector<char> RunStorage::get_serial_pars(int run_id)
{
	check_rec_id(run_id);
	std::int8_t r_status;

	vector<char> serial_data;
	serial_data.resize(run_par_byte_size);
	buf_stream.seekg(run_byte_size*run_id, ios_base::beg);
	buf_stream.seekg(get_stream_pos(run_id)+sizeof(r_status), ios_base::beg);
	buf_stream.read(serial_data.data(), serial_data.size());
	return serial_data;
}

Parameters RunStorage::get_parameters(int run_id)
{
	check_rec_id(run_id);

	size_t n_par = par_names.size();
	std::int8_t r_status;
	vector<double> par_data;
	par_data.resize(n_par);
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(par_data.data()), n_par*sizeof(double));
	Parameters pars;
	pars.update(par_names, par_data);

	return pars;
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
	if (serial_data.size() != run_data_byte_size)
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