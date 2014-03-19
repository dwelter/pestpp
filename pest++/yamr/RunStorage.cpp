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
#include <algorithm> 
#include "RunStorage.h"
#include "Serialization.h"
#include "Transformable.h"

using namespace std;

const double RunStorage::no_data = -9999.0;

RunStorage::RunStorage(const string &_filename) :filename(_filename), run_byte_size(0)
{
}

void RunStorage::reset(const vector<string> &_par_names, const vector<string> &_obs_names, const string &_filename)
{
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
	// calculate the number of bytes required to store parameter names
	vector<char> serial_pnames(Serialization::serialize(par_names));
	std::int64_t p_name_size_64 = serial_pnames.size() * sizeof(char);
	// calculate the number of bytes required to store observation names
	vector<char> serial_onames(Serialization::serialize(obs_names));
	std::int64_t o_name_size_64 = serial_onames.size() * sizeof(char);
	// calculate the number of bytes required to store a model run
	run_par_byte_size = par_names.size() * sizeof(double);
	run_data_byte_size = run_par_byte_size + obs_names.size() * sizeof(double);
	//compute the amount of memeory required to store a single model run
	// run_byte_size = size of run_status + size of info_txt + size of info_value + size of parameter oand observation data
	run_byte_size =  sizeof(std::int8_t) + 41*sizeof(char) * sizeof(double) + run_data_byte_size;
	std::int64_t  run_size_64 = run_byte_size;
	beg_run0 = 4 * sizeof(std::int64_t) + serial_pnames.size() + serial_onames.size();
	std::int64_t n_runs_64=0;
	// write header to file
	buf_stream.seekp(0, ios_base::beg);
	buf_stream.write((char*) &n_runs_64, sizeof(n_runs_64));
	buf_stream.write((char*) &run_size_64, sizeof(run_size_64));
	buf_stream.write((char*) &p_name_size_64, sizeof(p_name_size_64));
	buf_stream.write((char*) &o_name_size_64, sizeof(o_name_size_64));
	buf_stream.write(serial_pnames.data(), serial_pnames.size());
	buf_stream.write(serial_onames.data(), serial_onames.size());
}


void RunStorage::init_restart(const std::string &_filename)
{
	filename = _filename;
	par_names.clear();
	obs_names.clear();

	// a file needs to exist before it can be opened it with read and write 
	// permission.   So open it with write permission to crteate it, close 
	// and then reopen it with read and write permisssion.
	if (buf_stream.is_open())
	{
		buf_stream.close();
	}
	buf_stream.open(filename.c_str(), ios_base::out | ios_base::binary);
	buf_stream.close();
	buf_stream.open(filename.c_str(), ios_base::out | ios_base::in | ios_base::binary | ios_base::app);
	assert(buf_stream.good() == true);
	if (!buf_stream.good())
	{
		throw PestFileError(filename);
	}
	// read header
	buf_stream.seekg(0, ios_base::beg);

	std::int64_t n_runs_64;
	buf_stream.read((char*) &n_runs_64, sizeof(n_runs_64));
	
	std::int64_t  run_size_64;
	buf_stream.read((char*) &run_size_64, sizeof(run_size_64));
	run_byte_size = run_size_64;

	std::int64_t p_name_size_64;
	buf_stream.read((char*) &p_name_size_64, sizeof(p_name_size_64));

	std::int64_t o_name_size_64;
	buf_stream.read((char*) &o_name_size_64, sizeof(o_name_size_64));

	vector<char> serial_pnames;
	serial_pnames.resize(p_name_size_64);
	buf_stream.read(serial_pnames.data(), serial_pnames.size());
	Serialization::unserialize(serial_pnames, par_names);

	vector<char> serial_onames;
	serial_onames.resize(o_name_size_64);
	buf_stream.read(serial_onames.data(), serial_onames.size());
	Serialization::unserialize(serial_onames, obs_names);

	beg_run0 = 4 * sizeof(std::int64_t) + serial_pnames.size() + serial_onames.size();
	run_par_byte_size = par_names.size() * sizeof(double);
	run_data_byte_size = run_par_byte_size + obs_names.size() * sizeof(double);
	run_byte_size =  sizeof(std::int8_t) + run_data_byte_size;
}

int RunStorage::get_nruns()
{
	buf_stream.seekg(0, ios_base::beg);
	std::int64_t n_runs_64;
	buf_stream.read((char*) &n_runs_64, sizeof(n_runs_64));
	int n_runs = n_runs_64;
	return n_runs;
}

int RunStorage::increment_nruns()
{
	buf_stream.seekg(0, ios_base::beg);
	std::int64_t n_runs_64;
	buf_stream.read((char*) &n_runs_64, sizeof(n_runs_64));
	++n_runs_64;
	buf_stream.seekp(0, ios_base::beg);
	buf_stream.write((char*) &n_runs_64, sizeof(n_runs_64));
	int n_runs = n_runs_64;
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

 int RunStorage::add_run(const vector<double> &model_pars, const string &info_txt, double info_value)
 {
	std::int8_t r_status = 0;
	int run_id = increment_nruns() - 1;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');
	copy_n(info_txt.begin(), min(info_txt.size(), size_t(info_txt_length)-1) , info_txt_buf.begin());
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.write(reinterpret_cast<char*>(info_txt_buf.data()), sizeof(char)*info_txt_buf.size());
	buf_stream.write(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.write(reinterpret_cast<const char*>(&model_pars[0]), model_pars.size()*sizeof(double));
	return run_id;
 }

 int RunStorage::add_run(const Eigen::VectorXd &model_pars, const string &info_txt, double info_value)
 {
	std::int8_t r_status = 0;
	int run_id = increment_nruns() - 1;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');
	copy_n(info_txt.begin(), min(info_txt.size(), size_t(info_txt_length)-1) , info_txt_buf.begin());
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.write(reinterpret_cast<char*>(info_txt_buf.data()), sizeof(char)*info_txt_buf.size());
	buf_stream.write(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.write(reinterpret_cast<const char*>(&model_pars(0)), model_pars.size()*sizeof(model_pars(0)));
	return run_id;
 }


int RunStorage::add_run(const Parameters &pars, const string &info_txt, double info_value)
{
	vector<double> data(pars.get_data_vec(par_names));
	int run_id = add_run(data, info_txt, info_value);
	return run_id;
}

void RunStorage::update_run(int run_id, const Parameters &pars, const Observations &obs)
{
	//set run status flage to complete
	std::int8_t r_status = 1;
	check_rec_id(run_id);
	vector<double> par_data(pars.get_data_vec(par_names));
	vector<double> obs_data(obs.get_data_vec(obs_names));
	//write data
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	//skip over info_txt and info_value fields
	buf_stream.seekp(sizeof(char)*info_txt_length+sizeof(double), ios_base::cur);
	buf_stream.write(reinterpret_cast<char*>(par_data.data()), par_data.size() * sizeof(double));
	buf_stream.write(reinterpret_cast<char*>(obs_data.data()), obs_data.size() * sizeof(double));
}

void RunStorage::update_run(int run_id, const vector<char> serial_data)
{
	//set run status flage to complete
	std::int8_t r_status = 1;
	check_rec_size(serial_data);
	check_rec_id(run_id);
	//write data
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	//skip over info_txt and info_value fields
	buf_stream.seekp(sizeof(char)*info_txt_length+sizeof(double), ios_base::cur);
	buf_stream.write(serial_data.data(), serial_data.size());
}


void RunStorage::update_run_failed(int run_id)
{
	std::int8_t r_status = get_run_status_native(run_id);
	if (r_status < 1)
	{
		--r_status;
		check_rec_id(run_id);
		//update run status flag
		buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
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

void RunStorage::get_info(int run_id, int &run_status, string &info_txt, double &info_value)
{
	std::int8_t  r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');

	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));

	run_status = r_status;
	info_txt = info_txt_buf.data();
}

int RunStorage::get_run(int run_id, Parameters *pars, Observations *obs, string &info_txt, double &info_value)
{

	std::int8_t  r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');

	size_t n_par = par_names.size();
	size_t n_obs = obs_names.size();
	vector<double> par_data;
	vector<double> obs_data;
	par_data.resize(n_par);
	obs_data.resize(n_obs);

	check_rec_id(run_id);
	
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(&par_data[0]), n_par * sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(&obs_data[0]), n_obs * sizeof(double));
	pars->update(par_names, par_data);
	obs->update(obs_names, obs_data);
	int status = r_status;
	info_txt = info_txt_buf.data();
	return status;
}

int RunStorage::get_run(int run_id, Parameters *pars, Observations *obs)
{
	string info_txt;
	double info_value;
	return get_run(run_id, pars, obs, info_txt, info_value);
}

int RunStorage::get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs, string &info_txt, double &info_value)
{
	std::int8_t r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');

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
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(pars), p_size * sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(obs), o_size * sizeof(double));
	int status = r_status;
	info_txt = info_txt_buf.data();
	return status;
}

int RunStorage::get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs)
{
	string info_txt;
	double info_value;
	return get_run(run_id, pars, npars, obs, nobs, info_txt, info_value);
}

vector<char> RunStorage::get_serial_pars(int run_id)
{
	check_rec_id(run_id);
	std::int8_t r_status;

	vector<char> serial_data;
	serial_data.resize(run_par_byte_size);
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.seekg(sizeof(r_status)+sizeof(char)*info_txt_length+sizeof(double), ios_base::cur);
	buf_stream.read(serial_data.data(), serial_data.size());
	return serial_data;
}

Parameters RunStorage::get_parameters(int run_id)
{
	std::int8_t r_status;

	check_rec_id(run_id);

	size_t n_par = par_names.size();
	vector<double> par_data;
	par_data.resize(n_par);
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.seekg(sizeof(r_status)+sizeof(char)*info_txt_length+sizeof(double), ios_base::cur);
	buf_stream.read(reinterpret_cast<char*>(par_data.data()), n_par*sizeof(double));
	Parameters pars;
	pars.update(par_names, par_data);
	return pars;
}


vector<double> RunStorage::get_observations_vec(int run_id)
{
	std::int8_t r_status;

	check_rec_id(run_id);

	size_t n_par = par_names.size();
	size_t n_obs = obs_names.size();
	vector<double> obs_data;
	obs_data.resize(n_obs);
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.seekg(sizeof(r_status)+sizeof(char)*info_txt_length + sizeof(double)+n_par*sizeof(double), ios_base::cur);
	buf_stream.read(reinterpret_cast<char*>(obs_data.data()), n_obs*sizeof(double));
	return obs_data;
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

void RunStorage::check_rec_id(int run_id)
{
	int n_runs = get_nruns();
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
