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

#include <winsock2.h>
#include <ws2tcpip.h>
#include <string>
#include <sstream>
#include <memory>
#include <cassert>
#include "Serialization.h"
#include "Transformable.h"
#include "network_wrapper.h"
#include "utilities.h"

using namespace std;
using namespace pest_utils;


vector<char> Serialization::serialize(const Transformable &tr_data)
{
	vector<char> buf;
	unsigned long buf_sz = 0;
	unsigned long names_buf_sz = 0;;
	unsigned long data_buf_sz = 0;
	// calculate buffer size
 	for(auto &b : tr_data) 
	{
		names_buf_sz += b.first.size() + 1;
	}
	data_buf_sz = sizeof(double)*tr_data.size();
	buf_sz = sizeof(unsigned long) + names_buf_sz + sizeof(unsigned long) + data_buf_sz;
	// allocate space
	buf.resize(buf_sz,'\0');
	// build string with space deliminated names and array of numbers
	stringstream names;
	vector<double> values;
	for (auto &b : tr_data)
	{
		names <<  ' ' << (b.first);
		values.push_back(b.second);
	}
	unsigned long n_rec = values.size();
	//write information to buffer
	size_t i_start = 0;
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &names_buf_sz, sizeof(names_buf_sz));
	i_start += sizeof(names_buf_sz);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, names.str().c_str(), names_buf_sz);
	i_start +=names_buf_sz;
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &n_rec, sizeof(n_rec));
	i_start +=sizeof(n_rec);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &values[0], data_buf_sz);
	return buf;
}

vector<char> Serialization::serialize(const vector< const Transformable*> tr_vec)
{
	vector<char> buf;
	for (auto i : tr_vec)
	{
		vector<char> serial_data = serialize(*i);
		buf.insert(buf.end(), serial_data.begin(), serial_data.end());
	}
	return buf;
}

vector<char> Serialization::serialize(const std::vector<Transformable*> &tr_vec)
{
	vector<const Transformable*> const_data_vec;
	for (auto &i : tr_vec)
	{
		const_data_vec.push_back(i);
	}
	return serialize(const_data_vec);
}

vector<char> Serialization::serialize(const Parameters &pars, const Observations &obs)
{
	 vector<const Transformable*> tr_vec;
	 tr_vec.push_back(&pars);
	 tr_vec.push_back(&obs);
	 return serialize(tr_vec);
}

vector<char> Serialization::serialize(const vector<string> &string_vec)
{
	vector<char> serial_data;
	for (auto &i : string_vec)
	{
		serial_data.insert(serial_data.end(), i.begin(), i.end());
		serial_data.push_back('\0');
	}
	return serial_data;
}

unsigned Serialization::unserialize(const std::vector<char> buf, Transformable &tr_data, unsigned ser_data_loc)
{
	// delete all existing items
	tr_data.clear();
	unsigned long names_arg_sz;
	unsigned long n_rec;
	unsigned long bytes_read;
	size_t i_start = 0;
	// get size of names record
	i_start = ser_data_loc;
	w_memcpy_s(&names_arg_sz, sizeof(names_arg_sz), buf.data()+i_start, sizeof(names_arg_sz));
	i_start += sizeof(names_arg_sz);
	unique_ptr<char[]> names_buf(new char[names_arg_sz]);
	w_memcpy_s(names_buf.get(), sizeof(char)*names_arg_sz,  buf.data()+i_start, sizeof(char)*names_arg_sz);
	i_start += sizeof(char)*names_arg_sz;
	w_memcpy_s(&n_rec, sizeof(n_rec),  buf.data()+i_start, sizeof(n_rec));
	i_start += sizeof(n_rec);
	// build transformable data set
	double *value_ptr = (double *) ( buf.data()+i_start);
	vector<string> names_vec;
	tokenize(strip_cp(string(names_buf.get(), names_arg_sz)), names_vec);
	assert(names_vec.size() == n_rec);
	for (unsigned long i=0; i< n_rec; ++i)
	{
		tr_data.insert(names_vec[i], *(value_ptr+i));
	}
	bytes_read = i_start + n_rec * sizeof(double);
	return bytes_read;
}

unsigned Serialization::unserialize(const std::vector<char> ser_data, std::vector<Transformable*> &tr_vec)
{
	unsigned i_tr=0;
	unsigned i_char=0;
	unsigned long bytes_read = 0;
	unsigned long total_bytes_read = 0;


	while(i_tr < tr_vec.size() && i_char < ser_data.size())
	{
		bytes_read = Serialization::unserialize(ser_data, *tr_vec[i_tr], i_char);
		i_char +=  bytes_read;
		++i_tr;
		total_bytes_read += bytes_read;
	}
	return total_bytes_read;
}

unsigned Serialization::unserialize(const std::vector<char> data, Parameters &pars, Observations &obs)
{
	 unsigned total_bytes_read = 0;
	 vector<Transformable*> tr_vec;
	 tr_vec.push_back(&pars);
	 tr_vec.push_back(&obs);
	 total_bytes_read = unserialize(data, tr_vec);
	 return total_bytes_read;
}

static unsigned unserialize(const vector<char> ser_data, vector<string> &string_vec)
{
	unsigned total_bytes_read = 0;
	string tmp_str(ser_data.begin(), ser_data.end());
	string delm;
	delm.push_back('\0');
	strip_ip(tmp_str, "both", delm);
	tokenize(tmp_str, string_vec, delm);
	return ser_data.size();
}