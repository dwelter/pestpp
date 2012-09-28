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

#include <ostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <utility>
#include <cassert>
#include <memory>
#include "Transformable.h"
#include "Transformation.h"
#include "pest_error.h"
#include "pest_data_structs.h"
#include "network_wrapper.h"
#include "utilities.h"

using std::string;
using std::map;
using std::ostream;
using std::endl;
using namespace pest_utils;

const double Transformable::NO_DATA = -9.99E99;

Transformable::Transformable(const Transformable &copyin) : items(copyin.items)
{
}

Transformable::Transformable(const Transformable &copyin, const vector<string> &copy_names)
{
	for (vector<string>::const_iterator b=copy_names.begin(), e=copy_names.end(); b!=e; ++b) {
		items[*b] = copyin.get_rec(*b);
	}
}


const Transformable& Transformable::operator=(const Transformable &rhs)
{
	items = rhs.items;
	return *this;
}

Transformable& Transformable::operator+=(const Transformable &rhs)
{
	for(Transformable::iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		b->second += rhs.items.find(b->first)->second;
	}
    return *this;
 }


Transformable& Transformable::operator-=(const Transformable &rhs)
{
	for(Transformable::iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		b->second -= rhs.items.find(b->first)->second;
	}
    return *this;
 }

Transformable Transformable::operator-(const Transformable &rhs) const
{
	Transformable ret_val(*this);
	ret_val -= rhs;
	return ret_val;

	
}

double &Transformable::operator[](const string &name)
{
	return items[name];
}

pair<Transformable::iterator,bool> Transformable::insert(const string &name, double value)
{
	pair<string, double> rec(name, value);
	return items.insert(rec);
}

pair<Transformable::iterator,bool>  Transformable::insert(const pair<string, double> &x)
{
	return items.insert(x);
}

size_t Transformable::erase(const string &name)
{
	return items.erase(name);
}

Transformable::iterator Transformable::find(const string &name)
{
	return items.find(name);
}

Transformable::const_iterator Transformable::find(const string &name) const
{
	return items.find(name);
}

const double* Transformable::get_rec_ptr(const string &name) const
{
	const double *ret_val = 0;
	Transformable::const_iterator iter = find(name);
	if (iter != end()) {
		ret_val = &((*iter).second);
	}
	return ret_val;
}

const double Transformable::get_rec(const string &name) const
{
	double ret_val = -9999999;
	Transformable::const_iterator iter = find(name);
	if (iter != end()) {
		ret_val = (*iter).second;
	}
	else {
		throw(Transformable_value_error(name));
	}
	return ret_val;
}

void Transformable::update_rec(const string &name, double value)
{
	Transformable::iterator iter = find(name);
	assert(iter != end());
	if (iter != end()) {
		(*iter).second = value;
	}
	else {
		throw(Transformable_value_error(name));
	}
}

vector<double> Transformable::get_vector(const vector<string> &keys) const
{
	vector<double> v;
	v.resize(items.size(), 0.0);

	double value;
	int i = 0;

	for(vector<string>::const_iterator b=keys.begin(), e=keys.end();
		b!=e; ++b, ++i)
	{
		value = get_rec((*b));
		v[i] = value;
	}
	return v;
}

vector<string> Transformable::get_keys() const
{
	vector<string> ret_val;
	for(Transformable::const_iterator b=begin(), e=end();
		b!=e; ++b)
	{
		ret_val.push_back((*b).first);
	}
	return ret_val;
}


void Transformable::add_upgrade(const vector<string> &keys, const LaVectorDouble &del_values)
{
	int i = 0;
	Transformable::iterator found;
	Transformable::iterator not_found=items.end();

	for(vector<string>::const_iterator b=keys.begin(), e=keys.end(); b!=e; ++b, ++i)
	{
		found = items.find(*b);
		assert(found != not_found);
		(*found).second += del_values(i);
	}
}

ostream& operator<<(ostream& out, const Transformable &rhs)
{
	for(Transformable::const_iterator b=rhs.begin(), e=rhs.end();
		b!=e; ++b)
	{
		out << b->first << " = " << b->second << endl;
	}
	return out;
}

char* Transformable::serialize() const
{
	char *buf;
	unsigned long buf_sz = 0;
	unsigned long names_buf_sz = 0;;
	unsigned long data_buf_sz = 0;
	// calculate buffer size
 	for(auto b : *this) 
	{
		names_buf_sz += b.first.size() + 1;
	}
	data_buf_sz = sizeof(double)*items.size();
	buf_sz = sizeof(unsigned long) + names_buf_sz + sizeof(unsigned long) + data_buf_sz;
	// allocate space
	buf = new char[buf_sz];
	// build string with space deliminated names and array of numbers
	stringstream names;
	vector<double> values;
	for (auto b : *this)
	{
		names << " " << (b.first);
		values.push_back(b.second);
	}
	unsigned long n_rec = values.size();
	//write information to buffer
	size_t i_start = 0;
	w_memcpy_s(buf, buf_sz-i_start, &names_buf_sz, sizeof(names_buf_sz));
	i_start += sizeof(names_buf_sz);
	w_memcpy_s(buf+i_start, buf_sz-i_start, names.str().c_str(), names_buf_sz);
	i_start +=names_buf_sz;
	w_memcpy_s(buf+i_start, buf_sz-i_start, &n_rec, sizeof(n_rec));
	i_start +=sizeof(n_rec);
	w_memcpy_s(buf+i_start, buf_sz-i_start, &values[0], data_buf_sz);
	return buf;
}

void Transformable::unserialize(const char *buf)
{
	// delete all existing items
	clear();
	unsigned long names_arg_sz;
	unsigned long n_rec;
	size_t i_start = 0;
	// get size of names record
	i_start = 0;
	w_memcpy_s(&names_arg_sz, sizeof(names_arg_sz), buf+i_start, sizeof(names_arg_sz));
	i_start += sizeof(names_arg_sz);
	unique_ptr<char[]> names_buf(new char[names_arg_sz]);
	w_memcpy_s(names_buf.get(), sizeof(char)*names_arg_sz, buf+i_start, sizeof(char)*names_arg_sz);
	i_start += sizeof(char)*names_arg_sz;
	w_memcpy_s(&n_rec, sizeof(n_rec), buf+i_start, sizeof(n_rec));
	i_start += sizeof(n_rec);
	// build transformable data set
	double *value_ptr = (double *) (buf+i_start);
	vector<string> names_vec;
	tokenize(strip_cp(string(names_buf.get(), names_arg_sz)), names_vec);
	assert(names_vec.size() == n_rec);
	for (unsigned long i=0; i< n_rec; ++i)
	{
		insert(names_vec[i], *(value_ptr+i));
	}
}