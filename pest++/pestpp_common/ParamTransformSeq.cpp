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

#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "Transformable.h"
#include "Jacobian.h"

using namespace std;

//These routines handle reference count on the Transformation used by ParamTransformSeq.
// Using reference counting makes it possible to not store new copies of all the
// transformations each time a ParamTransformSeq is copied.

map<const Transformation*, int> ParamTransformSeq::tran_ref_count = map<const Transformation*, int>();

int ParamTransformSeq::tran_add_ref_count(const Transformation *new_tran)
{
	int ref_count = 0;
	map<const Transformation *,int>::iterator it;
	it = tran_ref_count.find(new_tran);
	if (it != tran_ref_count.end()) {
		(*it).second += 1;
		ref_count = (*it).second;
	}
	else {
		tran_ref_count.insert(pair<const Transformation *, int>(new_tran, 1));
		ref_count = 1;
	}
	return ref_count;
}

int ParamTransformSeq::tran_sub_ref_count(const Transformation *new_tran)
{
	int ref_count = 0;
	const Transformation *t_ptr;
	map<const Transformation *,int>::iterator it;
	it = tran_ref_count.find(new_tran);
	if (it != tran_ref_count.end()) {
		(*it).second -= 1;
		ref_count = (*it).second;
		if (ref_count<=0) {
			t_ptr = (*it).first;
			tran_ref_count.erase(it);
			delete t_ptr;
		}
	}
	return ref_count;
}


ParamTransformSeq::ParamTransformSeq(const ParamTransformSeq &rhs)
{
	copy(rhs);
}

ParamTransformSeq::ParamTransformSeq(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set)
{
	copy(rhs, deep_copy_tran_set);
}

ParamTransformSeq::~ParamTransformSeq()
{
	clear();
}

void ParamTransformSeq::clear()
{
	clear_tranSeq_ctl2model();
	clear_tranSeq_ctl2derivative();
	clear_tranSeq_derivative2numeric();
	default_deep_copy_tran_set.clear();
	name = "empty";
}

void ParamTransformSeq::clear_tranSeq_ctl2model()
{
	for(vector<Transformation*>::iterator i = tranSeq_ctl2model.begin(),
		e=tranSeq_ctl2model.end(); i != e; ++i)
	{
		tran_sub_ref_count(*i);
		default_deep_copy_tran_set.erase(*i);
	}
	ctl_offset_ptr = 0;
	ctl_scale_prt = 0;
	tranSeq_ctl2model.clear();
}

void  ParamTransformSeq::clear_tranSeq_ctl2derivative()
{
	for(vector<Transformation*>::iterator i = tranSeq_ctl2derivative.begin(),
		e=tranSeq_ctl2derivative.end(); i != e; ++i)
	{
		tran_sub_ref_count(*i);
		default_deep_copy_tran_set.erase(*i);
	}
	ctl_log10_ptr = 0;
	tranSeq_ctl2derivative.clear();
}


void ParamTransformSeq::clear_tranSeq_derivative2numeric()
{
	for(vector<Transformation*>::iterator i = tranSeq_derivative2numeric.begin(),
		e=tranSeq_derivative2numeric.end(); i != e; ++i)
	{
		tran_sub_ref_count(*i);
		default_deep_copy_tran_set.erase(*i);
	}
	svda_ptr = 0;
	tranSeq_derivative2numeric.clear();
}

void ParamTransformSeq::deep_copy(const ParamTransformSeq &rhs)
{
	set<Transformation *> deep_copy_set(rhs.tranSeq_ctl2model.begin(), rhs.tranSeq_ctl2model.end());
	deep_copy_set.insert(rhs.tranSeq_ctl2derivative.begin(), rhs.tranSeq_ctl2derivative.end());
	deep_copy_set.insert(rhs.tranSeq_derivative2numeric.begin(), rhs.tranSeq_derivative2numeric.end());
	copy(rhs, deep_copy_set);
}

ParamTransformSeq& ParamTransformSeq::operator=(const ParamTransformSeq &rhs)
{
	copy(rhs);
	return *this;
}

void ParamTransformSeq::copy(const ParamTransformSeq &rhs)
{
	copy(rhs, rhs.default_deep_copy_tran_set);
}

void ParamTransformSeq::copy(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set)
{
	clear();
	name = "copy of " + rhs.name;
	
	Transformation *t_ptr;
	for(vector<Transformation*>::const_iterator i = rhs.tranSeq_ctl2model.begin(),
		e=rhs.tranSeq_ctl2model.end(); i != e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_ctl2model(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
		if (*i == rhs.ctl_offset_ptr) ctl_offset_ptr = dynamic_cast<TranOffset*>(t_ptr);
		if (*i == rhs.ctl_scale_prt) ctl_scale_prt = dynamic_cast<TranScale*>(t_ptr);
	}
	for(vector<Transformation*>::const_iterator i = rhs.tranSeq_ctl2derivative.begin(),
		e=rhs.tranSeq_ctl2derivative.end(); i != e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_ctl2derivative(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
		if (*i == rhs.ctl_log10_ptr) ctl_log10_ptr = dynamic_cast<TranLog10*>(t_ptr);
	}
	for(vector<Transformation*>::const_iterator i = rhs.tranSeq_derivative2numeric.begin(),
		e=rhs.tranSeq_derivative2numeric.end(); i != e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_derivative2numeric(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
		if (*i == rhs.ctl_log10_ptr) ctl_log10_ptr = dynamic_cast<TranLog10*>(t_ptr);
		if (*i == rhs.svda_ptr) svda_ptr = dynamic_cast<TranSVD*>(t_ptr);
	}
	for (auto &itran_seq : rhs.custom_tran_seq ) 
	{
	  for (auto &itran : itran_seq.second )
	  { 
		if (deep_copy_tran_set.find(itran) ==  deep_copy_tran_set.end()) {
			t_ptr = itran;
		}
		else {
			t_ptr = itran->clone();
		}
		custom_tran_seq[itran_seq.first].push_back(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(itran) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	  }
	}
}

ParamTransformSeq &ParamTransformSeq::operator+=(const ParamTransformSeq &rhs)
{
	append(rhs, default_deep_copy_tran_set);
	return *this;
}

void ParamTransformSeq::append(const ParamTransformSeq &rhs)
{
	append(rhs, rhs.default_deep_copy_tran_set);
}

void ParamTransformSeq::append(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set )
{
	Transformation *t_ptr;

	for (vector<Transformation*>::const_iterator i=rhs.tranSeq_ctl2model.begin(), e=rhs.tranSeq_ctl2model.end();
		i!=e; ++i) 
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_ctl2model(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for (vector<Transformation*>::const_iterator i=rhs.tranSeq_ctl2derivative.begin(), e=rhs.tranSeq_ctl2derivative.end();
		i!=e; ++i) 
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_ctl2derivative(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for (vector<Transformation*>::const_iterator i=rhs.tranSeq_derivative2numeric.begin(), e=rhs.tranSeq_derivative2numeric.end();
		i!=e; ++i) 
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_derivative2numeric(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for (auto &itran_seq : rhs.custom_tran_seq ) 
	{
	  for (auto &itran : itran_seq.second )
	  { 
		if (deep_copy_tran_set.find(itran) ==  deep_copy_tran_set.end()) {
			t_ptr = itran;
		}
		else {
			t_ptr = itran->clone();
		}
		custom_tran_seq[itran_seq.first].push_back(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(itran) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	  }
	}
}

void ParamTransformSeq::push_back_ctl2model(Transformation *tr)
{
	tranSeq_ctl2model.push_back(tr);
	tran_add_ref_count(tr);
}

void ParamTransformSeq::push_back_ctl2derivative(Transformation *tr)
{
	tranSeq_ctl2derivative.push_back(tr);
	tran_add_ref_count(tr);
}


void ParamTransformSeq::push_back_derivative2numeric(Transformation *tr)
{
	tranSeq_derivative2numeric.push_back(tr);
	tran_add_ref_count(tr);
}


void ParamTransformSeq::ctl2model_ip(Parameters &data) const
{
	vector<Transformation*>::const_iterator iter, e;

	for(iter = tranSeq_ctl2model.begin(), e = tranSeq_ctl2model.end();
		iter != e; ++iter)
	{
		(*iter)->forward(data);
	}
}

Parameters ParamTransformSeq::ctl2model_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	ctl2model_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::ctl2derivative_ip(Parameters &data) const
{
	vector<Transformation*>::const_iterator iter, e;
	for(iter = tranSeq_ctl2derivative.begin(), e = tranSeq_ctl2derivative.end();
		iter != e; ++iter)
	{
		(*iter)->forward(data);
		
	}
}

Parameters ParamTransformSeq::ctl2derivative_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	ctl2derivative_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::ctl2numeric_ip(Parameters &data) const
{
	 ctl2derivative_ip(data);

	vector<Transformation*>::const_iterator iter, e;
	for(iter = tranSeq_derivative2numeric.begin(), e = tranSeq_derivative2numeric.end();
		iter != e; ++iter)
	{
		(*iter)->forward(data);
		
	}
}

Parameters ParamTransformSeq::ctl2numeric_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	ctl2numeric_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::model2ctl_ip(Parameters &data) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for(iter = tranSeq_ctl2model.rbegin(), e = tranSeq_ctl2model.rend();
		iter != e; ++iter)
	{
		(*iter)->reverse(data);
	}
}

Parameters ParamTransformSeq::model2ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	model2ctl_ip(ret_val);
	return ret_val;
}


void ParamTransformSeq::numeric2derivative_ip(Parameters &data) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for(iter = tranSeq_derivative2numeric.rbegin(), e = tranSeq_derivative2numeric.rend();
		iter != e; ++iter)
	{
		(*iter)->reverse(data);
	}
}


Parameters ParamTransformSeq::numeric2derivative_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	numeric2derivative_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::numeric2ctl_ip(Parameters &data) const
{
	numeric2derivative_ip(data);
	vector<Transformation*>::const_reverse_iterator iter, e;
	for(iter = tranSeq_ctl2derivative.rbegin(), e = tranSeq_ctl2derivative.rend();
		iter != e; ++iter)
	{
		(*iter)->reverse(data);
	}
}

Parameters ParamTransformSeq::numeric2ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	numeric2ctl_ip(ret_val);
	return ret_val;
}


void ParamTransformSeq::numeric2model_ip(Parameters &data) const
{
	numeric2ctl_ip(data);
	ctl2model_ip(data);
}

Parameters ParamTransformSeq::numeric2model_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	numeric2model_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::model2numeric_ip(Parameters &data) const
{
	model2ctl_ip(data);
	ctl2numeric_ip(data);
}

Parameters ParamTransformSeq::model2numeric_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	model2numeric_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::derivative2numeric_ip(Parameters &data) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for(iter = tranSeq_derivative2numeric.rbegin(), e = tranSeq_derivative2numeric.rend();
		iter != e; ++iter)
	{
		(*iter)->forward(data);
	}
}

Parameters ParamTransformSeq::derivative2numeric_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	derivative2numeric_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::derivative2ctl_ip(Parameters &data) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for(iter = tranSeq_ctl2derivative.rbegin(), e = tranSeq_ctl2derivative.rend();
		iter != e; ++iter)
	{
		(*iter)->reverse(data);
	}
}

Parameters ParamTransformSeq::derivative2ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	derivative2ctl_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::derivative2model_ip(Parameters &data) const
{
	derivative2ctl_ip(data);
	ctl2model_ip(data);
}

Parameters ParamTransformSeq::derivative2model_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	derivative2model_ip(ret_val);
	return ret_val;
}


bool ParamTransformSeq::is_one_to_one() const
{
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2model.begin(), e=tranSeq_ctl2model.end();
		b!=e; ++b) {
			if((*b)->is_one_to_one() == false) {
				return false;
			}
	}
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2derivative.begin(), e=tranSeq_ctl2derivative.end();
		b!=e; ++b) {
			if((*b)->is_one_to_one() == false) {
				return false;
			}
	}
	for (vector<Transformation*>::const_iterator b=tranSeq_derivative2numeric.begin(), e=tranSeq_derivative2numeric.end();
		b!=e; ++b) {
			if((*b)->is_one_to_one() == false) {
				return false;
			}
	}
	return true;
}

Transformation* ParamTransformSeq::get_transformation(const string &name)
{
	Transformation *t_ptr = 0;
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2model.begin(), e=tranSeq_ctl2model.end();
		b!=e; ++b) {
			if( (*b)->get_name() == name) {
				t_ptr = (*b);
			}
	}
	if (!t_ptr) {
		for (vector<Transformation*>::const_iterator b=tranSeq_ctl2derivative.begin(), e=tranSeq_ctl2derivative.end();
			b!=e; ++b) {
				if( (*b)->get_name() == name) {
					t_ptr = (*b);
			}	
		}
	}
	if (!t_ptr) {
		for (vector<Transformation*>::const_iterator b=tranSeq_derivative2numeric.begin(), e=tranSeq_derivative2numeric.end();
			b!=e; ++b) {
				if( (*b)->get_name() == name) {
					t_ptr = (*b);
			}	
		}
	}
	return t_ptr;
}

void ParamTransformSeq::print(ostream &os) const
{
	os << "ParamTransformSeq name = " << name << endl; 
	os << "Control file to model transformations" << endl;
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2model.begin(), e=tranSeq_ctl2model.end();
		b!=e; ++b) {
			(*b)->print(os);
	}
	os << "Control file to derivative transformations" << endl;
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2derivative.begin(), e=tranSeq_ctl2derivative.end();
		b!=e; ++b) {
			(*b)->print(os);
	}
	os << "Derivative to numeric transformations" << endl;
	for (vector<Transformation*>::const_iterator b=tranSeq_derivative2numeric.begin(), e=tranSeq_derivative2numeric.end();
		b!=e; ++b) {
			(*b)->print(os);
	}
}

void ParamTransformSeq::add_custom_tran_seq(const std::string &name,  const vector<Transformation*> &tran_seq)
{
	custom_tran_seq[name] = tran_seq;
}

const vector<Transformation*> ParamTransformSeq::get_custom_tran_seq(const string &name) const
{
	auto iter = custom_tran_seq.find(name);
	assert(iter != custom_tran_seq.end());
	return iter->second;
}

void ParamTransformSeq::custom_tran_seq_forward_ip(const std::string &name, Parameters &data) const
{
	auto iter= custom_tran_seq.find(name);
	assert(iter != custom_tran_seq.end());
	const vector<Transformation*> &tran_vec = iter->second;
	
	for(const auto &itran : tran_vec)
	{
		itran->forward(data);
	}

}

ostream& operator<< (ostream &os, const ParamTransformSeq& val)
{
	val.print(os);
	return os;
}

