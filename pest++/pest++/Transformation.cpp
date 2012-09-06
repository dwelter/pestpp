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
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <lapackpp.h>
#include <blas3pp.h>
#include <blas2pp.h>
#include <cassert>
#include "Transformation.h"
#include "Transformable.h"
#include "lapack_tools.h"
#include "Jacobian.h"
#include "QSqrtMatrix.h"
#include "pest_data_structs.h"
using namespace std;


///////////////// TranMapBase Methods /////////////////
void TranMapBase::insert(const string &item_name, double item_value)
{
	items[item_name] = item_value;
}

void TranMapBase::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranMapBase)" << endl; 
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  value = " << (*b).second << endl;   
	}
}

pair<bool, double> TranMapBase::get_value(const string &name) const
{
	pair<bool, double> ret_val(false, 0.0);
	map<string, double>::const_iterator it;
	
	it = items.find(name);
	if (it !=items.end()) {
		ret_val = pair<bool, double>(true, (*it).second);
	}
	return ret_val;
}

///////////////// TranSetBase Methods /////////////////
void TranSetBase::insert(const string &item_name)
{
	items.insert(item_name);
}

void TranSetBase::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranSetBase)" << endl;
	for (set<string>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << *b << endl;   
	}
}

bool TranSetBase::has_value(const string &name) const
{
	bool ret_val = false;
	set<string>::const_iterator it;
	
	it = items.find(name);
	if (it !=items.end()) {
		ret_val = true;
	}
	return ret_val;
}

///////////////// TranOffset Methods /////////////////
void TranOffset::forward(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second += (*b).second;
		}
	}
}

void TranOffset::reverse(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second -= b->second;
		}
	}
}

void TranOffset::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranOffset)" << endl; 
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  offset value = " << (*b).second << endl;   
	}
}



///////////////// TranScale Methods /////////////////
void TranScale::forward(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second *= b->second;
		}
	}
}


void TranScale::reverse(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second /= b->second;
		}
	}
}

void TranScale::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranScale)" << endl; 
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  scale value = " << (*b).second << endl;   
	}
}

///////////////// TranLog10 Methods /////////////////
void TranLog10::forward(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (set<string>::const_iterator b=items.begin(), e=items.end(); b!=e; ++b) 
	{
		data_iter = data.find(*b);
		if (data_iter != data_end)
		{
			(*data_iter).second = log10((*data_iter).second);
		}
	}
}

void TranLog10::reverse(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (set<string>::const_iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		data_iter = data.find(*b);
		if (data_iter != data_end)
		{
			(*data_iter).second = pow(10.0, (*data_iter).second);
		}
	}
}

void TranLog10::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranLog10)" << endl;
	for (set<string>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << *b << endl;   
	}
}

///////////////// TranInvLog10 Methods /////////////////
void TranInvLog10::forward(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (set<string>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(*b);
		if (data_iter != data_end)
		{
			(*data_iter).second = pow(10.0,(*data_iter).second);
		}
	}
}

void TranInvLog10::reverse(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (set<string>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(*b);
		if (data_iter != data_end)
		{
			(*data_iter).second = log10((*data_iter).second);
		}
	}
}

void TranInvLog10::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranInvLog10)" << endl;
	for (set<string>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << *b << endl;   
	}
}


///////////////// TranFixed Methods /////////////////
void TranFixed::forward(Transformable &data)
{
	for (map<string,double>::iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		data.erase(b->first);
	}
}

void TranFixed::reverse(Transformable &data)
{
	for (map<string,double>::iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		data.insert(b->first, b->second);
	}
}

void TranFixed::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranFixed)" << endl; 
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  imposed value = " << (*b).second << endl;   
	}
}


void TranTied::insert(const string &item_name, const pair<string, double> &item_value)
{
	items[item_name] = item_value;
}

void TranTied::reverse(Transformable &data)
{
	string const *base_name;
	double *factor;
	Transformable::iterator base_iter;
	for (map<string, pair_string_double>::iterator b=items.begin(), e=items.end();
		b!=e; ++b)
	{
		base_name = &(b->second.first);
		factor = &(b->second.second);
		base_iter = data.find(*base_name);
		if (base_iter != data.end())
		{
			data.insert(b->first, (*base_iter).second * (*factor));
		}
	}
}

void TranTied::forward(Transformable &data)
{
	for (map<string, pair_string_double>::iterator ii=items.begin(); ii!=items.end(); ++ii)
	{
		data.erase(ii->first);
	}
}
void TranTied::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranTied)" << endl; 
	for (map<string, pair_string_double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << "   tied to \"" << (*b).second.first << 
				"\" with factor " <<  (*b).second.second<<  endl;   
	}
}



void TranFrozen::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranFrozen)" << endl; 
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  frozen value = " << (*b).second << endl;   
	}
}

#ifndef NO_LAPACKPP
void TranSVD::update(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const Parameters &base_numeric_pars,
		int maxsing, double eigthresh, const vector<string> &par_names, const vector<string> &obs_names)
{
	stringstream sup_name;
	super_parameter_names.clear();
	base_parameter_names = par_names;

	LaGenMatDouble SqrtQ_J = Q_sqrt * jacobian.get_matrix(par_names, obs_names);
	Sigma.resize(min(SqrtQ_J.size(0), SqrtQ_J.size(1)));
	U.resize(SqrtQ_J.size(0), SqrtQ_J.size(0));
	Vt.resize(SqrtQ_J.size(1), SqrtQ_J.size(1));
	LaSVD_IP(SqrtQ_J, Sigma, U, Vt);
	// calculate the number of singluar values above the threshold
	LaGenMatDouble svd_mat_inv = SVD_inv(U, Sigma, Vt, maxsing, eigthresh, n_sing_val);
	super_parameter_names.clear();
	for(int i=0; i<n_sing_val; ++i) {
		sup_name.str("");
		sup_name << "sup_";
		sup_name << i+1;
		super_parameter_names.push_back(sup_name.str());
	}
	init_base_numeric_parameters = base_numeric_pars;
}





void TranSVD::reverse(Transformable &data)
{
	// Transform super-parameters to base parameters
	assert(Vt.cols() == base_parameter_names.size());
	int n_base = Vt.cols();
	vector<double> super_par_vec = data.get_vector(super_parameter_names);
	vector<double>::iterator it;
	for (it=super_par_vec.begin(); it!=super_par_vec.end(); ++it)
	{
		(*it) -= 10.0;
	}
	Transformable ret_base_pars;
	LaVectorDouble delta_base_mat(Vt.cols());
	Blas_Mat_Trans_Mat_Mult(Vt(LaIndex(0,n_sing_val-1), LaIndex(0, Vt.rows()-1)), stlvec2LaVec(super_par_vec), delta_base_mat, 1.0, 0.0);
	for (int i=0; i<n_base; ++i) {
		ret_base_pars.insert(base_parameter_names[i], delta_base_mat(i) + init_base_numeric_parameters.get_rec(base_parameter_names[i]))
;	}
	data = ret_base_pars;
}

void TranSVD::forward(Transformable &data)
{
	//Transform base parameters to super-parameters
	Transformable super_pars;
	LaVectorDouble value(Vt.cols());

	Transformable delta_data = data - init_base_numeric_parameters;
	Blas_Mat_Vec_Mult(Vt, stlvec2LaVec(delta_data.get_vector(base_parameter_names)), value, 1.0, 0.0);
	for (int i=0; i<n_sing_val; ++i) {
		super_pars.insert(super_parameter_names[i], value(i)+10.0);
	}
	data = super_pars;
}

ParameterGroupInfo TranSVD::build_par_group_info(const ParameterGroupInfo &base_pg_info)
{
	double derinc_sup;
	double derinc_par;
	int max_col;
	double max_val;
	ParameterGroupInfo pg_info;
	stringstream grp_name;

	for (int i_sup=0, n_sup=super_parameter_names.size(); i_sup < n_sup; ++i_sup)
	{
		get_LaGenMatDouble_row_abs_max(Vt, i_sup, &max_col, &max_val);
		//dew clean up
		derinc_par = base_pg_info.get_group_rec_ptr(base_parameter_names[max_col])->derinc;
		derinc_sup = derinc_par / (abs(max_val)*Sigma(i_sup));
		derinc_sup = 0.015;
		grp_name.str("");
		grp_name << "g_" << super_parameter_names[i_sup];
		ParameterGroupRec sup_rec(grp_name.str(), "ABSOLUTE", derinc_sup, 0.0, "SWITCH", 2.0, "PARABOLIC");
		//add new group
		pg_info.insert_group(grp_name.str(), sup_rec);
		// connect super parameter to new group
		pg_info.insert_parameter_link(super_parameter_names[i_sup], grp_name.str());
	}
	return pg_info;
}

void TranSVD::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranSVD)" << endl; 
	os << "  Singular Values = " << Sigma << endl;
}

#endif

void TranNormalize::forward(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,NormData>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second += b->second.offset;
			(*data_iter).second *= b->second.scale;
		}
	}
}


void TranNormalize::reverse(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,NormData>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second -= b->second.offset;
			(*data_iter).second /= b->second.scale;
		}
	}
}

void TranNormalize::insert(const string &item_name, double _offset, double _scale)
{
	items[item_name] = NormData(_offset, _scale);
}

void TranNormalize::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranNormalize)" << endl; 
	for (map<string,NormData>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  scale value = " << (*b).second.scale
				<< ";  offset value = " << (*b).second.offset <<endl;   
	}
}