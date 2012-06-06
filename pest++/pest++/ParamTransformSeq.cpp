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
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "Transformable.h"
#include "Jacobian.h"

using namespace std;

ParamTransformSeq::ParamTransformSeq(const ParamTransformSeq &rhs)
{
	*this = rhs;
}

void ParamTransformSeq::push_back_ctl2model(Transformation *tr)
{
	tranSeq_ctl2model.push_back(tr);
}

void ParamTransformSeq::push_back_ctl2numeric(Transformation *tr)
{
	tranSeq_ctl2numeric.push_back(tr);
}

void ParamTransformSeq::push_back_ctl2numeric_frozen()
{
	tranSeq_ctl2numeric.push_back(&ctl_frozen);
}

ParamTransformSeq::~ParamTransformSeq()
{
	//for(vector<Transformation*>::iterator i = tranSeq_ctl2model.begin(),
	//	e=tranSeq_ctl2model.end(); i != e; ++i)
	//{
	//	delete *i;
	//}
	//for(vector<Transformation*>::iterator i = tranSeq_ctl2numeric.begin(),
	//	e=tranSeq_ctl2numeric.end(); i != e; ++i)
	//{
	//	delete *i;
	//}
}

TranFrozen *ParamTransformSeq::get_frozen_ptr() {return &ctl_frozen;}

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

void ParamTransformSeq::ctl2numeric_ip(Parameters &data) const
{
	vector<Transformation*>::const_iterator iter, e;

	for(iter = tranSeq_ctl2numeric.begin(), e = tranSeq_ctl2numeric.end();
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

void ParamTransformSeq::numeric2ctl_ip(Parameters &data) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for(iter = tranSeq_ctl2numeric.rbegin(), e = tranSeq_ctl2numeric.rend();
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

bool ParamTransformSeq::is_one_to_one() const
{
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2model.begin(), e=tranSeq_ctl2model.end();
		b!=e; ++b) {
			if((*b)->is_one_to_one() == false) {
				return false;
			}
	}
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2numeric.begin(), e=tranSeq_ctl2numeric.end();
		b!=e; ++b) {
			if((*b)->is_one_to_one() == false) {
				return false;
			}
	}
	return true;
}

ParamTransformSeq &ParamTransformSeq::operator+=(const ParamTransformSeq &rhs)
{
	for (vector<Transformation*>::const_iterator b=rhs.tranSeq_ctl2model.begin(), e=rhs.tranSeq_ctl2model.end();
		b!=e; ++b) {
			tranSeq_ctl2model.push_back((*b)->clone());
	}
		for (vector<Transformation*>::const_iterator b=rhs.tranSeq_ctl2numeric.begin(), e=rhs.tranSeq_ctl2numeric.end();
		b!=e; ++b) {
			tranSeq_ctl2model.push_back((*b)->clone());
	}
	return *this;
}

const ParamTransformSeq& ParamTransformSeq::operator=(const ParamTransformSeq &rhs)
{
	name = "copy";
	Transformation *tran;
	tranSeq_ctl2model.clear();
	for (vector<Transformation*>::const_iterator b=rhs.tranSeq_ctl2model.begin(), e=rhs.tranSeq_ctl2model.end();
		b!=e; ++b) {
			tran = (*b)->clone();
			tranSeq_ctl2model.push_back(tran);
			if (*b == rhs.ctl_offset_ptr) ctl_offset_ptr = dynamic_cast<TranOffset*>(tran);
			if (*b == rhs.ctl_scale_prt) ctl_scale_prt = dynamic_cast<TranScale*>(tran);
	}
	ctl_frozen = rhs.ctl_frozen;
	tranSeq_ctl2numeric.clear();
	for (vector<Transformation*>::const_iterator b=rhs.tranSeq_ctl2numeric.begin(), e=rhs.tranSeq_ctl2numeric.end();
		b!=e; ++b) {
			if ((*b) == &(rhs.ctl_frozen))
			{
				tran = &ctl_frozen;
			}
			else {
				tran = (*b)->clone();
			}
			tranSeq_ctl2numeric.push_back(tran);
			if (*b == rhs.ctl_log10_ptr) ctl_log10_ptr = dynamic_cast<TranLog10*>(tran);
	}
	return *this;
}


void ParamTransformSeq::print(ostream &os) const
{
	os << "ParamTransformSeq name = " << name << endl; 
	os << "Control file to model transformations" << endl;
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2model.begin(), e=tranSeq_ctl2model.end();
		b!=e; ++b) {
			(*b)->print(os);
	}
	os << "Control file to numeric transformations" << endl;
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2numeric.begin(), e=tranSeq_ctl2numeric.end();
		b!=e; ++b) {
			(*b)->print(os);
	}
}


ostream& operator<< (ostream &os, const ParamTransformSeq& val)
{
	val.print(os);
	return os;
}

//void ParamTransformSeq::analytic_derivative_to_numeric(JacobianAnalytic &jac, const Parameters & model_pars) const
//{
//	const string *p_name;
//	const Parameters &pars = jac.get_base_numeric_parameters();
//	const Observations &obs = jac.get_base_sim_observations();
//	LaGenMatDouble &matrix = jac.get_matrix(pars.get_keys(), obs.get_keys());
//	pair<bool, double> val_pair;
//	double scale, offset;
//	double par_value;
//	double factor;
//	bool has_log;
//	int irow, icol;
//
//	icol = 0;
//	for (Parameters::const_iterator b=pars.begin(), e=pars.end();
//		b != e; ++b, ++icol)
//	{
//		p_name = &(*b).first;
//		par_value = model_pars.get_rec(*p_name);
//		val_pair = ctl_offset_ptr->get_value(*p_name);
//		if (val_pair.first == true) {
//			offset = val_pair.second;
//		}
//		else {
//			offset = 0.0; 
//		}
//		val_pair = ctl_scale_prt->get_value(*p_name);
//		if (val_pair.first == true) {
//			scale = val_pair.second;
//		}
//		else {
//			scale = 1.0;
//		}
//		has_log = ctl_log10_ptr->has_value(*p_name);
//
//		irow = 0;
//		if (has_log = true) {
//			factor = (par_value - offset) *  log(10.0);
//			for (Observations::const_iterator bo=obs.begin(), eo=obs.end();
//			bo != eo; ++bo, ++irow)
//			{
//				matrix(irow,icol)  *=  factor;
//			}
//		}
//		else if (scale != 1) {
//			factor =  par_value / scale;
//			for (Observations::const_iterator bo=obs.begin(), eo=obs.end();
//			bo != eo; ++bo, ++irow)
//			{
//				matrix(irow,icol)  *= factor;
//			}
//		}
//	}
//}