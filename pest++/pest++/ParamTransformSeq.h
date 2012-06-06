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
#ifndef ParamTransformSeq_H_
#define ParamTransformSeq_H_

/** @file
 @brief ParamTransformSeq class 
 
 This file defines the ParamTransformSeq.  This class is used to define the sequence of 
 transformations used to switch between model, control and numeric parameter representations
*/

#include <vector>
#include <string>
#include "Transformation.h"

class Parameters;
class JacobianAnalytic;

using namespace std;

/**
 @brief ParamTransformSeq class 
 
 This class is used to define the sequence of transformations used to switch
 between model, control and numeric parameter representations
*/
class ParamTransformSeq {
public:
	ParamTransformSeq(const string &_name="unnamed ParamTransformSeq") : name(_name), ctl_offset_ptr(0), ctl_scale_prt(0) {}
	ParamTransformSeq(const ParamTransformSeq &rhs);
	void push_back_ctl2model(Transformation *tr);
	void push_back_ctl2numeric(Transformation *tr);
	void push_back_ctl2numeric_frozen();
	void numeric2ctl_ip(Parameters &data) const;
	void numeric2model_ip(Parameters &data) const;
	void ctl2numeric_ip(Parameters &data) const;
	void ctl2model_ip(Parameters &data) const;
	void model2ctl_ip(Parameters &data) const;
	void model2numeric_ip(Parameters &data) const;
	Parameters numeric2ctl_cp(const Parameters &data) const;
	Parameters numeric2model_cp(const Parameters &data) const;
	Parameters ctl2numeric_cp(const Parameters &data) const;
	Parameters ctl2model_cp(const Parameters &data) const;
	Parameters model2ctl_cp(const Parameters &data) const;
	Parameters model2numeric_cp(const Parameters &data) const;
	void set_offset_ptr(TranOffset *ptr) {ctl_offset_ptr = ptr;}
	const TranOffset *get_offset_ptr() const{return ctl_offset_ptr;}
	void set_scale_ptr(TranScale *ptr) {ctl_scale_prt = ptr;}
	const TranScale *get_scale_ptr() const {return ctl_scale_prt;}
	void set_log10_ptr(TranLog10 *ptr) {ctl_log10_ptr = ptr;}
	const TranLog10 *get_log10_ptr() const {return ctl_log10_ptr;}
	TranFrozen *get_frozen_ptr();
	void analytic_derivative_to_numeric(JacobianAnalytic &jacobian, const Parameters & model_pars) const;
	bool is_one_to_one() const;
	void print(ostream &os) const;
	ParamTransformSeq &operator+=(const ParamTransformSeq &rhs);
	const ParamTransformSeq &operator=(const ParamTransformSeq &rhs);
	virtual ~ParamTransformSeq();
private:
	vector<Transformation*> tranSeq_ctl2model;
	vector<Transformation*> tranSeq_ctl2numeric;
	TranOffset *ctl_offset_ptr;
	TranScale *ctl_scale_prt;
	TranLog10 *ctl_log10_ptr;
	TranFrozen ctl_frozen;
	string name;
};

template <typename Parameters>
ostream& operator<< (ostream &os, const ParamTransformSeq& val);



#endif /* ParamTransformSeq_H_ */
