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
#ifndef REGULARIZATION_H_
#define REGULARIZATION_H_

class ModelRunAbstractBase;

/** @Regularization base class

The is the base class for regularization.  It provides default behavior that
always returns a weight factor of 1.0
*/
class Regularization
{
public:
	Regularization(void);
	virtual void update_weight(){}
	virtual double get_weight(ModelRunAbstractBase &model_run) const {return 1.0;}
	virtual ~Regularization(void);
};



class RegularizationPest : public Regularization
{
public:
	RegularizationPest(double _phi_m_lim, double _phi_m_accept,
		double _frac_phi_m, double _wf_min, double _wf_max,
		double _wffac, double _wftol);
	virtual double get_weight(ModelRunAbstractBase &model_run) const;
	virtual ~RegularizationPest(void){}
protected:
	double phi_m_lim;
	double phi_m_accept;
	double frac_phi_m;
	double wf_min;
	double wf_max;
	double wffac;
	double wftol;
};


#endif //REGULARIZATIONABSTRACT_H_