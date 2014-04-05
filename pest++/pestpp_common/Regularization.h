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

class ModelRun;

class DynamicRegularization
{
public:
	DynamicRegularization(double _phi_m_lim = 0, double _phi_m_accept = 0,
		double _frac_phi_m = 1, double _wf_min = 1e-10, double _wf_max = 1e10,
		double _wffac = 0, double _wftol = 1000, double _tikhonov_weight = 1);
	virtual double get_weight() const;
	virtual void set_weight(double _tikhonov_weight) {tikhonov_weight = _tikhonov_weight;}
	virtual ~DynamicRegularization(void){}
protected:
	double phi_m_lim;
	double phi_m_accept;
	double frac_phi_m;
	double wf_min;
	double wf_max;
	double wffac;
	double wftol;
	double tikhonov_weight;
};


#endif //REGULARIZATIONABSTRACT_H_