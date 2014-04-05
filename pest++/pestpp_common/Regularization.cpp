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
#include "Regularization.h"
#include "ObjectiveFunc.h"
#include "ModelRunPP.h"


DynamicRegularization::DynamicRegularization(double _phi_m_lim,
	double _phi_m_accept, double _frac_phi_m, double _wf_min, double _wf_max,
	double _wffac, double _wftol, double _tikhonov_weight)
	: phi_m_lim(_phi_m_lim), phi_m_accept(_phi_m_accept), frac_phi_m(_frac_phi_m),
	wf_min(_wf_min), wf_max(_wf_max), wffac(_wffac), wftol(_wftol), 
	tikhonov_weight(_tikhonov_weight)
{
}

double DynamicRegularization::get_weight() const
{
	return tikhonov_weight;
}