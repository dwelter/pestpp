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
#include "TerminationController.h"
#include <vector>
#include <algorithm>

using namespace std;

TerminationController::TerminationController(int _noptmax, double _phiredstp,
	int _nphistp, int _nphinoger, double _relparstp, int _nrelpar) 
	: noptmax(_noptmax), phiredstp(_phiredstp), nphistp(_nphistp),
	nphinoger(_nphinoger), relparstp(_relparstp), 
	nrelpar(_nrelpar), nphinoger_count(0), nrelpar_count(0), nopt_count(0)
{
}

void TerminationController::reset() {
	nopt_count = 0;
	nphinoger_count = 0;
	lowest_phi.clear();
}

//TerminationController& TerminationController::operator+=(const TerminationController &rhs)
//{
//	if (nphistp!=rhs.nphistp !!
//		phiredstp!=rhs.phiredstp
//		relparstp!=rhs.relparstp) {
//
//	}
//	// add number of iterations
//	nopt_count += rhs.nopt_count;
//	// keep track of NPHISTP lowest phi's
//	lowest_phi.insert(lowest_phi.end(), rhs.lowest_phi.begin(), rhs.lowest_phi.end());
//	sort(lowest_phi.begin(), lowest_phi.end());
//	if (lowest_phi.size() > nphistp) {
//		lowest_phi.resize(nphistp);
//	}
//	if (rhs.nrelpar_count == rhs.nopt_count) {
//		nrelpar_count += rhs.nrelpar_count;
//	}
//	else {
//		nrelpar_count = rhs.nrelpar_count;
//	}
//}

bool TerminationController::process_iteration(double phi, double relpar)
{
	++nopt_count;
	// keep track of number of iterations since lowest phi
	if (lowest_phi.size() == 0 || phi <= lowest_phi.front()) {
		nphinoger_count = 0;
	}
	else {
		++nphinoger_count;
	}

	// keep track of NPHISTP lowest phi's
	if (lowest_phi.size() < nphistp) {
		lowest_phi.push_back(phi);
	}
	else if (phi < lowest_phi.back()) {
		lowest_phi.back() = phi;
	}
	sort(lowest_phi.begin(), lowest_phi.end());

	// Check maximum relative parameter change
	if (relpar >= relparstp) {
		++nrelpar_count;
	}
	else {
		nrelpar_count = 0;
	}
	return check_last_iteration();
}
bool TerminationController::check_last_iteration()
{
	double min_phi_diff = lowest_phi.back() - lowest_phi.front();
	// Impose Termination Criteria
	if (nphinoger_count > nphinoger)
		return true;
	if (lowest_phi.size() >= nphistp && min_phi_diff <= phiredstp)
		return true;
	if (nrelpar_count > nrelpar)
		return true;
	if (nopt_count >= noptmax)
		return true;
	return false;
}

const TerminationController& TerminationController::operator=(const TerminationController &rhs)
{
	noptmax = rhs.noptmax;
	nphinoger = rhs.nphinoger;
	nphinoger_count = rhs.nphinoger_count;
	nrelpar = rhs.nrelpar;
	nrelpar_count = rhs.nrelpar_count;
	nphistp = rhs.nphistp;
	phiredstp = rhs.phiredstp;
	relparstp = rhs.relparstp;
	lowest_phi = rhs.lowest_phi;
	return *this;
}







TerminationController::~TerminationController(void)
{
}
