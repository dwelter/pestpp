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
#ifndef TERMINATIONCONTROLLER_H_
#define TERMINATIONCONTROLLER_H_

#include <vector>
#include <fstream>
#include <iostream>

class RestartController;

class TerminationController
{
	friend RestartController;
public:
	TerminationController(int _noptmax=0, double _phiredstp=0.0, int _nphistp=0, int _nphinored=0,
		double _relparstp=0.0, int _nrelpar=0);
	TerminationController(const TerminationController &rhs) {*this = rhs;}
	bool process_iteration(double phi, double relpar);
	void set_terminate(bool _terminate_code) { terminate_code = _terminate_code; }
	bool check_last_iteration();
	void reset();
	const TerminationController& operator=(const TerminationController &rhs);
	int get_iteration_number() {return nopt_count;}
	void save_state(std::ostream &fout);
	void read_state(const std::string &line);
	bool terminate() { return terminate_code; }
	void termination_summary(std::ostream &fout);
	~TerminationController(void);
private:
	int noptmax;
	int nopt_count;
	int nphinored;
	int nphinored_count;
	int nrelpar;
	int nrelpar_count;
	unsigned nphistp;
	double phiredstp;
	double relparstp;
	std::vector<double> lowest_phi;
	bool terminate_code;
	std::string termimate_reason;
};

#endif //TERMINATIONCONTROLLER_H_