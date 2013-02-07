/*  
    � Copyright 2012, David Welter
    
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
#include "PriorInformation.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "utilities.h"
#include "Transformable.h"

using namespace::std;
using namespace pest_utils;

PIAtom::PIAtom(const string &_par_name, bool _log_trans, double _factor)
	: par_name(_par_name), log_trans(_log_trans), factor(_factor)
{

}

PriorInformationRec::PriorInformationRec(double _pival, double _weight,
	const string &_group, const vector<PIAtom> _pi_atoms)
	: pival(_pival), weight(_weight), group(_group), pi_atoms(_pi_atoms)
{
}

PriorInformationRec::PriorInformationRec(const PriorInformationRec &rhs)
{
	*this = rhs;
}

const PriorInformationRec& PriorInformationRec::operator=(const PriorInformationRec &rhs)
{
	pi_atoms = rhs.pi_atoms;
	pival = rhs.pival;
	weight = rhs.weight;
	group = rhs.group;
	return *this;
}


double PriorInformationRec::calc_residual(const Parameters &pars) const
{
	double sim_value = 0;
	double par_value;
	Parameters::const_iterator p_iter;
	for (vector<PIAtom>::const_iterator b=pi_atoms.begin(), e=pi_atoms.end();
		b!=e; ++b) {
		par_value = pars.get_rec((*b).par_name);
		if ((*b).log_trans) {
			sim_value += (*b).factor * log10(par_value);
		}
		else {
			sim_value += (*b).factor * par_value;
		}
	}
	return (sim_value - pival);
}

double PriorInformationRec::calc_phi(const Parameters &pars) const
{
	return pow(calc_residual(pars), 2) * weight;
}


bool PriorInformationRec::is_regularization() const
{
	bool is_reg = false;
	int found = upper_cp(group).find("REGUL");
	if (found == 0) is_reg=true;
	return is_reg;
}




PriorInformationRec::~PriorInformationRec(void)
{
}


ostream& operator<< (ostream& out, const PriorInformationRec &rhs)
{
	struct Local {
		static string GetSign(double val) {
			if (val < 0) return "-";
			return "+";
		}
	};

	for(vector<PIAtom>::const_iterator b=rhs.pi_atoms.begin(), e=rhs.pi_atoms.end();
		b!=e; ++b) {
			out << " " << Local::GetSign((*b).factor);
			out << " " << abs((*b).factor);
			out << " * ";
			if ((*b).log_trans)
				out << "LOG(" << (*b).par_name << ")";
			else
				out << (*b).par_name;
			out << " = " << rhs.pival;
			out << "   " << rhs.weight;
			out << "   " << rhs.group;
			out << endl;
	}
	return out;
}

void PriorInformation::AddRecord(const string &pi_line)
{
	string par_name;
	double pifac;
	double pival;
	double weight;
	bool log_trans;
	string group;
	bool plus_sign;
	vector<PIAtom> pi_atoms;
	vector<string> tokens;
	tokenize(upper_cp(strip_cp(pi_line)), tokens);
	int n_tokens = tokens.size();
	group = tokens[n_tokens-1];
	convert_ip(tokens[n_tokens-2], weight);
	weight *= weight;
	convert_ip(tokens[n_tokens-3], pival);

	string prior_info_name = tokens[0];
	// process prior information equation
	plus_sign = true;
	for (int i=1; i< n_tokens-4; i+=4) 
	{
		convert_ip(tokens[i], pifac);
		// token i+1 = "*"
		par_name = tokens[i+2];
		if (par_name.find("LOG") == 0) // parameter name
		{ 
			log_trans = true;
			par_name = par_name.substr(4, par_name.size() - 5); // remove "log(" from front and ")" from rear of tokens[i+1]
		}
		else {
			log_trans = false;
		}
		if (!plus_sign)
			pifac *= -1.0;

		// add term to vector storing prior information equation
		pi_atoms.push_back(PIAtom(par_name, log_trans, pifac));


		if (tokens[i+3] == "+") {
			plus_sign = true;
		}
		else {
			plus_sign = false;
		}
	}
	prior_info_map[prior_info_name] = PriorInformationRec(pival, weight, group, pi_atoms);
}

vector<string> PriorInformation::get_keys() const
{
	vector<string> ret_val;

	for(std::map<std::string, PriorInformationRec>::const_iterator b=prior_info_map.begin(), e=prior_info_map.end();
		b!=e; ++b) 
	{
	ret_val.push_back((*b).first);
	}
	return ret_val;
}
