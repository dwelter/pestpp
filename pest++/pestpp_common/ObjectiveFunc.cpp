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
#include <ostream>
#include <list>
#include <iomanip>
#include "ObjectiveFunc.h"
#include "pest_data_structs.h"
#include "Transformable.h"
#include "PriorInformation.h"

using namespace std;


const PhiComponets& PhiComponets::operator=(const PhiComponets &rhs)
{
	meas = rhs.meas; 
	regul=rhs.regul;
	return * this;
}
double ObjectiveFunc::get_phi(const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg) const
{
	double phi;
	PhiComponets phi_comp;
	phi_comp = get_phi_comp(sim_obs, pars, dynamic_reg);
	phi = phi_comp.meas + phi_comp.regul;
	return phi;
}

PhiComponets ObjectiveFunc::get_phi_comp(const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg) const
{
	unordered_map<string, ObservationRec>::const_iterator info_iter;
	unordered_map<string, ObservationRec>::const_iterator info_end = obs_info_ptr->observations.end();
	Observations::const_iterator obs_iter;
	Observations::const_iterator obs_end = observations_ptr->end();


	PhiComponets phi;
	double tmp_phi;
	for(Observations::const_iterator sim_iter=sim_obs.begin(),
	sim_end=sim_obs.end(); sim_iter!=sim_end; ++sim_iter) {
		info_iter = obs_info_ptr->observations.find((*sim_iter).first);
		obs_iter = observations_ptr->find((*sim_iter).first);
		if (info_iter != info_end && obs_iter !=obs_end) {
			tmp_phi = pow(((*sim_iter).second - (*obs_iter).second), 2.0) * (*info_iter).second.weight;
			if( (*info_iter).second.is_regularization() ) {
				if (dynamic_reg.get_use_dynamic_reg()) tmp_phi *= dynamic_reg.get_weight();
				phi.regul += tmp_phi;
			}
			else {
				phi.meas += tmp_phi;
			}
		}
	}
	for(PriorInformation::const_iterator b=(*prior_info_ptr).begin(), e=(*prior_info_ptr).end();
		b!=e; ++b) {
			tmp_phi = (*b).second.calc_phi(pars);
			if( (*b).second.is_regularization() ) {
				phi.regul += tmp_phi;
			}
			else {
				phi.meas += tmp_phi;
			}
	}
	return phi;
}



map<string, double> ObjectiveFunc::get_group_phi(const Observations &sim_obs, const Parameters &pars) const
{
	map<string, double> group_phi;
	unordered_map<string, ObservationRec>::const_iterator info_iter;
	unordered_map<string, ObservationRec>::const_iterator info_end = obs_info_ptr->observations.end();
	Observations::const_iterator obs_iter;
	Observations::const_iterator obs_end = observations_ptr->end();
	double tmp_phi;
	const string *group;

	// first add all groups to group_phi
	for(unordered_map<string, ObservationGroupRec>::const_iterator grp_iter=obs_info_ptr->groups.begin(),
		grp_end=obs_info_ptr->groups.end(); grp_iter!=grp_end; ++grp_iter) {
			group = &(*grp_iter).first;
			group_phi[*group] = 0.0;
	}

	for(Observations::const_iterator sim_iter=sim_obs.begin(),
		sim_end=sim_obs.end(); sim_iter!=sim_end; ++sim_iter) {
			info_iter = (*obs_info_ptr).observations.find((*sim_iter).first);
			obs_iter = observations_ptr->find((*sim_iter).first);
			if (info_iter != info_end && obs_iter !=obs_end) {
				tmp_phi = pow(((*sim_iter).second - (*obs_iter).second),2.0) * (*info_iter).second.weight;
				group = &(*info_iter).second.group;
				group_phi[*group] += tmp_phi;
			}
	}
	for(PriorInformation::const_iterator b=(*prior_info_ptr).begin(), e=(*prior_info_ptr).end();
		b!=e; ++b) {
			group = (*b).second.get_group_ptr();
			tmp_phi = (*b).second.calc_phi(pars);
			group_phi[*group] += tmp_phi;
	}
	return group_phi;
}

PhiComponets ObjectiveFunc::phi_report(ostream &os, const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg) const
{
	map<string, double> group_phi;
	PhiComponets phi_comp = get_phi_comp(sim_obs, pars, dynamic_reg);
	double tikhonov_weight = dynamic_reg.get_weight();
	double total_phi = phi_comp.meas + phi_comp.regul;
	if (dynamic_reg.get_use_dynamic_reg())
	{
		os << "    Current regularization weight factor                      : " << tikhonov_weight << endl;
	}
	os << "    Starting phi for this iteration                     Total : " << total_phi << endl;
	group_phi = get_group_phi(sim_obs, pars);
	for (auto &igrp : group_phi)
	{
		if (ObservationGroupRec::is_regularization(igrp.first))
		{
			igrp.second *= tikhonov_weight;
		}
	}
	for(map<string, double>::const_iterator b=group_phi.begin(), e=group_phi.end();
		b!=e; ++b) {
			os << "    Contribution to phi from observation group ";
			os << setw(17) << setiosflags(ios::right) << "\"" + (*b).first + "\" : ";
			os << (*b).second << endl;
	}
	return phi_comp;
}


PhiComponets ObjectiveFunc::full_report(ostream &os, const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg) const
{
	map<string, double> group_phi;
	PhiComponets phi_comp = get_phi_comp(sim_obs, pars, dynamic_reg);
	double tikhonov_weight = dynamic_reg.get_weight();
	double total_phi = phi_comp.meas + phi_comp.regul;
	os << "    Phi                                                 Total : " << total_phi << endl;
	group_phi = get_group_phi(sim_obs, pars);
	for(map<string, double>::const_iterator b=group_phi.begin(), e=group_phi.end();
		b!=e; ++b) {
			os << "    Contribution to phi from observation group ";
			os << setw(17) << setiosflags(ios::right) << "\"" + (*b).first + "\" : ";
			os << (*b).second << endl;
	}
	os << endl;
	os << "     Parameter      Optimal"<< endl;
	os << "        Name         Value"<< endl;
	os << "    ------------  ------------" << endl;

	const string *p_name;
	double par_value;
	for( Parameters::const_iterator nb=pars.begin(), ne=pars.end();
		nb!=ne; ++nb) {
		p_name = &((*nb).first);
		par_value =  (*nb).second;
		os << left;
		os << "    " << setw(12) << *p_name
			;
		os << right;
		os << "  " << setw(12) << par_value << endl;
	}
	return phi_comp;
}


vector<double> ObjectiveFunc::get_residuals_vec(const Observations &sim_obs, const Parameters &pars, const vector<string> &obs_names) const
{
	vector<double> residuals_vec;
	residuals_vec.resize(obs_names.size(), 0.0);

	Observations::const_iterator found_obs;
	Observations::const_iterator not_found_obs=(*observations_ptr).end();
	PriorInformation::const_iterator found_prior_info;
	PriorInformation::const_iterator not_found_prior_info = prior_info_ptr->end();

	int i=0;
	for(vector<string>::const_iterator b=obs_names.begin(), e=obs_names.end(); b != e; ++b, ++i)
	{
		found_obs = observations_ptr->find(*b);
		found_prior_info = prior_info_ptr->find(*b);

		if (found_obs != not_found_obs)
		{
			residuals_vec[i] = sim_obs.get_rec(*b) - (*found_obs).second;
		}
		else if (found_prior_info != not_found_prior_info)
		{
			residuals_vec[i] = (*found_prior_info).second.calc_residual(pars);
		}
	}
	return residuals_vec;
}

const Observations *ObjectiveFunc::get_obs_ptr() const
{
	return observations_ptr;
}

const ObservationInfo* ObjectiveFunc::get_obs_info_ptr() const
{
	return obs_info_ptr;
}

const PriorInformation*  ObjectiveFunc::get_prior_info_ptr() const
{
	return prior_info_ptr;
}