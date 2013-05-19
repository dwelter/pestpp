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
#include <vector>
#include <map>
#include "utilities.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <utility>
#include "ModelRunPP.h"
#include "Transformable.h"
#include "system_variables.h"
#include "pest_error.h"
#include "Transformation.h"

using namespace std;
using namespace pest_utils;

bool ModelRun::Compare::operator()(ModelRun &run1, ModelRun &run2)
{
	bool ret_val = false;
	if (run1.frozen_ctl_pars != run2.frozen_ctl_pars)
	{
		ret_val = false;
	}
	if (field == NUMERIC_PAR) {
		ret_val = run1.get_numeric_pars().get_rec(name) < run2.get_numeric_pars().get_rec(name);
	}
	else if(field == MODEL_PAR) {
		ret_val = run1.get_model_pars().get_rec(name) < run2.get_model_pars().get_rec(name);
	}
	else if(field == CTL_PAR) {
		ret_val = run1.get_ctl_pars().get_rec(name) < run2.get_ctl_pars().get_rec(name);
	}
	else if (field == SIM_OBS) {
		ret_val = run1.get_obs().get_rec(name) < run2.get_obs().get_rec(name);
	}
	else {
		assert(false);
	}
	return ret_val;
}


ModelRun::ModelRun(const ObjectiveFunc *_obj_func_ptr, const ParamTransformSeq &_par_tran, const Observations &_sim_obs) 
	: obj_func_ptr(_obj_func_ptr), par_tran(_par_tran), sim_obs(_sim_obs), phi_comp(), 
	  numeric_pars_is_valid(false), model_pars_is_valid(false),
	obs_is_valid(false), phi_is_valid(false)
{
}

ModelRun::ModelRun(const ObjectiveFunc *_obj_func_ptr, const ParamTransformSeq &_par_tran, const Parameters &_numeric_pars, const Observations &_obs)
		: obj_func_ptr(_obj_func_ptr), par_tran(_par_tran),
		  numeric_pars(_numeric_pars), sim_obs(_obs), phi_comp(),
		  numeric_pars_is_valid(true), model_pars_is_valid(false),
		  obs_is_valid(true), phi_is_valid(false)
{}

ModelRun& ModelRun::operator=(const ModelRun &rhs)
{
	frozen_ctl_pars = rhs.frozen_ctl_pars;
	obj_func_ptr = rhs.obj_func_ptr;
	par_tran = rhs.par_tran;
	numeric_pars = rhs.numeric_pars;
	model_pars = rhs.model_pars;
	sim_obs = rhs.sim_obs; 
	phi_comp = rhs.phi_comp;
	numeric_pars_is_valid = rhs.numeric_pars_is_valid;
	model_pars_is_valid = rhs.model_pars_is_valid;
	obs_is_valid = rhs.obs_is_valid;
	phi_is_valid = rhs.phi_is_valid;
	return *this;
}


const Parameters& ModelRun::get_frozen_ctl_pars() const
{
	return frozen_ctl_pars;
}

void ModelRun::set_frozen_ctl_parameters(const Parameters &frz_pars)
{
	frozen_ctl_pars = frz_pars;
}

void ModelRun::add_frozen_ctl_parameters(const Parameters &frz_pars)
{
	for (auto &ipar : frz_pars)
	{
		frozen_ctl_pars[ipar.first] = ipar.second;
	}
}





void ModelRun::set_numeric_parameters(const Parameters &parameters)
{
	numeric_pars = parameters;
	numeric_pars_is_valid = true;
	model_pars_is_valid = false;
	obs_is_valid = false;
	phi_is_valid = false;
}

void ModelRun::set_ctl_parameters(const Parameters &pars)
{
	set_numeric_parameters(par_tran.ctl2numeric_cp(pars));
}

void ModelRun::set_model_parameters(const Parameters &parameters)
{
	model_pars = parameters;
	model_pars_is_valid = true;
	numeric_pars_is_valid = false;
	obs_is_valid = false;
	phi_is_valid = false;
}

void ModelRun::update(Parameters &model_pars, Observations &obs, PAR_UPDATE update_type)
{
	//Must set parameters before observations
	if(update_type == FORCE_PAR_UPDATE || get_par_tran().is_one_to_one())
	{
		// transform to numeric parameters
		Parameters numeric_pars = get_par_tran().model2numeric_cp(model_pars);
		set_numeric_parameters(numeric_pars);
	}

	// Process Observations
	set_observations(obs);
}

void ModelRun::set_observations(const Observations &observations)
{
	sim_obs = observations;
	obs_is_valid = true;
	phi_is_valid = false;
}

const Parameters &ModelRun::get_numeric_pars()
{
	if (numeric_pars_is_valid) {
	}
	else if (model_pars_is_valid && par_tran.is_one_to_one()) {
		numeric_pars = par_tran.model2numeric_cp(model_pars);
		numeric_pars_is_valid = true;
	}
	else if (model_pars_is_valid && !par_tran.is_one_to_one())
	{
		throw PestError("ModelRun::get_numeric_pars() - can not return numeric parameters.  Transformation is not one-to-one so numeric parameters can not be calculated from model parameters");
	}
	else {
		throw PestError("ModelRun::get_numeric_pars() - can not return numeric parameters.  Both model parameters and numeric parameters are undefined");
	}
	return numeric_pars;
}

Parameters ModelRun::get_ctl_pars()
{
	if (model_pars_is_valid) {
		return par_tran.model2ctl_cp(model_pars);
	}
	else if (numeric_pars_is_valid) {
 		return par_tran.numeric2ctl_cp(numeric_pars);
	}
	else {
		throw PestError("ModelRun::get_ctl_pars() - can not return control parameters.  Both model parameters and numeric parameters are undefined");
	}
	return Parameters();
}

const Parameters &ModelRun::get_model_pars() 
{
	if (model_pars_is_valid) {
	}
	else if (numeric_pars_is_valid) {
		model_pars = par_tran.numeric2model_cp(numeric_pars);
		model_pars_is_valid = true;
	}
	else {
		throw PestError("ModelRun::get_numeric_pars() - can not return model parameters.  Both model parameters and numeric parameters are undefined");
	}
	return model_pars;
}

const Observations &ModelRun::get_obs() const
{
	if( !obs_is_valid) {
		throw PestError("ModelRun::get_obs() - observations is invalid");
	}
	return sim_obs;
}

Observations ModelRun::get_obs_template() const
{
	Observations ret_val(sim_obs);
	for (Observations::iterator b=ret_val.begin(), e=ret_val.end();
		b!=e; ++b) {
			b->second = Observations::no_data;
	}
	return ret_val;
}

double ModelRun::get_phi(double regul_weight)
{
	PhiComponets phi_comp_temp = get_phi_comp();
	return phi_comp_temp.meas + regul_weight * phi_comp_temp.regul;
}


PhiComponets ModelRun::get_phi_comp()
{
	if( !obs_is_valid) {
		throw PestError("ModelRun::get_phi() - Simulated observations are invalid.  Can not calculate phi.");
	}
	if(!phi_is_valid) {
		phi_comp = obj_func_ptr->get_phi_comp(sim_obs, get_ctl_pars());
		phi_is_valid = true;
	}
	return phi_comp;
}


vector<double> ModelRun::get_residuals_vec(const vector<string> &obs_names)
{
	return obj_func_ptr->get_residuals_vec(get_obs(), get_ctl_pars(), obs_names);
}

void ModelRun::phi_report(ostream &os) 
{
	if( !obs_is_valid) {
		throw PestError("ModelRun::phi_report() - Simulated observations are invalid.  Can not produce phi report.");
	}
	phi_comp = obj_func_ptr->phi_report(os, sim_obs, get_ctl_pars());
	phi_is_valid = true;
}

void ModelRun::full_report(ostream &os) 
{
	if( !obs_is_valid) {
		throw PestError("ModelRun::full_report() - Simulated observations are invalid.  Can not produce full report.");
	}
	phi_comp = obj_func_ptr->full_report(os, sim_obs, get_ctl_pars());
	phi_is_valid = true;
}

bool ModelRun::phi_valid() const
{
	return phi_is_valid;
}

bool ModelRun::obs_valid() const
{
	return obs_is_valid;
}

bool ModelRun::numeric_pars_valid() const
{
	return numeric_pars_is_valid;
}

bool ModelRun::model_pars_valid() const
{
	return model_pars_is_valid;
}

ModelRun::~ModelRun()
{
}
