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
#include <fstream>
#include <iostream>
#include <lapackpp.h>
#include <iomanip>
#include <map>
#include "SVDASolver.h"
#include "ModelRunPP.h"
#include "QSqrtMatrix.h"
#include "lapack_tools.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "PriorInformation.h"

using namespace std;
using namespace pest_utils;
SVDASolver::SVDASolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_par_group_info_ptr, const ParameterInfo *_ctl_par_info_ptr,
		const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, const Parameters &_parameters, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, const Regularization *_regul_scheme)
		: SVDSolver(_ctl_info, _svd_info, _par_group_info_ptr, _ctl_par_info_ptr, _obs_info, 
		_file_manager, _observations, _obj_func, _par_transform, _parameters, _prior_info_ptr, _jacobian, 
		_regul_scheme, "super parameter solution")
{
}

map<string,double> SVDASolver::freeze_parameters(Parameters &cur_numeric_pars, const LaVectorDouble &svd_update_uvec, double svd_update_norm, const LaVectorDouble &grad_update_uvec, bool use_descent)
{
	//This routine has the side effect that it modifies par_transform transformtion
	// Test SVD parameter upgrade
	Parameters upgrade_pars = cur_numeric_pars;
	Parameters tmp_parameters;
	int ip = 0;
	map<string,double> tmp_svd;
	map<string,double> tmp_descent;
	for(Parameters::iterator b=upgrade_pars.begin(), e=upgrade_pars.end(); b!=e; ++b, ++ip)
	{
		b->second += svd_update_uvec(ip)*svd_update_norm;
	}
	tmp_svd = limit_parameters_ip(cur_numeric_pars, upgrade_pars);

	//// Test Greatest descent parameter upgrade
	//upgrade_pars = cur_solution.numeric_pars;
	//ip = 0;
	//for(Parameters::iterator b=upgrade_pars.begin(), e=upgrade_pars.end(); b!=e; ++b, ++ip)
	//{
	//	b->second += grad_update_uvec(ip)*svd_update_norm;
	//}
	//tmp_descent =  LimitParameters_one2one_IP(cur_solution.numeric_pars, upgrade_pars);
	//// remove entries from tmp_svd that are not also in tmp_descent
	//map<string,double>::iterator descent_end = tmp_descent.end();

	//for (map<string,double>::iterator b=tmp_svd.begin(), e=tmp_svd.end(); b!=e;) {
	//			if (tmp_descent.find((*b).first) != descent_end ) {
	//				tmp_svd.erase(b++);
	//			}
	//			else {
	//				++b;
	//			}
	//}
	for (map<string,double>::iterator b=tmp_svd.begin(), e=tmp_svd.end(); b!=e; ++b) {
		cur_numeric_pars.erase((*b).first);
		par_transform.get_frozen_ptr()->insert((*b).first, (*b).second);
	}
	return tmp_svd;
}



map<string, double> SVDASolver::limit_parameters_ip(const Parameters &init_numeric_pars, Parameters &upgrade_numeric_pars)
{
	const string *name;
	double p_init;
	double p_upgrade;
	double limit;
	map<string,double> tmp;

	pair<bool, double> par_limit;
	Parameters::iterator pu_iter;
	// TO DO this currently assumes SVD and that upgrades are orthogonal
	for(Parameters::const_iterator b=init_numeric_pars.begin(), e=init_numeric_pars.end(); b!=e; ++b)
	{
		par_limit = pair<bool, double>(false, 0.0);
		name = &(*b).first;  // parameter name
		p_init = (*b).second; // inital parameter value
		pu_iter = upgrade_numeric_pars.find(*name);
		p_upgrade = (*pu_iter).second;  // upgrade parameter value
		// only use relative change limits
		limit = abs((p_upgrade - p_init) / p_init);
		if(limit > ctl_info->relparmax)
		{
			(*pu_iter).second = p_init + sign(p_upgrade - p_init) * ctl_info->relparmax *  abs(p_init);
		}
		//cout << "limited update = " << (*pu_iter).second << endl;
	}
	return tmp;
}




SVDASolver::~SVDASolver(void)
{
}
