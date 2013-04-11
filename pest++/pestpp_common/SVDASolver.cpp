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
#include <iomanip>
#include <map>
#include "SVDASolver.h"
#include "ModelRunPP.h"
#include "QSqrtMatrix.h"
#include "eigen_tools.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "PriorInformation.h"

using namespace std;
using namespace pest_utils;
using namespace Eigen;

SVDASolver::SVDASolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_par_group_info_ptr, const ParameterInfo *_ctl_par_info_ptr,
		const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, const Regularization *_regul_scheme, int _n_rotation_fac)
		: SVDSolver(_ctl_info, _svd_info, _par_group_info_ptr, _ctl_par_info_ptr, _obs_info, 
		_file_manager, _observations, _obj_func, _par_transform, _prior_info_ptr, _jacobian, 
		_regul_scheme, _n_rotation_fac, "super parameter solution")
{
}

Parameters SVDASolver::freeze_parameters(Parameters &base_numeric_pars, Parameters &new_numeric_pars, Parameters &frozen_ctl_pars, bool freeze_facpar_relpar)
{
	//This routine has the side effect that it modifies par_transform transformtion
	// Test SVD parameter upgrade
	Parameters tmp_parameters;
	int ip = 0;
	Parameters tmp_svd;
	
	tmp_svd = limit_parameters_ip(base_numeric_pars, new_numeric_pars, frozen_ctl_pars);

	return tmp_svd;
}


Parameters SVDASolver::limit_parameters_ip(const Parameters &init_numeric_pars, Parameters &upgrade_numeric_pars, const Parameters &frozen_pars)
{
	const string *name;
	double val_init;
	double val_upgrade;
	double limit;
	Parameters tmp;

	pair<bool, double> par_limit;
	// TO DO this currently assumes SVD and that upgrades are orthogonal
	for(auto &ipar : init_numeric_pars)
	{
		par_limit = pair<bool, double>(false, 0.0);
		name = &(ipar.first);  // parameter name
		val_init = ipar.second; // inital parameter value
		auto iter_pu = upgrade_numeric_pars.find(*name);
		assert (iter_pu != upgrade_numeric_pars.end());
		val_upgrade = (*iter_pu).second;  // upgrade parameter value
		// only use relative change limits
		limit = abs((val_upgrade - val_init) / val_init);
		if(limit > ctl_info->relparmax)
		{
			(*iter_pu).second = val_init + sign(val_upgrade - val_init) * ctl_info->relparmax *  abs(val_init);
		}
		//cout << "limited update = " << (*pu_iter).second << endl;
	}
	return tmp;
}



SVDASolver::~SVDASolver(void)
{
}
