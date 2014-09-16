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
#ifndef SVDASOLVER_H_
#define SVDASOLVER_H_

#include <map>
#include <Eigen/Dense>
#include "Transformable.h"
#include "ParamTransformSeq.h"
#include "Jacobian.h"
#include "pest_data_structs.h"
#include "ModelRunPP.h"
#include "TerminationController.h"
#include "SVDSolver.h"
#include "PerformanceLog.h"


class SVDASolver : public SVDSolver
{
public:
	SVDASolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_base_parameter_group_info_ptr, const ParameterInfo *_ctl_par_info_ptr,
		const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, DynamicRegularization *_regul_scheme,
		OutputFileWriter &_output_file_writer, RestartController &_restart_controller, SVDSolver::MAT_INV _mat_inv, PerformanceLog *_performance_log, 
		const std::vector<double> &_base_lambda_vec, bool _phiredswh_flag, bool _splitswh_flag, int max_super_frz_iter);
	virtual Parameters limit_parameters_freeze_all_ip(const Parameters &init_active_ctl_pars,
		Parameters &upgrade_active_ctl_pars, const Parameters &prev_frozen_active_ctl_pars);
	virtual void calc_upgrade_vec(double i_lambda, Parameters &prev_frozen_active_ctl_pars, QSqrtMatrix &Q_sqrt, 
        const DynamicRegularization &regul, Eigen::VectorXd &residuals_vec,
		std::vector<std::string> &obs_names_vec, const Parameters &base_run_active_ctl_pars, Parameters &new_active_ctl_pars,
		MarquardtMatrix marquardt_type, bool scale_upgrade=false);
	void iteration_jac(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool calc_init_obs = false);
	ModelRun iteration_upgrd(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run);
	virtual const string &get_description(){return description;}
	virtual ParameterGroupInfo get_parameter_group_info() const { return super_parameter_group_info; }
	bool get_calc_jacobian() const { return calc_jacobian;}
	void set_calc_jacobian(bool tmp) { calc_jacobian=tmp; }
	virtual string get_solver_type() const { return svda_solver_type_name; }
	~SVDASolver(void);
private:
	const static string svda_solver_type_name;
	ParameterGroupInfo super_parameter_group_info;
	bool calc_jacobian;
	int max_super_frz_iter;
};


#endif //SVDASOLVER_H_
