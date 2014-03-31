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
#ifndef SVDSOLVER_H_
#define SVDSOLVER_H_

#include <map>
#include <iomanip>
#include <Eigen/Dense>
#include "Transformable.h"
#include "ParamTransformSeq.h"
#include "Jacobian.h"
#include "pest_data_structs.h"
#include "ModelRunPP.h"
#include "TerminationController.h"
#include "RunManagerAbstract.h"
#include "OutputFileWriter.h"
#include "RestartController.h"
#include "PerformanceLog.h"


class FileManager;
class ModelRun;
class QSqrtMatrix;
class PriorInformation;
class Regularization;
class SVDPackage;

class SVDSolver
{
public:
	enum class MAT_INV{ Q12J, JTQJ };
protected:
	enum class LimitType {NONE, LBND, UBND, REL, FACT};
	enum class MarquardtMatrix {IDENT, JTQJ};
public:
	SVDSolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_par_group_info_ptr, const ParameterInfo *_ctl_par_info_ptr,
		const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, 
		const Regularization *regul_scheme_ptr, OutputFileWriter &_output_file_writer, RestartController &_restart_controller, 
		SVDSolver::MAT_INV _mat_inv, PerformanceLog *_performance_log, const string &description = string("base parameter solution"));
	virtual ModelRun& solve(RunManagerAbstract &run_manager, TerminationController &termination_ctl, int max_iter, ModelRun &cur_run, ModelRun &optimum_run);
	virtual void iteration(RunManagerAbstract &run_manager, TerminationController &termination_ctl, bool calc_init_obs=false);
	ModelRun &cur_model_run() {return cur_solution;}
	virtual void set_svd_package(PestppOptions::SVD_PACK _svd_pack);
	virtual ParameterGroupInfo get_parameter_group_info() const { return *par_group_info_ptr; }
	virtual ~SVDSolver(void);
protected:
	class Upgrade {
	public:
		Eigen::VectorXd uvec;
		double norm;
		vector<string> par_name_vec;
		Parameters frozen_numeric_pars; 
	};

	SVDPackage *svd_package;
	MAT_INV mat_inv;
	const string description;
	const ControlInfo *ctl_info;
	SVDInfo svd_info;
	ObjectiveFunc *obj_func;
	ModelRun cur_solution;
	const ParameterInfo *ctl_par_info_ptr;
	const ParameterGroupInfo *par_group_info_ptr;
	ParamTransformSeq par_transform;
	const Observations *observations_ptr;
	const ObservationInfo *obs_info_ptr;
	const PriorInformation *prior_info_ptr;
	const Regularization *regul_scheme_ptr;
	FileManager &file_manager;
	Jacobian &jacobian;
	bool phiredswh_flag;
	bool save_next_jacobian;
	double prev_phi_percent;
	int num_no_descent;
	double best_lambda;
	OutputFileWriter &output_file_writer;
	RestartController &restart_controller;
	PerformanceLog *performance_log;

	//virtual void limit_parameters_ip(const Parameters &init_numeric_pars, 
	//	Parameters &upgrade_numeric_pars, LimitType &limit_type, 
	//	const Parameters &frozen_numeric_pars = Parameters());
	virtual Parameters limit_parameters_freeze_all_ip(const Parameters &init_ctl_pars, 
		Parameters &upgrade_ctl_pars, const Parameters &frozen_ctl_pars = Parameters());
	virtual const string &get_description(){return description;}
	virtual void iteration_update_and_report(ostream &os, ModelRun &upgrade, TerminationController &termination_ctl); 
	void param_change_stats(double p_old, double p_new, bool &have_fac, double &fac_change, bool &have_rel,
		double &rel_change);
	void calc_upgrade_vec(double i_lambda, Parameters &frozen_ctl_pars, QSqrtMatrix &Q_sqrt, Eigen::VectorXd &residuals_vec,
		vector<string> &obs_names_vec, const Parameters &base_run_ctl_pars, LimitType &limit_type,
		Parameters &new_ctl_pars, MarquardtMatrix marquardt_type);
	void calc_lambda_upgrade_vecQ12J(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_ctl_pars, const Parameters &freeze_ctl_pars,
		double lambda, Parameters &ctl_upgrade_pars, Parameters &upgrade_ctl_del_pars,
		Parameters &grad_ctl_del_pars, MarquardtMatrix marquardt_type);
	void calc_lambda_upgrade_vec_JtQJ(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_ctl_pars, const Parameters &freeze_ctl_pars,
		double lambda, Parameters &ctl_upgrade_pars, Parameters &upgrade_ctl_del_pars,
		Parameters &grad_ctl_del_pars, MarquardtMatrix marquardt_type);
	void check_limits(const Parameters &init_ctl_pars, const Parameters &upgrade_ctl_pars,
		map<string, LimitType> &limit_type_map, Parameters &ctl_parameters_at_limit);
	Eigen::VectorXd calc_residual_corrections(const Jacobian &jacobian, const Parameters &del_numeric_pars, 
							   const vector<string> obs_name_vec);
	bool par_heading_out_bnd(double org_par, double new_par, double lower_bnd, double upper_bnd);
	int check_bnd_par(Parameters &new_freeze_ctl_pars, const Parameters &current_ctl_pars, const Parameters &new_upgrade_ctl_pars, const Parameters &new_grad_ctl_pars=Parameters());
};

#endif /* SVDSOLVER_H_ */
