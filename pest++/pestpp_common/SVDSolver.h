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


class FileManager;
class ModelRun;
class QSqrtMatrix;
class PriorInformation;
class Regularization;
class SVDPackage;

class SVDSolver
{
	enum class LimitType {NONE, LBND, UBND, REL, FACT};
	enum class MarquardtMatrix {IDENT, JTQJ};
public:
	SVDSolver(const ControlInfo *_ctl_info, const SVDInfo &_svd_info, const ParameterGroupInfo *_par_group_info_ptr, const ParameterInfo *_ctl_par_info_ptr,
		const ObservationInfo *_obs_info, FileManager &_file_manager, const Observations *_observations, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, const PriorInformation *_prior_info_ptr, Jacobian &_jacobian, 
		const Regularization *regul_scheme_ptr, int _max_freeze_iter, OutputFileWriter &_output_file_writer, RestartController &_restart_controller, const string &description = "base parameter solution");
	virtual ModelRun& solve(RunManagerAbstract &run_manager, TerminationController &termination_ctl, int max_iter, ModelRun &cur_run, ModelRun &optimum_run);
	virtual void iteration(RunManagerAbstract &run_manager, TerminationController &termination_ctl, bool calc_init_obs=false);
	ModelRun &cur_model_run() {return cur_solution;}
	virtual void set_svd_package(PestppOptions::SVD_PACK _svd_pack);
	virtual ~SVDSolver(void);
protected:
	class Upgrade {
	public:
		Eigen::VectorXd uvec;
		double norm;
		vector<string> par_name_vec;
		Parameters frozen_numeric_pars;
		void print(ostream &os); 
	};

	SVDPackage *svd_package;
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
	int max_freeze_iter;
	double best_lambda;
	OutputFileWriter &output_file_writer;
	RestartController &restart_controller;

	virtual Parameters limit_parameters_ip(const Parameters &init_numeric_pars, 
		Parameters &upgrade_numeric_pars, LimitType &limit_type, 
		const Parameters &frozen_numeric_pars = Parameters());
	virtual Parameters limit_parameters_freeze_all_ip(const Parameters &init_numeric_pars, 
		Parameters &upgrade_numeric_pars, const Parameters &frozen_numeric_pars = Parameters());
	virtual const string &get_description(){return description;}
	virtual void iteration_update_and_report(ostream &os, ModelRun &upgrade, TerminationController &termination_ctl); 
	void param_change_stats(double p_old, double p_new, bool &have_fac, double &fac_change, bool &have_rel,
		double &rel_change);
	Upgrade calc_lambda_upgrade_vec(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
	const Eigen::VectorXd &Residuals, const vector<string> &par_name_vec, const vector<string> &obs_name_vec,
	const Parameters &base_numeric_pars, const Parameters &freeze_numeric_pars, int &tot_sing_val,
	double lambda, MarquardtMatrix marquardt_type=MarquardtMatrix::IDENT);
	Upgrade calc_upgrade_vec(const Upgrade &direction, 
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_numeric_pars, const Parameters &freeze_numeric_pars, double scale);
	Parameters apply_upgrade(const Parameters &init_numeric_pars,const Upgrade &upgrade, double scale = 1.0);
	void update_upgrade(Upgrade &upgrade, const Parameters &base_pars, const Parameters &new_pars, const Parameters &frozen_pars);
	void check_limits(const Parameters &init_derivative_pars, const Parameters &upgrade_derivative_pars,
		map<string, LimitType> &limit_type_map, Parameters &derivative_parameters_at_limit);
	Eigen::VectorXd calc_residual_corrections(const Jacobian &jacobian, const Parameters &del_numeric_pars, 
							   const vector<string> obs_name_vec);
};

#endif /* SVDSOLVER_H_ */