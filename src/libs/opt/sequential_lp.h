#pragma once
#ifndef SEQUENTIAL_LP_H
#define SEQUENTIAL_LP_H

#include "pest.h"
#include "Jacobian_1to1.h"
#include "ClpSimplex.hpp"
#include "RunManagerAbstract.h"
#include "TerminationController.h"
#include "covariance.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "Transformable.h"
#include "ModelRunPP.h"

class sequentialLP
{
	enum ConstraintSense {less_than,greater_than,equal_to};
public:
	sequentialLP(Pest &_pest_scenario, RunManagerAbstract* _run_mgr_ptr, 
		         Covariance &_parcov, 
				 FileManager &_file_mgr, OutputFileWriter* _out_wtr_ptr);
	void initialize_and_check();
	void solve();
	
	//ModelRun get_optimum_run() { return optimum_run; }

private:
	string obj_func_str;
	map<string, double> obj_func_coef_map;
	Observations obj_func_obs;
	ObservationInfo obj_func_info;
	double* ctl_ord_obj_func_coefs;
	int slp_iter;
	vector<double> iter_obj_values;
	map<string, ConstraintSense> constraint_sense_map;
	map <string, string> constraint_sense_name;
	vector<string> ctl_ord_dec_var_names;
	vector<string> ctl_ord_constraint_names;
	double* dec_var_lb;
	double* dec_var_ub;
	double* constraint_lb;
	double* constraint_ub;
	PriorInformation* null_prior = new PriorInformation();
	//ModelRun current_run;
	//ModelRun optimum_run;
	//ObjectiveFunc obj_func;
	Parameters all_pars_and_dec_vars;
	ParamTransformSeq par_trans;
	Observations constraints_obs;
	Observations constraints_sim;
	Pest pest_scenario;
	//Pest opt_scenario;
	RunManagerAbstract* run_mgr_ptr;
	//TerminationController* termination_ctl_ptr;
	Covariance parcov;
	FileManager file_mgr;
	OutputFileWriter* out_wtr_ptr;
	ClpSimplex solve_lp_problem(Jacobian_1to1 &jco);

	void initialize_obj_function();
	void initialize_dec_vars();
	void initialize_constraints();
	void initial_report();
	void presolve_constraint_report();
	pair<double,double> postsolve_decision_var_report(Parameters &upgrade_pars);
	void postsolve_constraint_report(Observations &upgrade_obs);
	void update(ClpSimplex &model);
	void update_and_report_decision_vars(ClpSimplex &model);
	void update_and_report_constraints(ClpSimplex &model);
	void separate_scenarios();
	void make_response_matrix_runs(Jacobian_1to1 &jco);
	bool make_upgrade_run(Parameters &upgrade_pars, Observations &upgrade_obs);
	void process_model(ClpSimplex &model);
	CoinPackedMatrix jacobian_to_coinpackedmatrix(Jacobian_1to1 &jco);
	void build_constraint_bound_arrays();
	void throw_sequentialLP_error(string message);
	void throw_sequentialLP_error(string message,const vector<string> &messages);
	void throw_sequentialLP_error(string message, const set<string> &messages);
	vector<double> get_constraint_residual_vec();
	vector<double> get_constraint_residual_vec(Observations &sim_vals);
	void build_obj_func_coef_array();

};




#endif