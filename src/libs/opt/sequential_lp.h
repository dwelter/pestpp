#pragma once
#ifndef SEQUENTIAL_LP_H
#define SEQUENTIAL_LP_H

#include "pest.h"
#include "Jacobian.h"
#include "ClpSimplex.hpp"
#include "RunManagerAbstract.h"
#include "TerminationController.h"
#include "covariance.h"
#include "FileManager.h"
#include "OutputFileWriter.h"

class sequentialLP
{
public:
	sequentialLP(Pest &_pest_scenario, Jacobian &_jco, RunManagerAbstract* _run_mgr_ptr, 
		         TerminationController* _termination_ctl_ptr, Covariance &_parcov, 
				 FileManager &_file_mgr, OutputFileWriter* _out_wtr_ptr);
	void initialize();
	void solve();
	
private:
	Pest pest_scenario;
	Pest opt_scenario;
	Jacobian jco;
	RunManagerAbstract* run_mgr_ptr;
	TerminationController* termination_ctl_ptr;
	Covariance parcov;
	FileManager file_mgr;
	OutputFileWriter* out_wtr_ptr;
	void solve_lp_problem();
	void update();
	void update_decision_vars();
	void update_constraints();
	void separate_scenarios();
	void make_runs();
	void process_runs();




};




#endif