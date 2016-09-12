#include "sequential_lp.h"
#include "pest.h"
#include "Jacobian_1to1.h"
#include "ClpSimplex.hpp"
#include "RunManagerAbstract.h"
#include "TerminationController.h"
#include "covariance.h"
#include "FileManager.h"
#include "OutputFileWriter.h"

sequentialLP::sequentialLP(Pest &_pest_scenario, RunManagerAbstract* _run_mgr_ptr,
	TerminationController* _termination_ctl_ptr, Covariance &_parcov, FileManager &_file_mgr, 
	OutputFileWriter* _out_wtr_ptr) : pest_scenario(_pest_scenario), run_mgr_ptr(_run_mgr_ptr),
	termination_ctl_ptr(_termination_ctl_ptr), parcov(_parcov), file_mgr(_file_mgr), out_wtr_ptr(_out_wtr_ptr)
{
	initialize();
}

void sequentialLP::initialize()
{
	//TODO: handle restart condition
	//TODO: handle base jco condition
	separate_scenarios();
	decision_vars = opt_scenario.get_ctl_parameters();
	constraints_obs = opt_scenario.get_ctl_observations();
	constraints_sim = Observations(constraints_obs);
	obj_func = ObjectiveFunc(&constraints_obs, &opt_scenario.get_ctl_observation_info(),
		       null_prior);
	optimum_run = ModelRun(&obj_func,constraints_sim);
	optimum_run.set_ctl_parameters(decision_vars);
	//TODO: error checking:
	//noptmax > 0
	
}

void sequentialLP::solve_lp_problem()
{

}

void sequentialLP::solve()
{
	//TODO: handle restart condition
	//TODO: handle base jco condition
	Jacobian_1to1 jco(file_mgr);
	jco.set_base_numeric_pars(decision_vars);
	jco.set_base_sim_obs(constraints_sim);

	for (int i = 0; i < pest_scenario.get_control_info().noptmax; i++)
	{
		make_runs(jco);
		solve_lp_problem();
	}
}

void sequentialLP::make_runs(Jacobian_1to1 &jco)
{
	ParamTransformSeq par_trans = opt_scenario.get_base_par_tran_seq();
	set<string> out_of_bounds;
	jco.build_runs(optimum_run, opt_scenario.get_ctl_ordered_par_names(), par_trans,
		opt_scenario.get_base_group_info(), opt_scenario.get_ctl_parameter_info(),
		*run_mgr_ptr, out_of_bounds);
	jco.make_runs(*run_mgr_ptr);
}

void sequentialLP::separate_scenarios()
{
	//if needed separate scenarios, otherwise, just set the opt_scenario to the pest_scenario.

	//this might be using the copy constructor, needs to be fixed
	opt_scenario = pest_scenario;
}

void sequentialLP::update()
{

}

void sequentialLP::update_decision_vars()
{

}

void sequentialLP::update_constraints()
{

}



