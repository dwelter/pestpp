#include "sequential_lp.h"
#include "pest.h"
#include "Jacobian.h"
#include "ClpSimplex.hpp"
#include "RunManagerAbstract.h"
#include "TerminationController.h"
#include "covariance.h"
#include "FileManager.h"
#include "OutputFileWriter.h"

sequentialLP::sequentialLP(Pest &_pest_scenario, Jacobian &_jco, RunManagerAbstract* _run_mgr_ptr,
	TerminationController* _termination_ctl_ptr, Covariance &_parcov, FileManager &_file_mgr, 
	OutputFileWriter* _out_wtr_ptr) : pest_scenario(_pest_scenario), jco(_jco), run_mgr_ptr(_run_mgr_ptr),
	termination_ctl_ptr(_termination_ctl_ptr), parcov(_parcov), file_mgr(_file_mgr), out_wtr_ptr(_out_wtr_ptr)
{
	initialize();
}

void sequentialLP::initialize()
{
	separate_scenarios();

}

void sequentialLP::solve_lp_problem()
{

}

void sequentialLP::solve()
{

}

void sequentialLP::make_runs()
{

}

void sequentialLP::separate_scenarios()
{

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



