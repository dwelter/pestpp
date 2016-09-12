#include "sequential_lp.h"
#include "pest.h"
#include "Jacobian_1to1.h"
#include "ClpSimplex.hpp"
#include "RunManagerAbstract.h"
#include "TerminationController.h"
#include "covariance.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include <Eigen/Sparse>

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

	//set decision vars attrib and ordered dec var name vec
	decision_vars = opt_scenario.get_ctl_parameters();
	for (auto &name : pest_scenario.get_ctl_ordered_par_names())
	{
		if (decision_vars.find(name) != decision_vars.end())
		{
			ctl_ord_dec_var_names.push_back(name);
		}
	}

	
	//set the two constraints attribs and ordered constraint name vec
	constraints_obs = opt_scenario.get_ctl_observations();
	constraints_sim = Observations(constraints_obs);
	for (auto &name : pest_scenario.get_ctl_ordered_obs_names())
	{
		if (constraints_obs.find(name) != constraints_obs.end())
		{
			ctl_ord_constraint_names.push_back(name);
		}
	}
	
	obj_func = ObjectiveFunc(&constraints_obs, &opt_scenario.get_ctl_observation_info(),
		       null_prior);
	optimum_run = ModelRun(&obj_func,constraints_sim);
	optimum_run.set_ctl_parameters(decision_vars);
	//TODO: error checking:
	//noptmax > 0
	
}

ClpSimplex sequentialLP::solve_lp_problem(Jacobian_1to1 &jco)
{
	//convert Jacobian_1to1 to CoinPackedMatrix
	CoinPackedMatrix matrix = jacobian_to_coinpackedmatrix(jco);

	ClpSimplex model;

	return model;

}

CoinPackedMatrix sequentialLP::jacobian_to_coinpackedmatrix(Jacobian_1to1 &jco)
{
	
	//for now, for debugging, just get an ordered jco. later, this can be more flexible and performant
	//Eigen::SparseMatrix<double> eig_ord_jco = jco.get_matrix(ctl_ord_dec_var_names,ctl_ord_constraint_names);
	
	Eigen::SparseMatrix<double> eig_ord_jco = jco.get_matrix();

	cout << "number of nonzero elements in jacobian: " << eig_ord_jco.nonZeros() << endl;


	//the triplet elements to pass to the coinpackedmatrix constructor
	int * row_idx = new int[eig_ord_jco.nonZeros()];
	int * col_idx = new int[eig_ord_jco.nonZeros()];
	double * elems = new double[eig_ord_jco.nonZeros()];

	int elem_count = 0;



	//iterate through the eigen sparse matrix
	for (int i = 0; i < eig_ord_jco.outerSize(); ++i)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(eig_ord_jco, i); it; ++it)
		{
			row_idx[elem_count] = it.row();
			col_idx[elem_count] = it.col();
			elems[elem_count] = it.value();
			elem_count++;
		}
	}
	if (elem_count != eig_ord_jco.nonZeros())
	{
		throw runtime_error("sequentialLP::jacobian_to_coinpackedmatrix() error: wrong number of triplet components");
	}

	if (elem_count == 0)
	{
		throw runtime_error("sequentialLP::jacobian_to_coinpackedmatrix() error: zero triplets found");
	}

	CoinPackedMatrix matrix(true,row_idx,col_idx,elems,elem_count);
	return matrix;
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
		ClpSimplex model = solve_lp_problem(jco);
		update(model);
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
	jco.process_runs(par_trans, opt_scenario.get_base_group_info(), *run_mgr_ptr, *null_prior, false);
}

void sequentialLP::separate_scenarios()
{
	//if needed separate scenarios, otherwise, just set the opt_scenario to the pest_scenario.

	//this might be using the copy constructor, needs to be fixed
	opt_scenario = pest_scenario;
}

void sequentialLP::update(ClpSimplex &model)
{
	update_decision_vars(model);
	update_constraints(model);
	return;
}

void sequentialLP::update_decision_vars(ClpSimplex &model)
{
	return;
}

void sequentialLP::update_constraints(ClpSimplex &model)
{
	return;
}



