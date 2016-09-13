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
#include "CoinFinite.hpp"
#include "ClpPresolve.hpp"
#include <iomanip>

sequentialLP::sequentialLP(Pest &_pest_scenario, RunManagerAbstract* _run_mgr_ptr,
	TerminationController* _termination_ctl_ptr, Covariance &_parcov, FileManager &_file_mgr, 
	OutputFileWriter* _out_wtr_ptr) : pest_scenario(_pest_scenario), run_mgr_ptr(_run_mgr_ptr),
	termination_ctl_ptr(_termination_ctl_ptr), parcov(_parcov), file_mgr(_file_mgr), out_wtr_ptr(_out_wtr_ptr)
{
	initialize_and_check();
}

void sequentialLP::throw_squentialLP_error(string message)
{
	string error_message = "error in sequentialLP process: " + message;
	file_mgr.rec_ofstream() << error_message << endl;
	file_mgr.close_file("rec");
	throw runtime_error(error_message);


}

void sequentialLP::initial_report()
{
	ofstream* f_rec = &file_mgr.rec_ofstream();
	*f_rec << endl << "  -------------------------------------------------------------" << endl;
	*f_rec << "  ---  sequential linear programming problem information  ---  " << endl;
	*f_rec << "  -------------------------------------------------------------" << endl << endl;
	
	*f_rec << "  ---  decision variables active in LP  ---  " << endl;
	for (auto &name : ctl_ord_dec_var_names)
		*f_rec << setw(20) << left << name << endl;
	*f_rec << " note: bound and initial value info reported in 'parameter' section" << endl << endl;

	*f_rec << "  ---  constraints in LP  ---  " << endl;
	constraint_report("rhs",constraints_obs);
	return;
}

void sequentialLP::constraint_report(string fieldname, Observations &constraints)
{
	ofstream* f_rec = &file_mgr.rec_ofstream();
	*f_rec << setw(20) << "name" << setw(20) << "sense" << setw(20) << fieldname << endl;
	for (auto &name : ctl_ord_constraint_names)
	{
		*f_rec << setw(20) << left << name;
		*f_rec << setw(20) << constraint_sense_name[int(constraint_sense_map[name])];
		*f_rec << setw(20) << constraints.get_rec(name) << endl;
	}

	
	return;

}

void sequentialLP::initialize_and_check()
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
	
	//TODO: override base objective func class with linear obj func class
	obj_func = ObjectiveFunc(&constraints_obs, &opt_scenario.get_ctl_observation_info(),
		       null_prior);
	

	//initialize the current and optimum model run instances
	optimum_run = ModelRun(&obj_func,constraints_sim);
	optimum_run.set_ctl_parameters(decision_vars);
	current_run = ModelRun(&obj_func, constraints_sim);
	current_run.set_ctl_parameters(decision_vars);


	//set the decision var lower and upper bound arrays
	dec_var_lb = new double[ctl_ord_dec_var_names.size()];
	dec_var_ub = new double[ctl_ord_constraint_names.size()];
	Parameters parlbnd = opt_scenario.get_ctl_parameter_info().get_low_bnd(ctl_ord_dec_var_names);
	Parameters parubnd = opt_scenario.get_ctl_parameter_info().get_up_bnd(ctl_ord_dec_var_names);
	for (int i = 0; i < ctl_ord_dec_var_names.size(); ++i)
	{
		dec_var_lb[i] = parlbnd.get_rec(ctl_ord_dec_var_names[i]);
		dec_var_ub[i] = parubnd.get_rec(ctl_ord_dec_var_names[i]);
	}

	//build map of constraint sense
	map<string, string> problem_constraints;
	for (auto &name : ctl_ord_constraint_names)
	{
		string group = opt_scenario.get_ctl_observation_info().get_observation_rec_ptr(name)->group;
		if (group == "L")
			constraint_sense_map[name] = ConstraintSense::less_than;
		else if (group == "G")
			constraint_sense_map[name] = ConstraintSense::greater_than;
		else if ((group == "E") || (group == "N"))
			constraint_sense_map[name] = ConstraintSense::equal_to;
		else
			problem_constraints[name] = group;
	}
	if (problem_constraints.size() > 0)
	{
		stringstream ss;
		ss << endl;
		for (auto &pc : problem_constraints)
		{
			ss << "name:" << pc.first << ", group:" << pc.second << endl;
		}
		throw_squentialLP_error("the following constraints do not have the correct group names {'l','g','e'}: " + ss.str());
	}
	
	//allocate the constraint bound arrays 
	constraint_lb = new double[ctl_ord_constraint_names.size()];
	constraint_ub = new double[ctl_ord_constraint_names.size()];

	//TODO: error checking:
	//noptmax > 0
	//no log transform for decision vars
	return;
}

void sequentialLP::build_constraint_bound_arrays()
{
	
	vector<double> residuals = current_run.get_residuals_vec(ctl_ord_constraint_names);

	for (int i = 0; i < ctl_ord_constraint_names.size(); ++i)
	{
		string name = ctl_ord_constraint_names[i];
		if (constraint_sense_map[name] == ConstraintSense::less_than)
			constraint_ub[i] = residuals[i];
		else
			constraint_ub[i] = COIN_DBL_MAX;

		if (constraint_sense_map[name] == ConstraintSense::greater_than)
			constraint_lb[i] = residuals[i];
		else
			constraint_lb[i] = -COIN_DBL_MAX;
		if (constraint_sense_map[name] == ConstraintSense::equal_to)
		{
			constraint_ub[i] = residuals[i];
			constraint_lb[i] = residuals[i];
		}
	}

	return;
}

ClpSimplex sequentialLP::solve_lp_problem(Jacobian_1to1 &jco)
{
	//convert Jacobian_1to1 to CoinPackedMatrix
	CoinPackedMatrix matrix = jacobian_to_coinpackedmatrix(jco);

	//set/update the constraint bound arrays
	build_constraint_bound_arrays();

	//temp obj function
	double* objective_func = new double[ctl_ord_dec_var_names.size()];
	for (int i = 0; i < ctl_ord_dec_var_names.size(); ++i)
	{
		objective_func[i] = 1.0;
	}
	
	//instantiate and load the linear simplex model
	ClpSimplex model;
	model.loadProblem(matrix, dec_var_lb, dec_var_ub, objective_func, constraint_lb, constraint_ub);

	//solve the linear program
	ClpPresolve presolve_info;
	ClpSimplex* presolved_model = presolve_info.presolvedModel(model);
	//if presolvedModel is Null, then it is primal infeasible, so 
	//try the dual
	if (!presolved_model)
	{
		presolved_model->dual();
	}
	//update the status arrays of both the presolve and original models
	presolve_info.postsolve(true);

	//check the solution
	model.checkSolution();

	//because of numerical tolerances, solve one more time
	model.primal(1);

	//return the model for reporting purposes
	return model;
}

CoinPackedMatrix sequentialLP::jacobian_to_coinpackedmatrix(Jacobian_1to1 &jco)
{
	
	Eigen::SparseMatrix<double> eig_ord_jco = jco.get_matrix(ctl_ord_constraint_names, ctl_ord_dec_var_names);

	if (eig_ord_jco.nonZeros() != jco.get_matrix_ptr()->nonZeros())
		throw_squentialLP_error("sequentialLP::jacobian_to_coinpackedmatrix() error: ordered jco nnz != org jco nnz");
	
	cout << "number of nonzero elements in response matrix: " << eig_ord_jco.nonZeros() << endl;


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
		throw_squentialLP_error("sequentialLP::jacobian_to_coinpackedmatrix() error: wrong number of triplet components");
	}

	if (elem_count == 0)
	{
		throw_squentialLP_error("sequentialLP::jacobian_to_coinpackedmatrix() error: zero triplets found");
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

	//get the base run and update current_run
	Parameters temp_pars;
	Observations temp_obs;
	run_mgr_ptr->get_run(0, temp_pars, temp_obs,false);
	par_trans.model2ctl_ip(temp_pars);
	current_run.update_ctl(temp_pars, temp_obs);

	//process the remaining responses
	jco.process_runs(par_trans, opt_scenario.get_base_group_info(), *run_mgr_ptr, *null_prior, false);
	return;
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



