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

vector<double> sequentialLP::get_constraint_residual_vec()
{
	vector<double> residuals_vec;
	residuals_vec.resize(ctl_ord_constraint_names.size(), 0.0);

	Observations::const_iterator found_obs;
	Observations::const_iterator not_found_obs = constraints_sim.end();
	PriorInformation::const_iterator found_prior_info;
	
	int i = 0;
	for (vector<string>::iterator b = ctl_ord_constraint_names.begin(), e = ctl_ord_constraint_names.end(); b != e; ++b, ++i)
	{
		found_obs = constraints_sim.find(*b);
		if (found_obs != not_found_obs)
		{
			residuals_vec[i] = constraints_obs.get_rec(*b) - (*found_obs).second;
		}
		
	}
	return residuals_vec;
}

void sequentialLP::initial_report()
{
	ofstream* f_rec = &file_mgr.rec_ofstream();
	*f_rec << endl << "  -------------------------------------------------------------" << endl;
	*f_rec << "  ---  sequential linear programming problem information  ---  " << endl;
	*f_rec << "  -------------------------------------------------------------" << endl << endl;
	
	*f_rec << "-->number of iterations of sequential linear programming (noptmax): " << opt_scenario.get_control_info().noptmax << endl;
	
	*f_rec << "  ---  decision variables active in SLP  ---  " << endl;
	for (auto &name : ctl_ord_dec_var_names)
		*f_rec << setw(20) << left << name << endl;
	*f_rec << " note: bound and initial value info reported in 'parameter' section" << endl << endl;

	*f_rec << "  ---  constraints in SLP  ---  " << endl;
	*f_rec << setw(20) << "name" << setw(20) << "sense" << setw(20) << "value" << endl;
	for (auto &name : ctl_ord_constraint_names)
	{
		*f_rec << setw(20) << left << name;
		*f_rec << setw(20) << constraint_sense_name[name];
		*f_rec << setw(20) << constraints_obs.get_rec(name) << endl;
	}
	return;
}

void sequentialLP::constraint_report()
{
	ofstream* f_rec = &file_mgr.rec_ofstream();
	vector<double> residuals = get_constraint_residual_vec();
	*f_rec << endl << "  constraint information at start of iteration " << slp_iter << endl;
	*f_rec << setw(20) << "name" << setw(10) << "sense" << setw(15) << "value";
	*f_rec << setw(15) << "residual" << setw(15) << "lower bound" << setw(15) << "upper bound" << endl;

	for (int i=0;i<ctl_ord_constraint_names.size();++i)
	{
		string name = ctl_ord_constraint_names[i];
		*f_rec << setw(20) << left << name;
		*f_rec << setw(10) << constraint_sense_name[name];
		*f_rec << setw(15) << constraints_obs.get_rec(name);
		*f_rec << setw(15) << residuals[i];
		*f_rec << setw(15) << constraint_lb[i];
		*f_rec << setw(15) << constraint_ub[i] << endl;

	}
	return;
}

void sequentialLP::build_obj_function_components()
{
	for (auto &name : ctl_ord_dec_var_names)
	{
		ObservationRec obs_rec;
		obj_func_obs.insert(name, decision_vars.get_rec(name));
		obs_rec.group = "SLP phi";
		obs_rec.weight = obj_func_coef_map[name];
		obj_func_info.observations[name] = obs_rec;
		
	}
	return;
}


void sequentialLP::initialize_and_check()
{
	ofstream *f_rec = &file_mgr.rec_ofstream();
	//TODO: handle restart condition
	//TODO: handle base jco condition
	separate_scenarios();

	if (opt_scenario.get_control_info().noptmax < 1)
		throw_squentialLP_error("noptmax must be greater than 0");

	//set decision vars attrib and ordered dec var name vec
	//and check for illegal parameter transformations
	decision_vars = opt_scenario.get_ctl_parameters();
	vector<string> problem_trans;
	for (auto &name : pest_scenario.get_ctl_ordered_par_names())
	{
		if (decision_vars.find(name) != decision_vars.end())
		{
			if (opt_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(name)->tranform_type != ParameterRec::TRAN_TYPE::NONE)
				problem_trans.push_back(name);
			ctl_ord_dec_var_names.push_back(name);
		}
	}

	//if any decision vars have a transformation that is not allowed
	if (problem_trans.size() > 0)
	{
		stringstream ss;
		for (auto &name : problem_trans)
			ss << ',' << name;
		throw_squentialLP_error("the following decision variables don't have 'none' type parameter transformation: " + ss.str());
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
	
	//initialize the objective function
	obj_func_str = opt_scenario.get_pestpp_options().get_opt_obj_func();
	if (empty(obj_func_str))
	{
		*f_rec << " warning: no ++opt_objective_function-->forming a generic objective function (1.0 coef for each decision var" << endl;
		for (auto &name : ctl_ord_dec_var_names)
			obj_func_coef_map[name] = 1.0;
	}
	else
	{
		//check if the obj_str is an observation
		if (pest_scenario.get_ctl_observations().find(obj_func_str) != pest_scenario.get_ctl_observations().end())
		{
			throw_squentialLP_error("observation-based objective function not implemented");
		}
		else if (pest_scenario.get_prior_info().find(obj_func_str) != pest_scenario.get_prior_info().end())
		{
			obj_func_coef_map = pest_scenario.get_prior_info().get_pi_rec_ptr(obj_func_str).get_atom_factors();
		}
		else
			throw_squentialLP_error("unrecognized ++opt_objective_function arg: " + obj_func_str);
	}

	//check that all obj_coefs are decsision vars
	vector<string> missing_vars;
	for (auto &coef : obj_func_coef_map)
		if (find(ctl_ord_dec_var_names.begin(), ctl_ord_dec_var_names.end(), coef.first) == ctl_ord_dec_var_names.end())
			missing_vars.push_back(coef.first);
	if (missing_vars.size() > 0)
	{
		stringstream ss;
		for (auto &name : missing_vars)
			ss << ', ' << name;
		throw_squentialLP_error("the following objective function components are not decision variables: " + ss.str());
	}


	//this is nasty...setup objective function components as obseravtions and obs info
	build_obj_function_components();
	///obj_func = ObjectiveFunc(&constraints_obs, &opt_scenario.get_ctl_observation_info(),
	//	       null_prior);
	
	//initialize the current and optimum model run instances
	//optimum_run = ModelRun(&obj_func,constraints_sim);
	//optimum_run.set_ctl_parameters(decision_vars);
	//current_run = ModelRun(&obj_func, constraints_sim);
	//current_run.set_ctl_parameters(decision_vars);


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
		{
			constraint_sense_map[name] = ConstraintSense::less_than;
			constraint_sense_name[name] = "less than";
		}
		else if (group == "G")
		{
			constraint_sense_map[name] = ConstraintSense::greater_than;
			constraint_sense_name[name] = "greater than";
		}
		else if ((group == "E") || (group == "N"))
		{
			constraint_sense_map[name] = ConstraintSense::equal_to;
			constraint_sense_name[name] = "equal to";
		}

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
	
	initial_report();
	return;
}

void sequentialLP::build_constraint_bound_arrays()
{
	
	vector<double> residuals = get_constraint_residual_vec();

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
	
	ofstream* f_rec = &file_mgr.rec_ofstream();

	//convert Jacobian_1to1 to CoinPackedMatrix
	CoinPackedMatrix matrix = jacobian_to_coinpackedmatrix(jco);

	//set/update the constraint bound arrays
	build_constraint_bound_arrays();

	constraint_report();

	//temp obj function
	double* objective_func = new double[ctl_ord_dec_var_names.size()];
	for (int i = 0; i < ctl_ord_dec_var_names.size(); ++i)
	{
		objective_func[i] = 1.0;
	}
	
	//instantiate and load the linear simplex model
	ClpSimplex model;
	model.loadProblem(matrix, dec_var_lb, dec_var_ub, objective_func, constraint_lb, constraint_ub);
	model.setLogLevel(opt_scenario.get_pestpp_options().get_opt_coin_loglev());

	*f_rec << "  ---  solving linear program for iteration " << slp_iter << "  ---  " << endl;
	cout << "  ---  solving linear program for iteration " << slp_iter << "  ---  " << endl;
	
	//solve the linear program
	ClpPresolve presolve_info;
	ClpSimplex* presolved_model = presolve_info.presolvedModel(model);
	
	//if presolvedModel is Null, then it is primal infeasible, so 
	//try the dual
	if (!presolved_model)
	{
		*f_rec << "  ---  primal presolve model infeasible, trying dual..." << endl;
		cout << "  ---  primal presolve model infeasible, trying dual..." << endl;
		presolved_model->dual();
	}
	
	//update the status arrays of both the presolve and original models
	presolve_info.postsolve(true);

	//check the solution
	model.checkSolution();

	if (model.isProvenOptimal())
	{
		*f_rec << " iteration " << slp_iter << " linear solution is proven optimal" << endl << endl;
		cout << " iteration " << slp_iter << " linear solution is proven optimal" << endl << endl;
	}
	else if ((model.isProvenPrimalInfeasible()) && (model.isProvenDualInfeasible()))
		throw_squentialLP_error("both primal and dual solutions are proven infeasible...cannot continue");
	//otherwise, try again...
	else
	{
		if (model.primalFeasible())
		{
			*f_rec << endl << "-->resolving primal in attempt to prove optimal..." << endl;
			model.primal(1);
		}
		else
		{
			*f_rec << endl << "-->resolving dual in attempt to prove optimal..." << endl;
			model.dual(1);
		}
		//check the solution
		model.checkSolution();
		if (model.isProvenOptimal())
		{
			*f_rec << "  ---  iteration " << slp_iter << " linear solution is proven optimal  ---  " << endl << endl;
			cout << "  ---  iteration " << slp_iter << " linear solution is proven optimal  ---  " << endl << endl;
		}
		else
			*f_rec << "  warning: iteration " << slp_iter << " linear solution is not optimal...continuing" << endl << endl;
	}


	*f_rec << endl << endl << "  ---  linear program solution complete for iteration " << slp_iter << "  ---  " << endl;
	cout << endl << endl << "  ---  linear program solution complete for iteration " << slp_iter << "  ---  " << endl;

	//return the model for updating and reporting purposes
	return model;
}

CoinPackedMatrix sequentialLP::jacobian_to_coinpackedmatrix(Jacobian_1to1 &jco)
{
	
	Eigen::SparseMatrix<double> eig_ord_jco = jco.get_matrix(ctl_ord_constraint_names, ctl_ord_dec_var_names);

	if (eig_ord_jco.nonZeros() != jco.get_matrix_ptr()->nonZeros())
		throw_squentialLP_error("sequentialLP::jacobian_to_coinpackedmatrix() error: ordered jco nnz != org jco nnz");
	
	file_mgr.rec_ofstream() << "number of nonzero elements in response matrix: " << eig_ord_jco.nonZeros() << " of " << eig_ord_jco.size() << endl;
	cout << "number of nonzero elements in response matrix: " << eig_ord_jco.nonZeros() << " of " << eig_ord_jco.size() << endl;


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
	ofstream* f_rec = &file_mgr.rec_ofstream();
	slp_iter = 1;
	while (true)
	{
		*f_rec << endl << endl << "  ---------------------------------" << endl;
		*f_rec <<         "  --- starting LP iteration " << slp_iter << "  ---  " << endl;
		*f_rec << "  ---------------------------------" << endl << endl << endl;
		cout << endl << endl << "  ----------------------------------" << endl;
		cout << "  --- starting LP iteration " << slp_iter << "  ---  " << endl;
		cout << "  ---------------------------------" << endl << endl << endl;
		make_runs(jco);
		ClpSimplex model = solve_lp_problem(jco);
		update(model);
		slp_iter++;
	}
}


void sequentialLP::make_runs(Jacobian_1to1 &jco)
{
	ofstream *f_rec = &file_mgr.rec_ofstream();
	*f_rec << "  ---  calculating response matrix for iteration " << slp_iter << "  ---  " << endl;
	cout << "  ---  calculating response matrix for iteration " << slp_iter << "  ---  " << endl;
	ParamTransformSeq par_trans = opt_scenario.get_base_par_tran_seq();
	set<string> out_of_bounds;
	jco.build_runs(decision_vars, constraints_sim, opt_scenario.get_ctl_ordered_par_names(), par_trans,
		opt_scenario.get_base_group_info(), opt_scenario.get_ctl_parameter_info(),
		*run_mgr_ptr, out_of_bounds);
	jco.make_runs(*run_mgr_ptr);

	//get the base run and update simulated constraint values
	Parameters temp_pars;
	Observations temp_obs;
	//run_mgr_ptr->get_run(0, temp_pars, temp_obs,false);
	run_mgr_ptr->get_run(0, temp_pars, constraints_sim, false);
	par_trans.model2ctl_ip(temp_pars);
	//current_run.update_ctl(temp_pars, temp_obs);

	//process the remaining responses
	jco.process_runs(par_trans, opt_scenario.get_base_group_info(), *run_mgr_ptr, *null_prior, false);
	//TODO: deal with failed runs


	stringstream ss;
	ss << slp_iter << ".jcb";
	string rspmat_file = file_mgr.build_filename(ss.str());
	*f_rec << endl << "saving iteration " << slp_iter << " reponse matrix to file: " << rspmat_file << endl;
	jco.save(ss.str());

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



