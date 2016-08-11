#include <random>
#include "DifferentialEvolution.h"
#include "RunManagerAbstract.h"
#include "ModelRunPP.h"


mt19937_64 DifferentialEvolution::rand_engine = mt19937_64(1);

DifferentialEvolution::ParameterInfoDE::ParameterInfoDE(double _lower_bnd, double _upper_bnd,
	bool _log_transform)
	: lower_bnd(_lower_bnd), upper_bnd(_upper_bnd), log_transform(_log_transform)
{
}

DifferentialEvolution::DifferentialEvolution(Pest &_pest_scenario,
	FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
	const ParamTransformSeq &_par_transform, OutputFileWriter &_output_file_writer,
	PerformanceLog *_performance_log,  unsigned int seed)
	: file_manager(_file_manager), obj_func_ptr(_obj_func_ptr), par_transform(_par_transform),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	gen_1(_file_manager.build_filename("de1"))
{
	// initialize random number generator
	rand_engine.seed(seed);

	par_list = _pest_scenario.get_ctl_ordered_par_names();
	for (const auto &i : par_list)
	{
		const ParameterRec *p_info = _pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(i);
		if (p_info->tranform_type == ParameterRec::TRAN_TYPE::LOG)
		{
			parameter_info[i] = ParameterInfoDE(p_info->lbnd, p_info->ubnd, true);
		}
		else if (p_info->tranform_type == ParameterRec::TRAN_TYPE::NONE)
		{
			parameter_info[i] = ParameterInfoDE(p_info->lbnd, p_info->ubnd, false);
		}
		else
		{
			// This is a fixed or tied parameter which will be handled by a transformations
		}

	}
}

void DifferentialEvolution::solve(RunManagerAbstract &run_manager,
	RestartController &restart_controller,
	int max_gen, double f, double cr, ModelRun &cur_run)
{
	ostream &os = file_manager.rec_ofstream();

	for (int iter = 0; iter < max_gen; ++iter)
	{
		Parameters tmp_pars;
		Observations tmp_obs;

		// write header for iteration
		cout << endl;
		output_file_writer.iteration_report(cout, iter+1, run_manager.get_total_runs(), "differntial evolution");
		os << endl;
		output_file_writer.iteration_report(os, iter + 1, run_manager.get_total_runs(), "differntial evolution");
		// write initial phi report for this iteration
		bool run_target_ok = gen_1.get_run(best_run_idx, tmp_pars, tmp_obs);
		par_transform.model2ctl_ip(tmp_pars);
		map<string, double> phi_comps = obj_func_ptr->phi_report(tmp_obs, tmp_pars, DynamicRegularization::get_unit_reg_instance());
		output_file_writer.phi_report(cout, iter, run_manager.get_nruns(), phi_comps, DynamicRegularization::get_unit_reg_instance().get_weight(), false);
		output_file_writer.phi_report(os, iter, run_manager.get_nruns(), phi_comps, DynamicRegularization::get_unit_reg_instance().get_weight(), false);

		run_manager.reinitialize();
		mutation(run_manager, f, cr);

		// make trial vector model runs
		cout << endl;
		cout << "  performing trial vector model runs... ";
		cout.flush();
		run_manager.run();
		os << endl;
		best_run_idx = recombination(run_manager);
		os << endl;

		// write initial phi report for this iteration
		run_target_ok = gen_1.get_run(best_run_idx, tmp_pars, tmp_obs);
		par_transform.model2ctl_ip(tmp_pars);
		phi_comps = obj_func_ptr->phi_report(tmp_obs, tmp_pars, DynamicRegularization::get_unit_reg_instance());
		output_file_writer.phi_report(cout, iter, run_manager.get_nruns(), phi_comps, DynamicRegularization::get_unit_reg_instance().get_weight(), true);
		cout << endl;
		output_file_writer.phi_report(os, iter, run_manager.get_nruns(), phi_comps, DynamicRegularization::get_unit_reg_instance().get_weight(), true);
		os << endl;
	}
}


void DifferentialEvolution::initialize_population(RunManagerAbstract &run_manager, int _d)
{
	d = _d;
	Parameters ctl_pars;
	for (int i = 0; i < d; ++i)
	{
		ctl_pars.clear();
		initialize_vector(ctl_pars);
		par_transform.ctl2model_ip(ctl_pars);
		run_manager.add_run(ctl_pars);
	}
	// make innitial population vector model runs
	cout << endl;
	cout << "  performing initial population model runs... ";
	cout.flush();
	run_manager.run();
	gen_1.copy(run_manager.get_runstorage_ref());

	// get the best_run to track phi
	int d = gen_1.get_nruns();
	int r_status;
	int n_par = par_list.size();
	double best_phi = std::numeric_limits<double>::max();

	Parameters tmp_pars;
	Observations tmp_obs;
	ModelRun tmp_run(obj_func_ptr);
	for (int i_run = 0; i_run < d; ++i_run)
	{
		bool r_status = gen_1.get_run(i_run, tmp_pars, tmp_obs);
		if (r_status > 0)
		{
			par_transform.model2ctl_ip(tmp_pars);
			tmp_run.update_ctl(tmp_pars, tmp_obs);
			double tmp_phi = tmp_run.get_phi(DynamicRegularization::get_unit_reg_instance());
			if (tmp_phi < best_phi)
			{
				best_phi = tmp_phi;
				best_run_idx = i_run;
			}
		}
	}

}

void DifferentialEvolution::initialize_vector(Parameters &ctl_pars)
{
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	
	for (const auto &i : parameter_info)
	{
		double p_val; 
		const string &p_name = i.first;
		double rnum = rand_engine();
		// use uniform distribution to initialize parameters
		if (i.second.log_transform)
		{
			double r_num = distribution(rand_engine);
			p_val = log10(i.second.lower_bnd) + r_num * (log10(i.second.upper_bnd) - log10(i.second.lower_bnd));
			p_val = pow(10.0, p_val);
			ctl_pars.insert(p_name, p_val);
		}
		else
		{
			p_val = i.second.lower_bnd + distribution(rand_engine) * (i.second.upper_bnd - i.second.lower_bnd);
			ctl_pars.insert(p_name, p_val);
		}
	}
}

void DifferentialEvolution::mutation(RunManagerAbstract &run_manager, double f, double cr)
{
	int d = gen_1.get_nruns();
	int r_status;
	int n_par = par_list.size();
	std::uniform_int_distribution<int> uni_par(0, n_par-1);
	std::uniform_real_distribution<double> cr_prob(0, 1);
	vector<int> successful_run_ids;
	// generate a vector of successful runs
	for (int i_run = 0; i_run < d; ++i_run)
	{
		r_status = gen_1.get_run_status(i_run);
		if (r_status > 0)
		{
			successful_run_ids.push_back(i_run);
		}
	}
	int d_ok = successful_run_ids.size();
	std::uniform_int_distribution<int> uni_run_ok(0, d_ok-1);

	Parameters xa;
	Parameters xb;
	Parameters xc;
	Parameters x_trial;
	for (int i_run = 0; i_run < d; ++i_run)
	{
		int xa_id = successful_run_ids[uni_run_ok(rand_engine)];
		int xb_id = successful_run_ids[uni_run_ok(rand_engine)];
		while (xa_id == xb_id)
		{
			xb_id = successful_run_ids[uni_run_ok(rand_engine)];
		}
		int xc_id = successful_run_ids[uni_run_ok(rand_engine)];
		//initialize trail vector with the traget vector
		gen_1.get_parameters(i_run, x_trial);
		gen_1.get_parameters(xa_id, xa);
		gen_1.get_parameters(xb_id, xb);
		gen_1.get_parameters(xc_id, xb);
		par_transform.model2ctl_ip(xa);
		par_transform.model2ctl_ip(xb);
		par_transform.model2ctl_ip(xc);

		int par_id_chg = uni_par(rand_engine);
		for (int idx=0; idx<n_par; ++idx)
		{
			const string &ipar = par_list[idx];
			double a = xa[ipar];
			double b = xb[ipar];
			double c = xc[ipar];
			double delta = a - b;
			double c_p = c + f * (delta);
			double p_min = parameter_info[ipar].lower_bnd;
			double p_max = parameter_info[ipar].upper_bnd;
			// reflect purturbation if parameter is outside it's bounds
			while (c_p < p_min || c_p > p_max)
			{
				if (c_p < p_min)
				{
					c_p += p_min - c_p;
				}
				if (c_p > p_max)
				{
					c_p -= c_p - p_max;
				}
			}
			// do cross over
			//the trial vector was initialized with the target vector
			// so the parameters from the parent don't need to be set 
			double rand_cr = cr_prob(rand_engine);
			if (rand_cr > cr || idx == par_id_chg)
			{
				x_trial[ipar] = c_p;
			}
		}
		par_transform.ctl2model_ip(x_trial);
		run_manager.add_run(x_trial);
	}
}

int DifferentialEvolution::recombination(RunManagerAbstract &run_manager)
{
	ostream &os = file_manager.rec_ofstream();

	int best_run_idx = 0;
	ModelRun run_target(obj_func_ptr);
	ModelRun run_canidate(obj_func_ptr);
	Parameters tmp_pars_targ;
	Observations tmp_obs_targ;
	Parameters tmp_pars_can;
	Observations tmp_obs_can;

	double best_phi = std::numeric_limits<double>::max();

	for (int i_run = 0; i_run < d; ++i_run)
	{
		double new_phi = std::numeric_limits<double>::max();
		bool run_target_ok = gen_1.get_run(i_run, tmp_pars_targ, tmp_obs_targ);
		bool  run_canidate_ok = run_manager.get_run(i_run, tmp_pars_can, tmp_obs_can);
		if (!run_canidate_ok && !run_canidate_ok)
		{
			//keep current target
			new_phi = std::numeric_limits<double>::max();

		}
		if (!run_canidate_ok)
		{
			//keep current target
			run_canidate.update_ctl(tmp_pars_targ, tmp_obs_targ);
			new_phi = run_target.get_phi(DynamicRegularization::get_unit_reg_instance());
		}
		else if (!run_target_ok)
		{
			gen_1.update_run(i_run, tmp_pars_can, tmp_obs_can);
			new_phi = run_canidate.get_phi(DynamicRegularization::get_unit_reg_instance());
		}
		else
		{
			par_transform.model2ctl_ip(tmp_pars_targ);
			run_target.update_ctl(tmp_pars_targ, tmp_obs_targ);
			double phi_target = run_target.get_phi(DynamicRegularization::get_unit_reg_instance());
			par_transform.model2ctl_ip(tmp_pars_can);
			run_canidate.update_ctl(tmp_pars_can, tmp_obs_can);
			double phi_canidate = run_canidate.get_phi(DynamicRegularization::get_unit_reg_instance());
			new_phi = min(phi_target, phi_canidate);
			os << "  id = " << i_run << ";  parent phi = " << phi_target << ";  child phi = " << phi_canidate << endl;
			if (phi_canidate < phi_target)
			{
				gen_1.update_run(i_run, tmp_pars_can, tmp_obs_can);
			}
		}
		if (new_phi < best_phi)
		{
			best_phi = new_phi;
			best_run_idx = i_run;
		}
	}
	return best_run_idx;
}

DifferentialEvolution::~DifferentialEvolution()
{
}
