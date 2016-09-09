
#include "RunManagerYAMR.h" //needs to be first because it includes winsock2.h
//#include <vld.h> // Memory Leak Detection using "Visual Leak Detector"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "config_os.h"
#include "Pest.h"
#include "Jacobian_1to1.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include "ModelRunPP.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "RunManagerGenie.h"
#include "RunManagerSerial.h"
#include "RunManagerExternal.h"
#include "OutputFileWriter.h"
#include "YamrSlave.h"
#include "Serialization.h"
#include "system_variables.h"
#include "pest_error.h"
#include "RestartController.h"
#include "PerformanceLog.h"
#include "debug.h"
#include "DifferentialEvolution.h"

#include "linear_analysis.h"
#include "logger.h"
#include "covariance.h"
#include "sequential_lp.h"

using namespace std;
using namespace pest_utils;

//using namespace pest_utils;

int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif
		string version = PESTPP_VERSION;
		cout << endl << endl;
		cout << "             pestpp-opt version " << version << endl << endl;
		cout << "                 by the PEST++ development team" << endl;
		// build commandline
		string commandline = "";
		for (int i = 0; i < argc; ++i)
		{
			commandline.append(" ");
			commandline.append(argv[i]);
		}

		vector<string> cmd_arg_vec(argc);
		copy(argv, argv + argc, cmd_arg_vec.begin());
		for (vector<string>::iterator it = cmd_arg_vec.begin(); it != cmd_arg_vec.end(); ++it)
		{
			transform(it->begin(), it->end(), it->begin(), ::tolower);
		}

		string complete_path;
		enum class RunManagerType { SERIAL, YAMR, GENIE, EXTERNAL };

		if (argc >= 2) {
			complete_path = argv[1];
		}
		else {
			cerr << "--------------------------------------------------------" << endl;
			cerr << "usage:" << endl << endl;
			cerr << "    serial run manager:" << endl;
			cerr << "        pestpp-opt pest_ctl_file.pst" << endl << endl;
			cerr << "    YAMR master:" << endl;
			cerr << "        pestpp-opt control_file.pst /H :port" << endl << endl;
			cerr << "    YAMR runner:" << endl;
			cerr << "        pestpp-opt /H hostname:port " << endl << endl;
			cerr << "    GENIE:" << endl;
			cerr << "        pestpp-opt control_file.pst /G hostname:port" << endl << endl;
			cerr << "    external run manager:" << endl;
			cerr << "        pestpp-opt control_file.pst /E" << endl;
			cerr << "--------------------------------------------------------" << endl;
			exit(0);
		}


		FileManager file_manager;
		string filename = complete_path;
		string pathname = ".";
		file_manager.initialize_path(get_filename_without_ext(filename), pathname);

		//by default use the serial run manager.  This will be changed later if another
		//run manger is specified on the command line.
		RunManagerType run_manager_type = RunManagerType::SERIAL;

		vector<string>::const_iterator it_find, it_find_next;
		string next_item;
		string socket_str = "";
		//Check for external run manager
		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/e");
		if (it_find != cmd_arg_vec.end())
		{
			run_manager_type = RunManagerType::EXTERNAL;
		}
		//Check for YAMR Slave
		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/h");
		next_item.clear();
		if (it_find != cmd_arg_vec.end() && it_find + 1 != cmd_arg_vec.end())
		{
			next_item = *(it_find + 1);
			strip_ip(next_item);
		}
		if (it_find != cmd_arg_vec.end() && !next_item.empty() && next_item[0] != ':')
		{
			// This is a YAMR Slave, start PEST++ as a YAMR Slave
			vector<string> sock_parts;
			vector<string>::const_iterator it_find_yamr_ctl;
			string file_ext = get_filename_ext(filename);
			tokenize(next_item, sock_parts, ":");
			try
			{
				if (sock_parts.size() != 2)
				{
					cerr << "YAMR slave requires the master be specified as /H hostname:port" << endl << endl;
					throw(PestCommandlineError(commandline));
				}
				YAMRSlave yam_slave;
				string ctl_file = "";
				try {
					string ctl_file;
					if (upper_cp(file_ext) == "YMR")
					{
						ctl_file = file_manager.build_filename("ymr");
						yam_slave.process_yamr_ctl_file(ctl_file);
					}
					else
					{
						// process traditional PEST control file
						ctl_file = file_manager.build_filename("pst");
						yam_slave.process_ctl_file(ctl_file);
					}
				}
				catch (PestError e)
				{
					cerr << "Error prococessing control file: " << ctl_file << endl << endl;
					cerr << e.what() << endl << endl;
					throw(e);
				}

				yam_slave.start(sock_parts[0], sock_parts[1]);
			}
			catch (PestError &perr)
			{
				cerr << perr.what();
				throw(perr);
			}
			cout << endl << "Simulation Complete..." << endl;
			exit(0);
		}
		//Check for YAMR Master
		else if (it_find != cmd_arg_vec.end())
		{
			// using YAMR run manager
			run_manager_type = RunManagerType::YAMR;
			socket_str = next_item;
		}

		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/g");
		next_item.clear();
		if (it_find != cmd_arg_vec.end() && it_find + 1 != cmd_arg_vec.end())
		{
			next_item = *(it_find + 1);
			strip_ip(next_item);

		}
		//Check for GENIE Master
		if (it_find != cmd_arg_vec.end())
		{
			//Using GENIE run manager
			run_manager_type = RunManagerType::GENIE;
			socket_str = next_item;
		}

		RestartController restart_ctl;

		//process restart and reuse jacobian directives
		vector<string>::const_iterator it_find_j = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/j");
		vector<string>::const_iterator it_find_r = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/r");
		bool restart_flag = false;
		bool save_restart_rec_header = true;

		debug_initialize(file_manager.build_filename("dbg"));
		if (it_find_j != cmd_arg_vec.end())
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::REUSE_JACOBIAN;
			file_manager.open_default_files();
		}
		else if (it_find_r != cmd_arg_vec.end())
		{
			ifstream &fin_rst = file_manager.open_ifile_ext("rst");
			restart_ctl.process_rst_file(fin_rst);
			file_manager.close_file("rst");
			restart_flag = true;
			file_manager.open_default_files(true);
			ofstream &fout_rec_tmp = file_manager.rec_ofstream();
			fout_rec_tmp << endl << endl;
			if (run_manager_type == RunManagerType::EXTERNAL)
			{
				save_restart_rec_header = false;
			}
			else
			{
				fout_rec_tmp << "Restarting pestpp-opt ....." << endl << endl;
				cout << "    Restarting pestpp-opt ....." << endl << endl;
			}
		}
		else
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
			file_manager.open_default_files();
		}

		ofstream &fout_rec = file_manager.rec_ofstream();
		PerformanceLog performance_log(file_manager.open_ofile_ext("pfm"));

		if (!restart_flag || save_restart_rec_header)
		{
			fout_rec << "             pestpp-opt version " << version << endl << endl;
			fout_rec << "                 by the petspp development team" << endl;
			fout_rec << endl;
			fout_rec << "using control file: \"" << complete_path << "\"" << endl << endl;
		}

		cout << endl;
		cout << "using control file: \"" << complete_path << "\"" << endl << endl;

		// create pest run and process control file to initialize it
		Pest pest_scenario;
		pest_scenario.set_defaults();

		try {
			performance_log.log_event("starting to process control file", 1);
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"));
			file_manager.close_file("pst");
			performance_log.log_event("finished processing control file");
		}
		catch (PestError e)
		{
			cerr << "Error prococessing control file: " << filename << endl << endl;
			cerr << e.what() << endl << endl;
			throw(e);
		}
		pest_scenario.check_inputs(fout_rec);

		//if base jco arg read from control file, reset restart controller
		if (!pest_scenario.get_pestpp_options().get_basejac_filename().empty())
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::REUSE_JACOBIAN;
		}

		//Initialize OutputFileWriter to hadle IO of suplementary files (.par, .par, .svd)
		//bool save_eign = pest_scenario.get_svd_info().eigwrite > 0;	
		OutputFileWriter output_file_writer(file_manager, pest_scenario, restart_flag);
		output_file_writer.set_svd_output_opt(pest_scenario.get_svd_info().eigwrite);
		if (!restart_flag)
		{
			output_file_writer.scenario_report(fout_rec);
		}
		if (pest_scenario.get_pestpp_options().get_iter_summary_flag())
		{
			output_file_writer.write_par_iter(0, pest_scenario.get_ctl_parameters());
		}
		RunManagerAbstract *run_manager_ptr;
		if (run_manager_type == RunManagerType::YAMR)
		{
			string port = socket_str;
			strip_ip(port);
			strip_ip(port, "front", ":");
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerYAMR(
				file_manager.build_filename("rns"), port,
				file_manager.open_ofile_ext("rmr"),
				pest_scenario.get_pestpp_options().get_max_run_fail(),
				pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_fac());
		}
		else if (run_manager_type == RunManagerType::GENIE)
		{
			strip_ip(socket_str);
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerGenie(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), socket_str);
		}
		else if (run_manager_type == RunManagerType::EXTERNAL)
		{
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerExternal(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), file_manager.build_filename("ext"),
				file_manager.build_filename("exi"),
				pest_scenario.get_pestpp_options().get_max_run_fail());
		}
		else
		{
			performance_log.log_event("starting basic model IO error checking", 1);
			cout << "checking model IO files...";
			pest_scenario.check_io();
			performance_log.log_event("finished basic model IO error checking");
			cout << "done" << endl;
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerSerial(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), pathname,
				pest_scenario.get_pestpp_options().get_max_run_fail());
		}

		//setup the parcov, if needed
		Covariance parcov;
		//if (pest_scenario.get_pestpp_options().get_use_parcov_scaling())
		double parcov_scale_fac = pest_scenario.get_pestpp_options().get_parcov_scale_fac();
		if (parcov_scale_fac > 0.0)
		{
			parcov.try_from(pest_scenario, file_manager);
		}
		const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();

		ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
		Jacobian *base_jacobian_ptr = new Jacobian_1to1(file_manager);

		TerminationController termination_ctl(pest_scenario.get_control_info().noptmax, pest_scenario.get_control_info().phiredstp,
			pest_scenario.get_control_info().nphistp, pest_scenario.get_control_info().nphinored, pest_scenario.get_control_info().relparstp,
			pest_scenario.get_control_info().nrelpar, pest_scenario.get_regul_scheme_ptr()->get_use_dynamic_reg(),
			pest_scenario.get_regul_scheme_ptr()->get_phimaccept(), pest_scenario.get_pestpp_options().get_reg_frac());

		//if we are doing a restart, update the termination_ctl
		if (restart_flag)
		{
			restart_ctl.update_termination_ctl(termination_ctl);
		}

		Parameters cur_ctl_parameters = pest_scenario.get_ctl_parameters();
		//Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
		//Neither of these will change over the course of the simulation

		if (restart_ctl.get_restart_option() == RestartController::RestartOption::RESUME_JACOBIAN_RUNS)
		{
			run_manager_ptr->initialize_restart(file_manager.build_filename("rnj"));
		}
		else
		{
			run_manager_ptr->initialize(base_trans_seq.ctl2model_cp(cur_ctl_parameters), pest_scenario.get_ctl_observations());
		}


		ModelRun optimum_run(&obj_func, pest_scenario.get_ctl_observations());


		// if noptmax=0 make one run with the intital parameters
		if (pest_scenario.get_control_info().noptmax == 0) {
			Parameters init_model_pars = base_trans_seq.ctl2model_cp(cur_ctl_parameters);
			optimum_run.set_ctl_parameters(init_model_pars);
			run_manager_ptr->reinitialize();
			run_manager_ptr->add_run(init_model_pars);
			try
			{
				run_manager_ptr->run();
			}
			catch (exception &e)
			{
				cout << "Model run failed.  No results were recorded." << endl << e.what() << endl;
				fout_rec << "Model run failed.  No results were recorded." << endl << e.what() << endl;
				exit(1);
			}
			Parameters tmp_pars;
			Observations tmp_obs;
			bool success = run_manager_ptr->get_run(0, tmp_pars, tmp_obs);
			base_trans_seq.model2ctl_ip(tmp_pars);
			//termination_ctl.set_terminate(true);

			if (success)
			{
				termination_ctl.check_last_iteration();
				optimum_run.update_ctl(tmp_pars, tmp_obs);
				// save parameters to .par file
				output_file_writer.write_par(file_manager.open_ofile_ext("par"), optimum_run.get_ctl_pars(), *(base_trans_seq.get_offset_ptr()),
					*(base_trans_seq.get_scale_ptr()));
				file_manager.close_file("par");
				// save new residuals to .rei file
				output_file_writer.write_rei(file_manager.open_ofile_ext("rei"), 0,
					*(optimum_run.get_obj_func_ptr()->get_obs_ptr()),
					optimum_run.get_obs(), *(optimum_run.get_obj_func_ptr()),
					optimum_run.get_ctl_pars());
				file_manager.close_file("rei");
				run_manager_ptr->free_memory();
			}
			else
			{
				cout << "Model run failed.  No results were recorded." << endl << endl;
				fout_rec << "Model run failed.  No results were recorded." << endl << endl;
				exit(1);
			}
			termination_ctl.set_terminate(true);
		}


		//// Differential Evolution
		//if (pest_scenario.get_pestpp_options().get_global_opt() == PestppOptions::OPT_DE)
		//{
		//	int rand_seed = 1;
		//	int np = pest_scenario.get_pestpp_options().get_de_npopulation();
		//	int max_gen = pest_scenario.get_pestpp_options().get_de_max_gen();
		//	double f = pest_scenario.get_pestpp_options().get_de_f();
		//	double cr = pest_scenario.get_pestpp_options().get_de_cr();
		//	bool dither_f = pest_scenario.get_pestpp_options().get_de_dither_f();
		//	ModelRun init_run(&obj_func, pest_scenario.get_ctl_observations());
		//	Parameters cur_ctl_parameters = pest_scenario.get_ctl_parameters();
		//	run_manager_ptr->reinitialize();
		//	DifferentialEvolution de_solver(pest_scenario, file_manager, &obj_func,
		//		base_trans_seq, output_file_writer, &performance_log, rand_seed);
		//	de_solver.initialize_population(*run_manager_ptr, np);
		//	de_solver.solve(*run_manager_ptr, restart_ctl, max_gen, f, cr, dither_f, init_run);
		//	run_manager_ptr->free_memory();
		//	exit(1);
		//}


		//Define model Run for Base Parameters (uses base parameter tranformations)
		ModelRun cur_run(&obj_func, pest_scenario.get_ctl_observations());

		cur_run.set_ctl_parameters(cur_ctl_parameters);
		//If this is a restart we need to get the latest ctl parameters
		if (restart_ctl.get_restart_option() != RestartController::RestartOption::NONE)
		{
			Parameters restart_pars = restart_ctl.get_restart_parameters(file_manager.build_filename("parb"), file_manager.build_filename("par"));
			if (restart_pars.size() > 0)
			{
				cur_run.set_ctl_parameters(restart_pars);
			}
		}
		if (!restart_flag || save_restart_rec_header)
		{
			fout_rec << "   -----    Starting Optimization Iterations    ----    " << endl << endl;
		}
		while (!termination_ctl.terminate())
		{
			//base parameter iterations
			try
			{
				
			}
			catch (exception &e)
			{
				cout << endl << endl;;
				cout << e.what() << endl;
				fout_rec << endl << endl;
				fout_rec << e.what() << endl;
				cout << "FATAL ERROR: base parameter run failed." << endl;
				fout_rec << "FATAL ERROR: base parameter run failed." << endl;
				exit(1);
			}
			
		}
		cout << endl;
		termination_ctl.termination_summary(cout);
		cout << endl;
		termination_ctl.termination_summary(fout_rec);
		fout_rec << endl;
		cout << "FINAL OPTIMISATION RESULTS" << endl << endl;
		fout_rec << "FINAL OPTIMISATION RESULTS" << endl << endl;

		fout_rec << "  Optimal parameter values  " << endl;
		output_file_writer.par_report(fout_rec, optimum_run.get_ctl_pars());

		fout_rec << endl << "  Observations with optimal model-simulated equivalents and residuals" << endl;
		output_file_writer.obs_report(fout_rec, *obj_func.get_obs_ptr(), optimum_run.get_obs(), obj_func);

		fout_rec << endl << "Final composite objective function " << endl;
		map<string, double> phi_report = obj_func.phi_report(optimum_run.get_obs(), optimum_run.get_ctl_pars(), *(pest_scenario.get_regul_scheme_ptr()));
		output_file_writer.phi_report(fout_rec, termination_ctl.get_iteration_number() + 1, run_manager_ptr->get_total_runs(), phi_report, 0.0, true);
		output_file_writer.phi_report(cout, termination_ctl.get_iteration_number() + 1, run_manager_ptr->get_total_runs(), phi_report, 0.0, true);
		fout_rec << endl << endl;
		fout_rec << "Number of forward model runs performed during optimiztion: " << run_manager_ptr->get_total_runs() << endl;

		// clean up
		fout_rec.close();
		delete base_jacobian_ptr;
		delete run_manager_ptr;
		cout << endl << endl << "Simulation Complete..." << endl;
		cout << flush;
#ifndef _DEBUG
	}
	catch (exception &e)
	{
		cout << "Error condition prevents further execution: " << endl << e.what() << endl;
		//cout << "press enter to continue" << endl;
		//char buf[256];
		//OperSys::gets_s(buf, sizeof(buf));
	}
#endif
}
