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
#include "RunManagerYAMR.h" //needs to be first because it includes winsock2.h
//#include <vld.h> // Memory Leak Detection using "Visual Leak Detector"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "Pest.h"
#include "Jacobian_1to1.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include "ModelRunPP.h"
#include "SVDASolver.h"
#include  "QSqrtMatrix.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "RunManagerGenie.h"
#include "RunManagerSerial.h"
#include "RunManagerExternal.h"
#include "SVD_PROPACK.h"
#include "OutputFileWriter.h"
#include "YamrSlave.h"
#include "Serialization.h"
#include "system_variables.h"
#include "pest_error.h"
#include "RestartController.h"
#include "PerformanceLog.h"
#include "debug.h"


using namespace std;
using namespace pest_utils;

//using namespace pest_utils;

int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif
		string version = "3.0.0_rc3";
		cout << endl << endl;
		cout << "             PEST++ Version " << version << endl << endl;
		cout << "                 by Dave Welter" << endl;
		cout << "     Computational Water Resource Engineering" << endl << endl << endl;
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
			cerr << "        pest++ pest_ctl_file.pst" << endl << endl;
			cerr << "    YAMR master:" << endl;
			cerr << "        pest++ control_file.pst /H :port" << endl << endl;
			cerr << "    YAMR runner:" << endl;
			cerr << "        pest++ /H hostname:port " << endl << endl;
			cerr << "    GENIE:" << endl;
			cerr << "        pest++ control_file.pst /G hostname:port" << endl << endl;
			cerr << "    external run manager:" << endl;
			cerr << "        pest++ control_file.pst /E" << endl;
			cerr << "--------------------------------------------------------" << endl;
			exit(0);
		}

		//by default use the serial run manager.  This will be changed later if another
		//run manger is specified on the command line.
		RunManagerType run_manager_type = RunManagerType::SERIAL;

		vector<string>::const_iterator it_find, it_find_next;
		string next_item;
		string socket_str = "";
		//Check for external run manager
		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/e");
		if (it_find != cmd_arg_vec.end() )
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
			tokenize(next_item, sock_parts, ":");
			try
			{
				if (sock_parts.size() != 2)
				{
					cerr << "YAMR slave requires the master be specified as /H hostname:port" << endl << endl;
					throw(PestCommandlineError(commandline));
				}
				YAMRSlave yam_slave;
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

		//string filename = get_filename(complete_path);
		string filename = complete_path;
		filename = remove_file_ext(filename); // remove .pst extension
		//string pathname = get_pathname(complete_path);
		//if (pathname.empty()) pathname = ".";
		string pathname = ".";

		RestartController restart_ctl;
		FileManager file_manager;
		//process restart and reuse jacobian directives
		vector<string>::const_iterator it_find_j = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/j");
		vector<string>::const_iterator it_find_r = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/r");
		bool restart_flag = false;
		bool save_restart_rec_header = true;

		file_manager.initialize_path(filename, pathname);
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
				fout_rec_tmp << "Restarting PEST++ ....." << endl << endl;
				cout << "    Restarting PEST++ ....." << endl << endl;
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
			fout_rec << "             PEST++ Version " << version << endl << endl;
			fout_rec << "                 by Dave Welter" << endl;
			fout_rec << "     Computational Water Resource Engineering" << endl << endl << endl;
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
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager);
			file_manager.close_file("pst");
			performance_log.log_event("finished processing control file");
		}
		catch (PestError e)
		{
			cerr << "Error prococessing control file: " << filename << endl << endl;
			cerr << e.what() << endl << endl;
			throw(e);
		}
		pest_scenario.check_inputs();
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
			run_manager_ptr = new RunManagerYAMR(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec,
				exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), port,
				file_manager.open_ofile_ext("rmr"),
				pest_scenario.get_pestpp_options().get_max_run_fail());
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
			//pest_scenario.check_par_obs();
			performance_log.log_event("finished basic model IO error checking");
			cout << "done" << endl;
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerSerial(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), pathname,
				pest_scenario.get_pestpp_options().get_max_run_fail());
		}

		const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();

		ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
		Jacobian *base_jacobian_ptr = new Jacobian_1to1(file_manager);

		TerminationController termination_ctl(pest_scenario.get_control_info().noptmax, pest_scenario.get_control_info().phiredstp,
			pest_scenario.get_control_info().nphistp, pest_scenario.get_control_info().nphinored, pest_scenario.get_control_info().relparstp,
			pest_scenario.get_control_info().nrelpar, pest_scenario.get_regul_scheme_ptr()->get_use_dynamic_reg(), pest_scenario.get_regul_scheme_ptr()->get_phimaccept());

		//if we are doing a restart, update the termination_ctl
		if (restart_flag)
		{
			restart_ctl.update_termination_ctl(termination_ctl);
		}

		SVDSolver::MAT_INV mat_inv = SVDSolver::MAT_INV::JTQJ;
		if (pest_scenario.get_pestpp_options().get_mat_inv() == PestppOptions::Q12J) mat_inv = SVDSolver::MAT_INV::Q12J;
		SVDSolver base_svd(&pest_scenario.get_control_info(), pest_scenario.get_svd_info(), &pest_scenario.get_base_group_info(),
			&pest_scenario.get_ctl_parameter_info(), &pest_scenario.get_ctl_observation_info(), file_manager,
			&pest_scenario.get_ctl_observations(), &obj_func, base_trans_seq, pest_scenario.get_prior_info_ptr(),
			*base_jacobian_ptr, pest_scenario.get_regul_scheme_ptr(),
			output_file_writer, mat_inv, &performance_log, pest_scenario.get_pestpp_options().get_base_lambda_vec(), 
			"base parameter solution", pest_scenario.get_pestpp_options().get_der_forgive());

		base_svd.set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
		//Build Super-Parameter problem
		Jacobian *super_jacobian_ptr = new Jacobian(file_manager);
		ParamTransformSeq trans_svda;
		// method must be involked as pointer as the transformation sequence it is added to will
		// take responsibility for destroying it
		TranSVD *tran_svd = new TranSVD(pest_scenario.get_pestpp_options().get_max_n_super(),
			pest_scenario.get_pestpp_options().get_super_eigthres(), "SVD Super Parameter Tranformation");

		if (pest_scenario.get_pestpp_options().get_svd_pack() == PestppOptions::PROPACK)
		{
			tran_svd->set_SVD_pack_propack();
		}

		TranFixed *tr_svda_fixed = new TranFixed("SVDA Fixed Parameter Transformation");
		trans_svda = base_trans_seq;
		trans_svda.push_back_ctl2active_ctl(tr_svda_fixed);
		trans_svda.add_custom_tran_seq(string("svda_derv2basenumeric"), trans_svda.get_ctl2active_ctl_tranformations());
		trans_svda.push_back_active_ctl2numeric(tran_svd);

		ControlInfo svd_control_info = pest_scenario.get_control_info();
		svd_control_info.relparmax = pest_scenario.get_pestpp_options().get_super_relparmax();
		// Start Solution iterations
		cout << endl << endl;
		int n_base_iter = pest_scenario.get_pestpp_options().get_n_iter_base();
		int n_super_iter = pest_scenario.get_pestpp_options().get_n_iter_super();
		int max_n_super = pest_scenario.get_pestpp_options().get_max_n_super();
		double super_eigthres = pest_scenario.get_pestpp_options().get_super_eigthres();

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
			fout_rec << "   -----    Starting PEST++ Iterations    ----    " << endl << endl << endl;
		}
		while (!termination_ctl.terminate())
		{
			//base parameter iterations
			try
			{

				if (restart_ctl.get_restart_option() != RestartController::RestartOption::NONE  && restart_ctl.get_iteration_type() == RestartController::IterationType::SUPER)
				{
					base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, false, file_manager.build_filename("jcb"));
				}
				else if (restart_ctl.get_restart_option() == RestartController::RestartOption::REUSE_JACOBIAN && n_base_iter < 0)
				{
					base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, false, file_manager.build_filename("jco"));
				}
				else if (n_base_iter < 0)
				{
					bool restart_runs = (restart_ctl.get_restart_option() == RestartController::RestartOption::REUSE_JACOBIAN);
					cur_run = base_svd.compute_jacobian(*run_manager_ptr, termination_ctl, cur_run, restart_runs);
				}
				else if (pest_scenario.get_control_info().noptmax < 0)
				{
					bool restart_runs = (restart_ctl.get_restart_option() == RestartController::RestartOption::RESUME_JACOBIAN_RUNS);
					cur_run = base_svd.compute_jacobian(*run_manager_ptr, termination_ctl, cur_run, restart_runs);
					if (restart_runs) restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
					optimum_run = cur_run;
					termination_ctl.set_terminate(true);
				}
				else if (restart_ctl.get_restart_option() == RestartController::RestartOption::REUSE_JACOBIAN)
				{
					bool calc_first_jacobian = false;
					base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, false, file_manager.build_filename("jco"));
					cur_run = base_svd.solve(*run_manager_ptr, termination_ctl, n_base_iter, cur_run, optimum_run, restart_ctl, calc_first_jacobian);
					termination_ctl.check_last_iteration();
				}
				else
				{
					if (restart_ctl.get_restart_option() == RestartController::RestartOption::RESUME_UPGRADE_RUNS)
					{
						base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, false, file_manager.build_filename("jcb"));
					}
					bool calc_first_jacobian = true;
					cur_run = base_svd.solve(*run_manager_ptr, termination_ctl, n_base_iter, cur_run, optimum_run, restart_ctl, calc_first_jacobian);
					termination_ctl.check_last_iteration();
				}
				cur_ctl_parameters = cur_run.get_ctl_pars();
				if (termination_ctl.terminate())  break;
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
			// Build Super Parameter or SVDA problem
			try
			{
				const vector<string> &nonregul_obs = pest_scenario.get_nonregul_obs();
				Parameters base_numeric_pars = base_trans_seq.ctl2numeric_cp(cur_ctl_parameters);
				const vector<string> &pars = base_numeric_pars.get_keys();
				QSqrtMatrix Q_sqrt(&(pest_scenario.get_ctl_observation_info()), &pest_scenario.get_prior_info());
				if (restart_ctl.get_restart_option() == RestartController::RestartOption::RESUME_JACOBIAN_RUNS)
				{
					//read previously computed super parameter transformation
					ifstream &fin_rst = file_manager.open_ifile_ext("rtj", ios_base::in | ios_base::binary);
					(*tran_svd).read(fin_rst);
					file_manager.close_file("rtj");
				}
				else
				{
					(*tran_svd).update_reset_frozen_pars(*base_jacobian_ptr, Q_sqrt, base_numeric_pars, max_n_super, super_eigthres, pars, nonregul_obs, cur_run.get_frozen_ctl_pars());
					(*tr_svda_fixed).reset((*tran_svd).get_frozen_derivative_pars());
				}
				SVDASolver super_svd(&svd_control_info, pest_scenario.get_svd_info(), &pest_scenario.get_base_group_info(), &pest_scenario.get_ctl_parameter_info(),
					&pest_scenario.get_ctl_observation_info(), file_manager, &pest_scenario.get_ctl_observations(), &obj_func,
					trans_svda, &pest_scenario.get_prior_info(), *super_jacobian_ptr, pest_scenario.get_regul_scheme_ptr(),
					output_file_writer, mat_inv, &performance_log, pest_scenario.get_pestpp_options().get_base_lambda_vec(), false, base_svd.get_phiredswh_flag(), base_svd.get_splitswh_flag(),
					pest_scenario.get_pestpp_options().get_max_super_frz_iter());
				super_svd.set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
				//use base jacobian to compute first super jacobian if there was not a super upgrade
				bool calc_first_jacobian = true;
				if (n_base_iter == -1)
				{ 
					//transform base jacobian to super jacobian
					super_svd.get_jacobian() = base_svd.get_jacobian();
					super_svd.get_jacobian().transform(base_trans_seq, &ParamTransformSeq::jac_numeric2active_ctl_ip);
					super_svd.get_jacobian().transform(trans_svda, &ParamTransformSeq::jac_active_ctl_ip2numeric_ip);
					//rerun base run to account for round off error in super parameters
					cur_run = super_svd.update_run(*run_manager_ptr, cur_run);
					calc_first_jacobian = false;
				}
				cur_run = super_svd.solve(*run_manager_ptr, termination_ctl, n_super_iter, cur_run, optimum_run, restart_ctl, calc_first_jacobian);
				cur_ctl_parameters = cur_run.get_ctl_pars();
				base_svd.set_phiredswh_flag(super_svd.get_phiredswh_flag());
				base_svd.set_splitswh_flag(super_svd.get_splitswh_flag());
				if (super_svd.local_iteration_terminatated() && n_base_iter == -1)
				{
					n_base_iter = 1;
				}
			}
			catch (exception &e)
			{
				cout << e.what() << endl;
				fout_rec << e.what() << endl;
				cout << "WARNING: super parameter run failed.  Switching to base parameters" << endl << endl;
				fout_rec << "WARNING: super parameter run failed.  Switching to base parameters" << endl << endl;
				ofstream &fout_restart = file_manager.get_ofstream("rst");
				RestartController::write_start_failed_super(fout_restart);
				restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
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
		map<string, double> phi_report = obj_func.phi_report(optimum_run.get_obs(), optimum_run.get_ctl_pars(), DynamicRegularization(false));
		output_file_writer.phi_report(fout_rec, termination_ctl.get_iteration_number() + 1, run_manager_ptr->get_total_runs(), phi_report, 0.0, true);
		output_file_writer.phi_report(cout, termination_ctl.get_iteration_number() + 1, run_manager_ptr->get_total_runs(), phi_report, 0.0, true);
		fout_rec << endl << endl;
		fout_rec << "Number of forward model runs performed during optimiztion: " << run_manager_ptr->get_total_runs() << endl;
		fout_rec.close();
		// clean up
		delete base_jacobian_ptr;
		delete super_jacobian_ptr;
		delete run_manager_ptr;
		cout << endl << endl << "Simulation Complete..." << endl;
		cout << flush;
#ifndef _DEBUG
	}
	catch (exception &e)
	{
		cout << "Error condition prevents further execution: " << endl << e.what() << endl;
		cout << "press enter to continue" << endl;
		char buf[256];
		OperSys::gets_s(buf, sizeof(buf));
	}
#endif
}
