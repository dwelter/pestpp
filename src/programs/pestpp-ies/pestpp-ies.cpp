// pestpp-ies.cpp : Defines the entry point for the console application.
//

#include "RunManagerYAMR.h" //needs to be first because it includes winsock2.h
//#include <vld.h> // Memory Leak Detection using "Visual Leak Detector"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "config_os.h"
#include "Pest.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include "ModelRunPP.h"
#include "FileManager.h"
#include "RunManagerSerial.h"
#include "OutputFileWriter.h"
#include "YamrSlave.h"
#include "Serialization.h"
#include "system_variables.h"
#include "pest_error.h"
#include "RestartController.h"
#include "PerformanceLog.h"
#include "debug.h"
#include "logger.h"
#include "Ensemble.h"

using namespace std;
using namespace pest_utils;


int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif
		cout << endl << endl;
		cout << "             pestpp-ies.exe - a GLM iterative ensemble smoother" << endl;
		cout << "                     for PEST(++) datasets " << endl << endl;
		cout << "                 by the PEST++ developement team" << endl << endl << endl;
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
			cerr << "        pestpp-ies control_file.pst" << endl << endl;
			cerr << "    YAMR master:" << endl;
			cerr << "        pestpp-ies control_file.pst /H :port" << endl << endl;
			cerr << "    YAMR runner:" << endl;
			cerr << "        pestpp-ies control_file.pst /H hostname:port " << endl << endl;
			cerr << "control file pest++ options:" << endl;
			cerr << "    ++ies_par_csv(pars_file.csv)" << endl;
			cerr << "        - csv file with each row as a parameter realization" << endl;
			cerr << "    ++ies_obs_csv(obs_file.csv)" << endl;
			cerr << "        - csv file with each row as an observation realization" << endl;
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
			throw runtime_error("External run manager not supported by pestpp-ies");
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
			throw runtime_error("GENIE not supported by pestpp-ies");
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
			throw runtime_error("/j option not supported by pestpp-ies");
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
			fout_rec_tmp << "Restarting sweep ....." << endl << endl;
			cout << "    Restarting sweep ....." << endl << endl;

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
			fout_rec << "             pestpp-ies.exe - a GLM iterative Ensemble Smoother" << endl << "for PEST(++) datasets " << endl << endl;
			fout_rec << "                 by the PEST++ developement team" << endl << endl << endl;
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
			fout_rec << "Error prococessing control file: " << filename << endl << endl;
			fout_rec << e.what() << endl << endl;
			fout_rec.close();
			throw(e);
		}
		pest_scenario.check_inputs(fout_rec);

		PestppOptions ppopt = pest_scenario.get_pestpp_options();

		fout_rec << "    pestpp-ies parameter csv file = " << left << setw(50) << ppopt.get_ies_par_csv() << endl;
		fout_rec << "    pestpp-ies observation csv file = " << left << setw(50) << ppopt.get_ies_obs_csv() << endl;
		
		//Initialize OutputFileWriter to handle IO of suplementary files (.par, .par, .svd)
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

		ParameterEnsemble pe(base_trans_seq, pest_scenario, file_manager, output_file_writer, &performance_log);
		pe.from_csv(pest_scenario.get_pestpp_options().get_ies_par_csv());

		cout << *pe.get_mean_diff().eptr() << endl;


		ObservationEnsemble oe(&obj_func, pest_scenario, file_manager, output_file_writer, &performance_log);
		oe.from_csv(pest_scenario.get_pestpp_options().get_ies_obs_csv());

		ObservationEnsemble oe_org = oe;

		EnsemblePair epair(pe, oe);

		epair.run(run_manager_ptr);

		
		/*vector<string>obs_names = pest_scenario.get_ctl_ordered_obs_names();
		string last = obs_names[obs_names.size()-1];
		obs_names.pop_back();
		obs_names.pop_back();
		obs_names.insert(obs_names.begin(), last);
		obs_names.insert(obs_names.begin(), last);
		obs_names.insert(obs_names.begin(), last);*/

		/*vector<string> obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
		cout << oe.get_eigen(pe.get_real_names(),obs_names) << endl;
		cout << oe.get_reals() << endl;
		cout << oe.get_eigen(epair.get_oe_ptr()->get_real_names(), obs_names) - oe.get_eigen(oe.get_real_names(), obs_names) << endl;
*/
		//cout << epair.get_active_oe_eigen() << endl;
		//cout << epair.get_active_pe_eigen() << endl;

		// clean up
		fout_rec.close();
		delete run_manager_ptr;
		cout << endl << endl << "pestpp-ies analysis complete..." << endl;
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