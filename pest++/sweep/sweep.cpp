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

using namespace std;
using namespace pest_utils;


vector<string> prepare_parameter_csv(Parameters pars, ifstream &csv, bool forgive)
{
	if (!csv.good())
	{
		throw runtime_error("ifstream not good");
	}

	//process the header
	//any missing header labels will be marked to ignore those columns later
	string line;
	vector<string> tokens;
	if (!getline(csv, line))
		throw runtime_error("error reading header (first) line from csv file :");
	upper_ip(line);
	tokenize(line, tokens, ",", false);
	vector<string> header_tokens = tokens;

	// check for parameter names that in the pest control file but that are missing from the csv file
	vector<string> missing_names;
	string name;
	for (auto &p : pars)
	if (find(header_tokens.begin(), header_tokens.end(), p.first) == header_tokens.end())
		missing_names.push_back(p.first);

	if (missing_names.size() > 0)
	{
		stringstream ss;
		ss << " the following pest control file parameter names were not found in the parameter csv file:" << endl;
		for (auto &n : missing_names) ss << n << endl;
		if (!forgive)
			throw runtime_error(ss.str());
		else
			cout << ss.str() << endl << "continuing anyway..." << endl;
	}

	return header_tokens;
}

vector<Parameters> load_parameters_from_csv(vector<string> header_tokens, vector<string> ctl_pnames, ifstream &csv, int chunk)
{
	

	//process each parameter value line in the csv file
	int lcount = 1;
	vector<Parameters> sweep_pars;
	double val;
	string line;
	vector<string> tokens;
	Parameters pars;
	while (getline(csv, line))
	{
		tokens.clear();
		tokenize(line, tokens, ",", false);
		if (tokens.size() != header_tokens.size())
		{
			stringstream ss;
			ss << "error parsing csv file on line " << lcount << ": wrong number of tokens, ";
			ss << "expecting " << header_tokens.size() << ", found " << tokens.size();
			throw runtime_error(ss.str());
		}

		for (int i = 0; i < header_tokens.size(); i++)
		{
			// if the par name of this column is in the passed-in par names, then replace the
			// value in pars
			if (find(ctl_pnames.begin(), ctl_pnames.end(), header_tokens[i]) != ctl_pnames.end())
			{
				try
				{
					val = convert_cp<double>(tokens[i]);
				}
				catch (exception &e)
				{
					stringstream ss;
					ss << "error converting '" << tokens[i] << "' to double on line " << lcount << ": " << endl << e.what();
					throw runtime_error(ss.str());
				}
				pars[header_tokens[i]] = val;
			}
		}
		//make a cope of pars and store it
		sweep_pars.push_back(pars);
		lcount++;
		if (lcount > chunk)
			break;
	}
	//csv.close();
	return sweep_pars;
}


ofstream prep_sweep_output_file(Pest &pest_scenario)
{
	ofstream csv(pest_scenario.get_pestpp_options().get_sweep_output_csv_file());
	if (!csv.good())
	{
		throw runtime_error("could not open sweep_output_csv_file for writing: " +
			pest_scenario.get_pestpp_options().get_sweep_output_csv_file());
	}
	csv.precision(numeric_limits<double>::digits10);
	csv << "run_id";
	csv << ",phi,meas_phi,regul_phi";
	for (auto &ogrp : pest_scenario.get_ctl_ordered_obs_group_names())
	{
		csv << ',' << ogrp;
	}
	for (auto &oname : pest_scenario.get_ctl_ordered_obs_names())
		csv << ',' << oname;
	csv << endl;
	csv.flush();
	return csv;

}


void process_sweep_runs(ofstream &csv, Pest &pest_scenario, RunManagerAbstract* run_manager_ptr, vector<int> run_ids, ObjectiveFunc obj_func)
{
	Parameters pars;
	Observations obs;
	double fail_val = -1.0E+10;
	for (auto &run_id : run_ids)
	{
		csv << run_id;
		// if the run was successful
		if (run_manager_ptr->get_run(run_id, pars, obs))
		{
			map<string, double> phi_report = obj_func.phi_report(obs, pars, *(pest_scenario.get_regul_scheme_ptr()));

			csv << ',' << phi_report.at("TOTAL");
			csv << ',' << phi_report.at("MEAS");
			map<string, double>::const_iterator iregul = phi_report.find("REGUL");
			double val = 0.0;
			if (iregul != phi_report.end())
			{
				val = phi_report.at("REGUL");
			}
			csv << ',' << val;
			for (auto &obs_grp : pest_scenario.get_ctl_ordered_obs_group_names())
			{
				csv << ',' << phi_report.at(obs_grp);
			}
			for (auto oname : pest_scenario.get_ctl_ordered_obs_names())
			{
				csv << ',' << obs[oname];
			}
			csv << endl;
		}
		//if the run bombed
		else
		{
			csv << ',,,';
			for (auto &ogrp : pest_scenario.get_ctl_ordered_obs_group_names())
			{
				csv << ',';
			}
			for (int i = 0; i < pest_scenario.get_ctl_ordered_obs_names().size(); i++)
			{
				csv << ',' << fail_val;
			}
			csv << endl;
		}
	}
}



int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif
		cout << endl << endl;
		cout << "             sweep.exe - a parameteric sweep utility" << endl;
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
			cerr << "        sweep control_file.pst" << endl << endl;
			cerr << "    YAMR master:" << endl;
			cerr << "        sweep control_file.pst /H :port" << endl << endl;
			cerr << "    YAMR runner:" << endl;
			cerr << "        sweep control_file.pst /H hostname:port " << endl << endl;
			cerr << "control file pest++ options:" << endl;
			cerr << "    ++sweep_parameter_csv_file(pars_file.csv)" << endl;
			cerr << "        - csv file with each row as a par set" << endl;
			cerr << "    ++sweep_forgive(true)" << endl;
			cerr << "        - forgive control file pars missing from csv file" << endl;
			cerr << "    ++sweep_output_csv_file(output.csv)" << endl;
			cerr << "        - the csv to save run results to" << endl;
			cerr << "    ++sweep_chunk(500)" << endl;
			cerr << "        - number of runs to process in a single batch" << endl;
			cerr << "    ++sweep_base_run(true)" << endl;
			cerr << "        - run the parameter values in control file" << endl;
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
			throw runtime_error("External run manager not supported by sweep");
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
			throw runtime_error("GENIE not supported by sweep");
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
			throw runtime_error("/j option not supported by sweep");
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
			fout_rec << "             sweep.exe - a parameteric sweep utility" << endl << "for PEST(++) datasets " << endl << endl;
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


		// process the parameter csv file
		if (pest_scenario.get_pestpp_options().get_sweep_parameter_csv_file().empty())
		{
			throw runtime_error("control file pest++ type argument 'SWEEP_PARAMETER_CSV_FILE' is required for sweep");
		}

		string par_csv_file = pest_scenario.get_pestpp_options().get_sweep_parameter_csv_file();
		ifstream par_stream(par_csv_file);
		if (!par_stream.good())
		{
			throw runtime_error("could not open parameter csv file " + par_csv_file);
		}

		vector<string> header_tokens = prepare_parameter_csv(pest_scenario.get_ctl_parameters(),
			par_stream, pest_scenario.get_pestpp_options().get_sweep_forgive());
		


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

		
		// prepare the output file
		ofstream obs_stream = prep_sweep_output_file(pest_scenario);

		int chunk = pest_scenario.get_pestpp_options().get_sweep_chunk();
		vector<int> run_ids;
		vector<Parameters> sweep_pars;
		
		//if desired, add the base run to the list of runs
		if (pest_scenario.get_pestpp_options().get_sweep_base_run())
		{
			sweep_pars.push_back(pest_scenario.get_ctl_parameters());
		}
		int total_runs_done = 0;
		while (true)
		{
			//read some realizations
			sweep_pars.clear();
			try{
				performance_log.log_event("starting to read parameter csv file", 1);
				sweep_pars = load_parameters_from_csv(header_tokens, pest_scenario.get_ctl_parameters().get_keys(), 
					par_stream, chunk);
				performance_log.log_event("finished reading parameter csv file");
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << "error processing parameter csv file: " << e.what();
				performance_log.log_event(ss.str());
				fout_rec << endl << ss.str() << endl;
				fout_rec.close();

				throw runtime_error(ss.str());
			}

			// if there are no parameters to run, we are done
			if (sweep_pars.size() == 0)
				break;

			cout << "starting runs " << total_runs_done << " --> " << total_runs_done + sweep_pars.size() << endl;

			// queue up some runs
			run_ids.clear();
			for (auto &par : sweep_pars)
			{
				run_ids.push_back(run_manager_ptr->add_run(base_trans_seq.ctl2model_cp(par)));
			}

			//make some runs
			run_manager_ptr->run();

			//process the runs
			process_sweep_runs(obs_stream, pest_scenario, run_manager_ptr, run_ids, obj_func);

			total_runs_done += sweep_pars.size();

		}

		// clean up
		fout_rec.close();
		obs_stream.close();
		delete run_manager_ptr;
		cout << endl << endl << "Sweep Complete..." << endl;
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
