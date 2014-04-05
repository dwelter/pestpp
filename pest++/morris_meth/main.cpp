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
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MorrisMethod.h"
#include "sobol.h"
#include "Pest.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include "ModelRunPP.h"
#include "FileManager.h"
#include "RunManagerGenie.h"
#include "RunManagerSerial.h"
#include "OutputFileWriter.h"
#include "YamrSlave.h"
#include "Serialization.h"
#include "system_variables.h"



using namespace std;
using namespace pest_utils;
using Eigen::MatrixXd;
using Eigen::VectorXd;


int main(int argc, char* argv[])
{
	string version = "2.1.0";
	cout << endl << endl;
	cout << "             GSA++ Version " << version << endl << endl;
	cout << "                 by Dave Welter" << endl;
	cout << "     Computational Water Resource Engineering"<< endl << endl << endl;
	// build commandline
	string commandline = "";
	for(int i=0; i<argc; ++i)
	{
		commandline.append(" ");
		commandline.append(argv[i]);
	}
	string complete_path;
	enum class RunManagerType {SERIAL, YAMR, GENIE};
	string socket_str;

	if (argc >=2) {
		complete_path = argv[1];
	}
	else {
		cerr << "--------------------------------------------------------" << endl;
		cerr << "usage:" << endl << endl;
		cerr << "    serial run manager:" << endl;
		cerr << "        gsa pest_ctl_file.pst" << endl << endl;
		cerr << "    YAMR master:" << endl;
		cerr << "        gsa control_file.pst /H :port" << endl; 
		cerr << "    YAMR runner:" << endl;
		cerr << "        gsa /H hostname:port " << endl << endl;
		cerr << "    GENIE:" << endl;
		cerr << "        gsa control_file.pst /G hostname:port" << endl;
		cerr << "--------------------------------------------------------" << endl;
		exit(0);
	}

	// testing
	ifstream fin_junk;
	fin_junk.open("C:\\Users\\dwelter\\Desktop\\junk.txt");
	MorrisMethod::process_pooled_var_file(fin_junk);


	// This is a YAMR Slave, start PEST++ as a YAMR Slave
	if (argc >=3 && upper(argv[1]) == "/H") {
		string socket_str = argv[2];
		strip_ip(socket_str);
		vector<string> sock_parts;
		tokenize(socket_str, sock_parts, ":");
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
		}
		cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
		char buf[256];
		OperSys::gets_s(buf, sizeof(buf));
		exit(0);
	}

	RunManagerType run_manager_type = RunManagerType::SERIAL;
	// Start PEST++ using YAMR run manager
	if (argc >=4 &&  upper(argv[2]) == "/H") {
		run_manager_type = RunManagerType::YAMR;
	}
	else if (argc >=4 &&  upper(argv[2]) == "/G") {
		run_manager_type = RunManagerType::GENIE;
	}


	string filename = get_filename(complete_path);
	filename = remove_file_ext(filename); // remove .pst extension
	string pathname = get_pathname(complete_path);
	if (pathname.empty()) pathname = ".";
	FileManager file_manager(filename, pathname);

	ofstream &fout_rec = file_manager.rec_ofstream();
	fout_rec << "             GSA++ Version " << version << endl << endl;
	fout_rec << "                 by Dave Welter" << endl;
	fout_rec << "     Computational Water Resource Engineering"<< endl << endl << endl;

	// create pest run and process control file to initialize it
	Pest pest_scenario;
	pest_scenario.set_defaults();
	try {
		pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager);
		file_manager.close_file("pst");
	}
	catch(PestError e)
	{
		cerr << "Error prococessing control file: " << filename << endl << endl;
		cerr << e.what() << endl << endl;
		throw(e);
	}
	pest_scenario.check_inputs();

	RunManagerAbstract *run_manager_ptr;
	if (run_manager_type == RunManagerType::YAMR)
	{
		string port = argv[3];
		strip_ip(port);
		strip_ip(port, "front", ":");
		const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
		run_manager_ptr = new RunManagerYAMR (exi.comline_vec,
		exi.tplfile_vec, exi.inpfile_vec,
		exi.insfile_vec, exi.outfile_vec,
		file_manager.build_filename("rns"), port,
		file_manager.open_ofile_ext("rmr"));
	}
	else if (run_manager_type == RunManagerType::GENIE)
	{
		string socket_str = argv[3];
		strip_ip(socket_str);
		const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
		run_manager_ptr = new RunManagerGenie(exi.comline_vec,
		exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
		file_manager.build_filename("rns"), socket_str);
	}
	else
	{
		const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
		run_manager_ptr = new RunManagerSerial(exi.comline_vec,
		exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
		file_manager.build_filename("rns"), pathname);
	}

	cout << endl;
	fout_rec << endl;
	cout << "using control file: \"" <<  complete_path << "\"" << endl;
	fout_rec << "using control file: \"" <<  complete_path << "\"" << endl;

	map<string, string> gsa_opt_map;
	//process .gsa file
	try
	{
		gsa_opt_map = GsaAbstractBase::process_gsa_file(file_manager.open_ifile_ext("gsa"), file_manager);
		file_manager.close_file("gsa");
	}
	catch(PestError e)
	{
		cerr << "Error prococessing .gsa file: " << filename << endl << endl;
		cerr << e.what() << endl << endl;
		throw(e);
	}

	//Build Transformation with ctl_2_numberic
	ParamTransformSeq base_partran_seq(pest_scenario.get_base_par_tran_seq());
	Parameters ctl_par = pest_scenario.get_ctl_parameters();

	//Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
	//Neither of these will change over the course of the simulation
	run_manager_ptr->initialize(base_partran_seq.ctl2model_cp(ctl_par), pest_scenario.get_ctl_observations());

	// Get the lower and upper bounds of the parameters
	//remove fixed parameters from ctl_par
	vector<string> adj_par_name_vec;
	Parameters fixed_pars;
	for (auto &i : pest_scenario.get_ctl_ordered_par_names())
	{
		if(pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(i)->tranform_type != ParameterRec::TRAN_TYPE::FIXED)
		{
			adj_par_name_vec.push_back(i);
		}
		else
		{
			fixed_pars.insert(i, pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(i)->init_value);
		}
	}


	Parameters lower_bnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(ctl_par.get_keys());
	Parameters upper_bnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(ctl_par.get_keys());

	//Build Transformation with ctl_2_numberic
	ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
	ModelRun model_run(&obj_func, pest_scenario.get_ctl_observations());
	const set<string> &log_trans_pars = base_partran_seq.get_log10_ptr()->get_items();
	auto method = gsa_opt_map.find("METHOD");

	GsaAbstractBase* gsa_method = nullptr;
	if (method != gsa_opt_map.end() && method->second == "MORRIS")
	{
		int morris_r = 4;
		int morris_p = 5;
		auto morris_r_it = gsa_opt_map.find("MORRIS_R");
		if (morris_r_it != gsa_opt_map.end())
		{
			convert_ip(morris_r_it->second, morris_r);
		}
		auto morris_p_it = gsa_opt_map.find("MORRIS_P");
		if (morris_p_it != gsa_opt_map.end())
		{
			convert_ip(morris_p_it->second, morris_p);
		}

		gsa_method = new MorrisMethod(adj_par_name_vec, fixed_pars, lower_bnd, upper_bnd, log_trans_pars,
			morris_p, morris_r, run_manager_ptr, &base_partran_seq, pest_scenario.get_ctl_ordered_obs_names(), &file_manager);
		auto  obs_type_grps = gsa_method->process_obt_file(file_manager.open_ifile_ext("obt"), file_manager);
		file_manager.close_file("obt");
	}
	else if (method != gsa_opt_map.end() && method->second == "SOBOL")
	{
		int n_sample = 100;
		auto morris_r_it = gsa_opt_map.find("SOBOL_SAMPLES");
		if (morris_r_it != gsa_opt_map.end())
		{
			convert_ip(morris_r_it->second, n_sample);
		}

		gsa_method = new Sobol(adj_par_name_vec, fixed_pars, lower_bnd, upper_bnd, n_sample, 
			run_manager_ptr, &base_partran_seq, pest_scenario.get_ctl_ordered_obs_names(), &file_manager);
	}
	else
	{
		throw PestError("A valid method for computing the sensitivity must be specified in the *.gsa file");
	}


	// make model runs
	Parameters model_pars = base_partran_seq.ctl2model_cp(ctl_par);
	run_manager_ptr->reinitialize();
	gsa_method->assemble_runs();
	run_manager_ptr->run();

	gsa_method->calc_sen(model_run, file_manager.open_ofile_ext("srw"), file_manager.open_ofile_ext("msn"), file_manager.open_ofile_ext("orw"));
	file_manager.close_file("srw");
	file_manager.close_file("msn");
	file_manager.close_file("orw");
	delete run_manager_ptr;
	cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
	char buf[256];
	OperSys::gets_s(buf, sizeof(buf));
}
