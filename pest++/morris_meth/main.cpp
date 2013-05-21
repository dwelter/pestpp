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
#include "MorrisMethod.h"
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
	string version = "1.0.0";
	string complete_path;
	enum class RunManagerType {SERIAL, YAMR, GENIE};
	string socket_str;

	if (argc >=2) {
		complete_path = argv[1];
	}
	else {
		cerr << endl<< "Morris++ Version " << version << endl << endl;
		cerr << "Usage: morris pest_ctl_file" << endl << endl;
		exit(1);
	}


	// This is a YAMR Slave, start PEST++ as a YAMR Slave
	if (argc >=3 && upper(argv[1]) == "/H") {
		string socket_str = argv[2];
		strip_ip(socket_str);
		vector<string> sock_parts;
		tokenize(socket_str, sock_parts, ":");
		YAMRSlave yam_slave;
		yam_slave.start(sock_parts[0], sock_parts[1]);
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

	 int morris_r = 7;

	string filename = get_filename(complete_path);
	filename = remove_file_ext(filename); // remove .pst extension
	string pathname = get_pathname(complete_path);
	if (pathname.empty()) pathname = ".";
	FileManager file_manager(filename, pathname);

	ofstream &fout_rec = file_manager.rec_ofstream();
	cout << "Morris Version " << version << endl << endl;
	fout_rec << "Morris Version " << version << endl << endl;

	cout << "using control file: \"" <<  complete_path << "\"" << endl << endl;
	fout_rec << "Control file = " <<  complete_path  << "\"" << endl << endl;

	// create pest run and process control file to initialize it
	Pest pest_scenario;
	pest_scenario.set_defaults();
	try {
		pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager);
		file_manager.close_file("pst");
	}
	catch(PestFileError)
	{
		exit(1);
	}
	pest_scenario.check_inputs();

	RunManagerAbstract *run_manager_ptr;
	if (run_manager_type == RunManagerType::YAMR)
	{
		cout << "initializing YAMR run manager" << endl;
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
		cout << "initializing Genie run manager" << endl;
		string socket_str = argv[3];
		strip_ip(socket_str);
		const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
		run_manager_ptr = new RunManagerGenie(exi.comline_vec,
		exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
		file_manager.build_filename("rns"), socket_str);
	}
	else
	{
		cout << "initializing serial run manager" << endl;
		const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
		run_manager_ptr = new RunManagerSerial(exi.comline_vec,
		exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
		file_manager.build_filename("rns"), pathname);
	}

	// Get the lower bounds of the parameters
	Parameters ctl_par = pest_scenario.get_ctl_parameters();
	//remove fixed parameters from ctl_par
	vector<string> par_name_vec;
	Parameters fixed_pars;
	for (auto &i : pest_scenario.get_ctl_ordered_par_names())
	{
		if(pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(i)->tranform_type != ParameterRec::TRAN_TYPE::FIXED)
		{
			par_name_vec.push_back(i);
		}
		else
		{
			fixed_pars.insert(i, pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(i)->init_value);
		}
	}


	Parameters lower_bnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(ctl_par.get_keys());
	Parameters upper_bnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(ctl_par.get_keys());

	//Build Transformation with ctl_2_numberic
	ParamTransformSeq base_partran_seq(pest_scenario.get_base_par_tran_seq());
	ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
	ModelRun model_run(&obj_func, pest_scenario.get_ctl_observations());

	MorrisMethod morris(par_name_vec, fixed_pars, lower_bnd, upper_bnd, 8); //8 levels for each parameters

	// make model runs
	Parameters model_pars = base_partran_seq.ctl2model_cp(ctl_par);
	run_manager_ptr->reinitialize();
	for (int i=0; i<morris_r; ++i)
	{
		morris.assemble_runs(*run_manager_ptr, base_partran_seq);
	}
	run_manager_ptr->run();

	morris.calc_sen(*run_manager_ptr, model_run, file_manager.open_ofile_ext("msn"));
	file_manager.close_file("msn");
	delete run_manager_ptr;
	cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
	char buf[256];
	OperSys::gets_s(buf, sizeof(buf));
}
