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

#include <iostream>
#include <fstream>
#include "MorrisMethod.h"
#include <Eigen/Dense>
#include <utilities.h>
#include <FileManager.h>
#include <Pest.h>
#include <RunManagerGenie.h>
#include <RunManagerSerial.h>


using namespace std;
using namespace pest_utils;
using Eigen::MatrixXd;



int main(int argc, char* argv[])
{
	cout << "Starting Program" << endl;
		string version = "1.0.1";
	string complete_path;
	if (argc >=2) {
		complete_path = argv[1];
	}
	else {
		cerr << "Usage: gsen pest_ctl_file" << endl;
		exit(1);
	}

	string filename = get_filename(complete_path);
	string pathname = get_pathname(complete_path);
	if (pathname.empty()) pathname = ".";
	FileManager file_manager(filename, pathname);

	ofstream &fout_rec = file_manager.rec_ofstream();
	cout << "gsen Version " << version << endl << endl;
	fout_rec << "gsen Version " << version << endl << endl;

	cout << "using control file: \"" <<  complete_path << "\"" << endl << endl;
	fout_rec << "Control file = " <<  complete_path  << "\"" << endl << endl;

	// create to pest run and process control file to initialize it initialize
	Pest pest_scenario;
	pest_scenario.set_defaults();
	pest_scenario.process_ctl_file(file_manager.ctl_filename(), file_manager);
	pest_scenario.check_inputs();

	// Get the lower bounds of the parameters
	Parameters ctl_par = pest_scenario.get_ctl_parameters();
	Parameters lower_bnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(ctl_par.get_keys());
	Parameters upper_bnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(ctl_par.get_keys());
	Parameters delta;
	//delta = upper_bnd - lower_bnd;
	cout << lower_bnd << endl;
	//Build Transformation with ctl_2_numberic
	ParamTransformSeq base_partran_seq(pest_scenario.get_base_par_tran_seq());
	

	RunManagerAbstract *run_manager_ptr;
	if (pest_scenario.get_pestpp_options().get_gman_socket().empty())
	{
		cout << "initializing serial run manager" << endl;
		run_manager_ptr = new RunManagerSerial(pest_scenario.get_model_exec_info(), pathname);
	}
	//else {
	//	cout << "initializing Genie run manager" << endl;
	//	run_manager_ptr = new RunManagerGenie (pest_scenario.get_model_exec_info(),  pest_scenario.get_pestpp_options().get_gman_socket());
	//}

	MatrixXd b_star_mat;
	MorrisMethod morris(8); //8 levels for each parameters

	b_star_mat = morris.create_P_star_mat(ctl_par.size()); 
	cout << b_star_mat << endl << endl;
	// make model runs
	ParamTransformSeq ctl2model_tran = pest_scenario.get_base_par_tran_seq();
	Parameters model_pars = ctl2model_tran.ctl2model_cp(ctl_par);
	//run_manager_ptr->allocate_memory(model_pars, pest_scenario.get_ctl_observations(), b_star_mat.rows());
	//run_manager_ptr->add_run(ctl2model_tran.ctl2model_cp(ctl_par));


	cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
	char buf[256];
    gets_s(buf, sizeof(buf));
}