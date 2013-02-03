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
#include "Pest.h"
#include "Jacobian_1to1.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include <lapackpp.h>
#include "ModelRunPP.h"
#include "SVDASolver.h"
#include  "QSqrtMatrix.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "RunManagerGenie.h"
#include "RunManagerSerial.h"
#include "SVD_PROPACK.h"
#include "OutputFileWriter.h"
#include "YamrSlave.h"
#include "Serialization.h"

using namespace std;
using namespace pest_utils;

//using namespace pest_utils;

int main(int argc, char* argv[])
{
	string version = "1.2.0";
	string complete_path;
	enum class RunManagerType {SERIAL, YAMR, GENIE};
	string socket_str;

	if (argc >=2) {
		complete_path = argv[1];
	}
	else {
		cerr << endl<< "PEST++ Version " << version << endl << endl;
		cerr << "Usage: pest++ pest_ctl_file" << endl << endl;
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



	string filename = get_filename(complete_path);
	filename = remove_file_ext(filename); // remove .pst extension
	string pathname = get_pathname(complete_path);
	if (pathname.empty()) pathname = ".";
	FileManager file_manager(filename, pathname);

	ofstream &fout_rec = file_manager.rec_ofstream();
	cout << "PEST++ Version " << version << endl << endl;
	fout_rec << "PEST++ Version " << version << endl << endl;

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

	const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();
	
	ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
	Jacobian *base_jacobian_ptr = new Jacobian_1to1(file_manager);

	SVDSolver base_svd(&pest_scenario.get_control_info(), pest_scenario.get_svd_info(), &pest_scenario.get_base_group_info(), &pest_scenario.get_ctl_parameter_info(),
		&pest_scenario.get_ctl_observation_info(), file_manager, &pest_scenario.get_ctl_observations(), &obj_func,
		base_trans_seq, pest_scenario.get_prior_info_ptr(), *base_jacobian_ptr, pest_scenario.get_regul_scheme_ptr());

	base_svd.set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
	//Build Super-Parameter problem
	Jacobian *super_jacobian_ptr = new Jacobian(file_manager);
	ParamTransformSeq trans_svda;
	// method be involked as pointer as the transformation sequence it is added to will
	// take responsibility for destroying it
	TranSVD *tran_svd = new TranSVD("SVD Super Parameter Tranformation");
	trans_svda = base_trans_seq;
	trans_svda.push_back_derivative2numeric(tran_svd);

	ParameterGroupInfo sup_group_info;
	//
	ControlInfo svd_control_info = pest_scenario.get_control_info();
	svd_control_info.relparmax = 0.1;
	// Start Solution iterations
	cout << endl << endl;
	int n_base_iter = pest_scenario.get_pestpp_options().get_n_iter_base();
	int n_super_iter = pest_scenario.get_pestpp_options().get_n_iter_super();
	int super_nmax = pest_scenario.get_pestpp_options().get_super_nmax();
	double super_eigthres = pest_scenario.get_pestpp_options().get_super_eigthres();
	TerminationController termination_ctl(pest_scenario.get_control_info().noptmax, pest_scenario.get_control_info().phiredstp,
			pest_scenario.get_control_info().nphistp, pest_scenario.get_control_info().nphinored, pest_scenario.get_control_info().relparstp,
			pest_scenario.get_control_info().nrelpar);

	Parameters cur_ctl_parameters = pest_scenario.get_ctl_parameters();
	//strat junk
	//char *serial = cur_ctl_parameters.serialize();
	//Parameters junk;
	//junk.unserialize(serial);
	//end junk
	//Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
	//Niether of these will change over the course of the simulation
	run_manager_ptr->initialize(base_trans_seq.ctl2model_cp(cur_ctl_parameters), pest_scenario.get_ctl_observations());


	ModelRun optimum_run(&obj_func, base_trans_seq, pest_scenario.get_ctl_observations());
	// if noptmax=0 make one run with the intital parameters
	if (pest_scenario.get_control_info().noptmax == 0) {
		Parameters init_model_pars = base_trans_seq.ctl2model_cp(cur_ctl_parameters);
		optimum_run.set_ctl_parameters(init_model_pars);
		run_manager_ptr->reinitialize();
		run_manager_ptr->add_run(init_model_pars);
		run_manager_ptr->run();

        Parameters tmp_pars;
        Observations tmp_obs;
        run_manager_ptr->get_run(0, tmp_pars, tmp_obs);
		optimum_run.update(tmp_pars, tmp_obs, ModelRun::FORCE_PAR_UPDATE);
		// save parameters to .par file
		OutputFileWriter::write_par(file_manager.open_ofile_ext("par"), optimum_run.get_ctl_pars(), *(optimum_run.get_par_tran().get_offset_ptr()), 
			*(optimum_run.get_par_tran().get_scale_ptr()));
		file_manager.close_file("par");
		// save new residuals to .rei file

		OutputFileWriter::write_rei(file_manager.open_ofile_ext("rei"), 0, 
			*(optimum_run.get_obj_func_ptr()->get_obs_ptr()), 
			optimum_run.get_obs(), *(optimum_run.get_obj_func_ptr()),
			optimum_run.get_ctl_pars());
		file_manager.close_file("rei");
		run_manager_ptr->free_memory();
	}

	//Define model Run for Base Parameters (uses base parameter tranformations)
	for (int i_iter = 0; i_iter<pest_scenario.get_control_info().noptmax; ++i_iter)
	{
		ModelRun *cur_run;
		cur_run = &(base_svd.solve(*run_manager_ptr, termination_ctl, n_base_iter, cur_ctl_parameters, optimum_run));
		cur_ctl_parameters = base_svd.cur_model_run().get_ctl_pars();
		if(termination_ctl.check_last_iteration()) break;
		// Build Super Parameter or SVDA problem

	try
		{
			const vector<string> &nonregul_obs = pest_scenario.get_nonregul_obs();
			const vector<string> &pars = base_svd.cur_model_run().get_numeric_pars().get_keys();
			QSqrtMatrix Q_sqrt(pest_scenario.get_ctl_observation_info(), nonregul_obs, &pest_scenario.get_prior_info(), 1.0);
			(*tran_svd).update(*base_jacobian_ptr, Q_sqrt, base_svd.cur_model_run().get_numeric_pars(), super_nmax, super_eigthres, pars, nonregul_obs);
			sup_group_info = (*tran_svd).build_par_group_info(pest_scenario.get_base_group_info());		
			SVDASolver super_svd(&svd_control_info, pest_scenario.get_svd_info(), &sup_group_info, &pest_scenario.get_ctl_parameter_info(),
				&pest_scenario.get_ctl_observation_info(),  file_manager, &pest_scenario.get_ctl_observations(), &obj_func,
				trans_svda, &pest_scenario.get_prior_info(), *super_jacobian_ptr, pest_scenario.get_regul_scheme_ptr());
			super_svd.set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
			cur_run = &(super_svd.solve(*run_manager_ptr, termination_ctl, n_super_iter, cur_ctl_parameters, optimum_run));
			cur_ctl_parameters = super_svd.cur_model_run().get_ctl_pars();
		}
		catch(...)
		{
			cout << "WARNING: super parameter run failed.  Switching to base parameters" << endl;
			fout_rec << "WARNING: super parameter run failed.  Switching to base parameters" << endl;
		}
	}
	cout << endl;
	cout << "FINAL OPTIMISATION RESULTS" << endl << endl;
	fout_rec << "FINAL OPTIMISATION RESULTS" << endl << endl;
	optimum_run.full_report(cout);
	optimum_run.full_report(fout_rec);
	fout_rec.close();
	// clean up
	delete base_jacobian_ptr;
	delete run_manager_ptr;
	cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
	char buf[256];
    gets_s(buf, sizeof(buf));
}
