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

//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "Pest.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include <lapackpp.h>
#include "ModelRunPP.h"
#include "SVDASolver.h"
#include "utilities.h"
#include  "QSqrtMatrix.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "RunManagerGenie.h"
#include "RunManagerSerial.h"
#include "SVD_PROPACK.h"

using namespace std;
using namespace pest_utils;

//using namespace pest_utils;

int main(int argc, char* argv[])
{
	//SVD_PROPACK pro_test;
	//pro_test.test();
	//cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
	//char bufb[256];
 //   gets(bufb);
	//exit(0);

	string version = "1.1.1";
	string complete_path;
	if (argc >=2) {
		complete_path = argv[1];
	}
	else {
		cerr << "Usage: pest++ pest_ctl_file" << endl;
		exit(1);
	}

	string filename = get_filename(complete_path);
	string pathname = get_pathname(complete_path);
	if (pathname.empty()) pathname = ".";
	FileManager file_manager(filename, pathname);

	ofstream &fout_rec = file_manager.rec_ofstream();
	cout << "PEST++ Version " << version << endl << endl;
	fout_rec << "PEST++ Version " << version << endl << endl;

	cout << "using control file: \"" <<  complete_path << "\"" << endl << endl;
	fout_rec << "Control file = " <<  complete_path  << "\"" << endl << endl;

	// create to pest run and process control file to initialize it initialize
	Pest pest_scenario;
	pest_scenario.set_defaults();
	pest_scenario.process_ctl_file(file_manager.ctl_filename(), file_manager);
	pest_scenario.check_inputs();

	RunManagerAbstract *run_manager_ptr;
	if (pest_scenario.get_pestpp_options().get_gman_socket().empty())
	{
		cout << "initializing serial run manager" << endl;
		run_manager_ptr = new RunManagerSerial(pest_scenario.get_model_exec_info(), pathname);
	}
	else {
		cout << "initializing Genie run manager" << endl;
		run_manager_ptr = new RunManagerGenie (pest_scenario.get_model_exec_info(),  pest_scenario.get_pestpp_options().get_gman_socket());
	}

	const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();
	
	ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
	Jacobian *base_jacobian_ptr = new Jacobian(file_manager);

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
	trans_svda.push_back_ctl2numeric(tran_svd);

	ParameterGroupInfo sup_group_info;
	//
	ControlInfo svd_control_info = pest_scenario.get_control_info();
	svd_control_info.relparmax = 0.2;
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
	//Define model Run for Base Parameters (uses base parameter tranformations)
	ModelRun optimum_run(&obj_func, base_trans_seq, pest_scenario.get_ctl_observations());
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
			(*tran_svd).update(*base_jacobian_ptr, Q_sqrt, super_nmax, super_eigthres, pars, nonregul_obs);
			sup_group_info = (*tran_svd).build_par_group_info(pest_scenario.get_base_group_info());		
			SVDASolver super_svd(&svd_control_info, pest_scenario.get_svd_info(), &sup_group_info, &pest_scenario.get_ctl_parameter_info(),
				&pest_scenario.get_ctl_observation_info(),  file_manager, &pest_scenario.get_ctl_observations(), &obj_func,
				trans_svda, &pest_scenario.get_prior_info(), *super_jacobian_ptr, pest_scenario.get_regul_scheme_ptr());
			super_svd.set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
			cout << endl;
			cur_run = &(super_svd.solve(*run_manager_ptr, termination_ctl, n_super_iter, cur_ctl_parameters, optimum_run));
			cur_ctl_parameters = super_svd.cur_model_run().get_ctl_pars();
		}
		catch(...)
		{
			cout << "WARNING: super parameter run failed.  Switching to base parameters" << endl;
		}
	}
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
