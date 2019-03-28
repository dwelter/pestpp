#include "RunManagerPanther.h" //needs to be first because it includes winsock2.h
#include "pch.h"
#include <iostream>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <FileManager.h>
#include "utilities.h"
#include "Pest.h"

using namespace std;
using namespace pest_utils;



int main(int argc, char* argv[])
{
	char endl[] = "\n";
	std::cout << endl;
	std::cout << "'RunManagerUser' test project started." << endl << endl;


	//Get a FileManager. This is handy for working with PEST paths and files.
	FileManager file_manager;
	file_manager.initialize_path(get_filename_without_ext("control_file.pst"), ".");


	//Create pest run and process control file to initialize it
	Pest pest_scenario;
	pest_scenario.set_defaults();
	pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"));
	file_manager.close_file("pst");


	//Get the port from the socket_str
	string port = "3801";
	strip_ip(port);
	strip_ip(port, "front", ":");
	cout << "using socket: \"" << port << "\"" << endl;


	//Get a PANTHER run manager instance
	RunManagerPanther *run_manager_ptr;
	run_manager_ptr = new RunManagerPanther(
		file_manager.build_filename("rns"), port,
		file_manager.open_ofile_ext("rmr"),
		pest_scenario.get_pestpp_options().get_max_run_fail(),
		pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
		pest_scenario.get_pestpp_options().get_overdue_giveup_fac()); 


	//Build Transformation with ctl_2_numberic
	ParamTransformSeq base_partran_seq(pest_scenario.get_base_par_tran_seq());
	Parameters ctl_par = pest_scenario.get_ctl_parameters();
	run_manager_ptr->initialize(base_partran_seq.ctl2model_cp(ctl_par), pest_scenario.get_ctl_observations());


	//Configure PANTHER with file transfer settings info from pest_scenario
	run_manager_ptr->set_transfer_file_names(pest_scenario.get_transferfile_vec());
	run_manager_ptr->set_transfer_security(pest_scenario.get_security_method(), pest_scenario.get_security_key());


	//Do a single run.
	int first_run_id = run_manager_ptr->add_run(ctl_par, 1);
	run_manager_ptr->run();


	//If the run is still available request files 0 and 1 from that slave 
	//and save them as 4 and 5.
	if (run_manager_ptr->is_run_last(first_run_id))
	{
		run_manager_ptr->queue_file_transfer_from_worker(0, 4, first_run_id);
		run_manager_ptr->queue_file_transfer_from_worker(1, 5, first_run_id);
		run_manager_ptr->run();
	}


	//Send files 2 and 3 to all slaves and save them as 4 and 5.
	run_manager_ptr->queue_file_transfer_to_workers(2, 6);
	run_manager_ptr->queue_file_transfer_to_workers(3, 7);
	run_manager_ptr->run();


	//Do another run for fun.
	int second_run_id = run_manager_ptr->add_run(ctl_par, 1);
	run_manager_ptr->run();


	//End
	cout << endl << endl;
	system("PAUSE");
	return 0;
}