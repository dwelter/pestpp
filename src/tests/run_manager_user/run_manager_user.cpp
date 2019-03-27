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


	//Configure PANTHER with file transfer info from pest_scenario
	auto transfer_files = pest_scenario.get_transferfile_vec();
	auto security_key = pest_scenario.get_security_key();
	auto security_method = RunManagerPanther::SecurityMethod::HMAC; //should get this from pest_scenario 
	run_manager_ptr->set_transfer_file_names(transfer_files);
	run_manager_ptr->set_transfer_security(security_method, security_key);


	//Add a run
	//Start the manager, allow slaves to join and complete the run.
	std::cout << endl;
	std::cout << "HHHHHHHHHHHHHHHHHHHHHHH Doing a run." << endl << endl;
	int first_run_id = run_manager_ptr->add_run(ctl_par, 1);
	run_manager_ptr->run();


	//Test 
	std::cout << endl;
	std::cout << "HHHHHHHHHHHHHHHHHHHHHHH Doing an empty run." << endl << endl;
	run_manager_ptr->run();


	if (true)
	{
		//Now we want to check whether the first run is still available.
		bool available = run_manager_ptr->is_run_last(first_run_id);
		if (available)
		{
			std::cout << endl;
			std::cout << "HHHHHHHHHHHHHHHHHHHHHHH Doing an empty run to retrieve a file from the slave that did run1." << endl << endl;

			//We still have the first run available.
			//Instruct the run manager to retrieve a file from the slave who performed that run. 
			//Call run() to complete the transfer.
			int file_to_retrieve_index_on_manager = 0;
			int file_to_retrieve_index_on_worker = 0;
			run_manager_ptr->transfer_file_from_worker(file_to_retrieve_index_on_worker, file_to_retrieve_index_on_manager, first_run_id); //this could be called multiple times for different files
			run_manager_ptr->run();
		}
	}


	if (true)
	{
		std::cout << endl;
		std::cout << "HHHHHHHHHHHHHHHHHHHHHHH Doing an empty run to send a file to the slave(s)." << endl << endl;

		//Send a file to all slaves.
		int file_to_send_index_on_manager = 1;
		int file_to_send_index_on_worker = 1;
		run_manager_ptr->transfer_file_to_all_workers(file_to_send_index_on_worker, file_to_send_index_on_manager); //this could be called multiple times for different files
		run_manager_ptr->run();
	}


	//Do another run runs.
	std::cout << endl;
	std::cout << "HHHHHHHHHHHHHHHHHHHHHHH Doing a run." << endl << endl;
	int second_run_id = run_manager_ptr->add_run(ctl_par, 1);
	run_manager_ptr->run();


	//End
	cout << endl << endl;
	system("PAUSE");
	return 0;
}