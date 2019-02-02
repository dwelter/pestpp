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


	vector<string> cmd_arg_vec(argc);
	copy(argv, argv + argc, cmd_arg_vec.begin());
	for (vector<string>::iterator it = cmd_arg_vec.begin(); it != cmd_arg_vec.end(); ++it)
	{
		transform(it->begin(), it->end(), it->begin(), ::tolower);
	}



	//Get the control file name.
	string complete_path;
	if (argc >= 2) {
		complete_path = argv[1];
		std::cout << "current path: " << std::experimental::filesystem::current_path() << endl;
		cout << "using control file: \"" << complete_path << "\"" << endl;
	}
	else {
		cerr << "--------------------------------------------------------" << endl;
		cerr << "usage:" << endl << endl;
		cerr << "    PANTHER master:" << endl;
		cerr << "        run_manager_user.exe control_file.pst /H :port" << endl << endl;
		cerr << "--------------------------------------------------------" << endl;
		exit(0);
	}


	//Get a FileManager. This is handy for working with PEST paths and files.
	FileManager file_manager;
	string filename = complete_path;
	string pathname = ".";
	file_manager.initialize_path(get_filename_without_ext(filename), pathname);


	//Get the socket from the commandline.
	vector<string>::const_iterator it_find;
	it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/h");
	string socket_str = "";
	string next_item = "";
	next_item = *(it_find + 1);
	strip_ip(next_item);
	socket_str = next_item;
	cout << "using socket: \"" << socket_str << "\"" << endl;


	//Create pest run and process control file to initialize it
	Pest pest_scenario;
	pest_scenario.set_defaults();
	pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"));
	file_manager.close_file("pst");


	//Get the port from the socket_str
	string port = socket_str;
	strip_ip(port);
	strip_ip(port, "front", ":");


	//Get a PANTHER run manager instance
	RunManagerPanther *run_manager_ptr;
	run_manager_ptr = new RunManagerPanther(
		file_manager.build_filename("rns"), port,
		file_manager.open_ofile_ext("rmr"),
		pest_scenario.get_pestpp_options().get_max_run_fail(),
		pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
		pest_scenario.get_pestpp_options().get_overdue_giveup_fac()); 



	//Somewhere here, the softare using the PANTHER manager may want to send or recieve files to the slaves.
	//Are these transfers assoicated with specific runs? Or specific slaves? 
	//
	//I expected to be able to start the run manager and have slaves connect and wait for my instruction 
	//without terminating. But the run manager just wants to get on with runs, so I'm not sure how the file 
	//transfer functionality should fit within the current PANTHER API/ethos. To demonstrate what I have done
	//I have written a new function inside RunManagerPanther, called RunManagerPanther::file_transfer_demonstration()
	//which can be armed by setting the following flag.
	run_manager_ptr->demo_file_transfer = true;								//This flag means that when it a slave connects, the manager will demonstrate file transfer.

	//Build Transformation with ctl_2_numberic
	ParamTransformSeq base_partran_seq(pest_scenario.get_base_par_tran_seq());
	Parameters ctl_par = pest_scenario.get_ctl_parameters();

	//Add a run and start the run manager.
	run_manager_ptr->initialize(base_partran_seq.ctl2model_cp(ctl_par), pest_scenario.get_ctl_observations());
	int run_id = run_manager_ptr->add_run(ctl_par, 1);
	//run_manager_ptr->run_until(RunManagerAbstract::RUN_UNTIL_COND::TIME, 0, 10);
	run_manager_ptr->run();



	cout << endl << endl;
	system("PAUSE");
	return 0;
}