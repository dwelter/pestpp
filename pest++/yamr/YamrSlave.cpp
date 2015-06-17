#include "YamrSlave.h"
#include "utilities.h"
#include "Serialization.h"
#include "system_variables.h"
#include <cassert>
#include <cstring>
#include <algorithm>
#include <thread>
#include "system_variables.h"
#include "iopp.h"
#include "utilities.h"

using namespace pest_utils;

int  linpack_wrap(void);

void YAMRSlave::init_network(const string &host, const string &port)
{
	w_init();


	int status;
	struct addrinfo hints;
	struct addrinfo *servinfo;

	memset(&hints, 0, sizeof hints);
	//Use this for IPv4 aand IPv6
	//hints.ai_family = AF_UNSPEC;
	//Use this just for IPv4;
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;

	status = w_getaddrinfo(host.c_str(), port.c_str(), &hints, &servinfo);
	w_print_servinfo(servinfo, cout);
	cout << endl;
	// connect
	addrinfo* connect_addr = nullptr;
	while  (connect_addr == nullptr)
	{
		connect_addr = w_connect_first_avl(servinfo, sockfd);
		if (connect_addr == nullptr) {
			cerr << endl;
			cerr << "failed to connect to master" << endl;
		}
	}
	cout << "connection to master succeeded on socket: " << w_get_addrinfo_string(connect_addr) << endl << endl;
	freeaddrinfo(servinfo);

	fdmax = sockfd;
	FD_ZERO(&master);
	FD_SET(sockfd, &master);
	// send run directory to master
}


YAMRSlave::~YAMRSlave()
{
	w_close(sockfd);
	w_cleanup();
}


void YAMRSlave::process_ctl_file(const string &ctl_filename)
{
	ifstream fin;
	long lnum;
	long sec_begin_lnum;
	long sec_lnum;
	string section("");
	string line;
	string line_upper;
	vector<string> tokens;

	int num_par;
	int num_tpl_file;

	comline_vec.clear();
	tplfile_vec.clear();
	inpfile_vec.clear();
	insfile_vec.clear();
	outfile_vec.clear();
	fin.open(ctl_filename);
	try {
		for (lnum = 1, sec_begin_lnum = 1; getline(fin, line); ++lnum)
		{
			strip_ip(line);
			line_upper = upper_cp(line);
			tokens.clear();
			tokenize(line_upper, tokens);
			sec_lnum = lnum - sec_begin_lnum;
			if (tokens.empty())
			{
				//skip blank line
			}
			else if (line_upper.substr(0, 2) == "++")
			{
			}

			else if (line_upper[0] == '*')
			{
				section = upper_cp(strip_cp(line_upper, "both", " *\t\n"));
				sec_begin_lnum = lnum;
			}
			else if (section == "CONTROL DATA")
			{
				if (sec_lnum == 2)
				{
					convert_ip(tokens[0], num_par);
				}
				if (sec_lnum == 3)
				{
					convert_ip(tokens[0], num_tpl_file);
				}
			}
			else if (section == "MODEL COMMAND LINE")
			{
				comline_vec.push_back(line);
			}
			else if (section == "MODEL INPUT/OUTPUT")
			{
				vector<string> tokens_case_sen;
				tokenize(line, tokens_case_sen);
				if (sec_lnum <= num_tpl_file)
				{
					tplfile_vec.push_back(tokens_case_sen[0]);
					inpfile_vec.push_back(tokens_case_sen[1]);
				}
				else
				{
					insfile_vec.push_back(tokens_case_sen[0]);
					outfile_vec.push_back(tokens_case_sen[1]);
				}
			}
		}
	}
	catch (PestConversionError &e) {
		std::stringstream out;
		out << "Error parsing \"" << ctl_filename << "\" on line number " << lnum << endl;
		out << e.what() << endl;
		e.add_front(out.str());
		e.raise();
	}
	fin.close();
}



void YAMRSlave::process_yamr_ctl_file(const string &ctl_filename)
{
	ifstream fin;
	long lnum;
	long sec_begin_lnum;
	long sec_lnum;
	string section("");
	string line;
	string line_upper;
	vector<string> tokens;

	int num_par;
	int num_tpl_file;

	comline_vec.clear();
	tplfile_vec.clear();
	inpfile_vec.clear();
	insfile_vec.clear();
	outfile_vec.clear();
	fin.open(ctl_filename);
	try {
		for (lnum = 1, sec_begin_lnum = 1; getline(fin, line); ++lnum)
		{
			strip_ip(line);
			line_upper = upper_cp(line);
			tokens.clear();
			tokenize(line_upper, tokens);
			sec_lnum = lnum - sec_begin_lnum;
			if (tokens.empty())
			{
				//skip blank line
			}
			else if (line_upper.substr(0, 2) == "++")
			{
			}

			else if (line_upper[0] == '*')
			{
				section = upper_cp(strip_cp(line_upper, "both", " *\t\n"));
				sec_begin_lnum = lnum;
			}
			else if (section == "MODEL COMMAND LINE")
			{
				comline_vec.push_back(line);
			}
			else if (section == "MODEL INPUT")
			{
				vector<string> tokens_case_sen;
				tokenize(line, tokens_case_sen);
				tplfile_vec.push_back(tokens_case_sen[0]);
				inpfile_vec.push_back(tokens_case_sen[1]);
			}
			else if (section == "MODEL OUTPUT")
			{
				vector<string> tokens_case_sen;
				tokenize(line, tokens_case_sen);
				insfile_vec.push_back(tokens_case_sen[0]);
				outfile_vec.push_back(tokens_case_sen[1]);
			}
		}
	}
	catch (PestConversionError &e) {
		std::stringstream out;
		out << "Error parsing \"" << ctl_filename << "\" on line number " << lnum << endl;
		out << e.what() << endl;
		e.add_front(out.str());
		e.raise();
	}
	fin.close();
}


int YAMRSlave::recv_message(NetPackage &net_pack, struct timeval *tv)
{
	fd_set read_fds;
	int err = -1;
	int recv_fails = 0;
	while (recv_fails < max_recv_fails && err != 1)
	{
		read_fds = master; // copy master
		int result = w_select(fdmax + 1, &read_fds, NULL, NULL, tv);
		if (result == -1)
		{
			cerr << "fatal network error while receiving messages. ERROR: select() failure";
			return -990;
		}
		if (result == 0)
		{
			// no messages available for reading
			if (tv == NULL)
			{
				cerr << "fatal network error while receiving messages. ERROR: blocking select() call failure";
				return -990;
			}
			else
			{
				return 2;
			}
		}
		for (int i = 0; i <= fdmax; i++) {
			if (FD_ISSET(i, &read_fds)) { // got message to read
				err = net_pack.recv(i); // error or lost connection
				if (err == -2) {
					vector<string> sock_name = w_getnameinfo_vec(i);
					cerr << "recieved corrupt message from to master: " << sock_name[0] << ":" << sock_name[1] << endl;
					w_close(i); // bye!
					FD_CLR(i, &master); // remove from master set
					err = -999;
					return err;
				}
				else if (err < 0) {
					recv_fails++;
					vector<string> sock_name = w_getnameinfo_vec(i);
					cerr << "receive from master failed: " << sock_name[0] << ":" << sock_name[1] << endl;
					err = -1;
				}
				else if(err == 0) {
					vector<string> sock_name = w_getnameinfo_vec(i);
					cerr << "lost connection to master: " << sock_name[0] << ":" << sock_name[1] << endl;
						w_close(i); // bye!
						FD_CLR(i, &master); // remove from master set
						err = -999;
						return err;
				}
				else
				{
					// received data sored in net_pack return to calling routine to process it
					err = 1;
					return err;
				}
			}
		}
	}
	cerr << "send to master failed " << max_send_fails << " times, exiting..." << endl;
	return err;
	// returns -1  receive error
	//         -990  error in call to select()
	//         -991  connection closed
	//          1  message recieved
	//          2  no message recieved
}

int YAMRSlave::recv_message(NetPackage &net_pack, long  timeout_seconds, long  timeout_microsecs)
{
	int err = -1;
	int result = 0;
	struct timeval tv;
	tv.tv_sec = timeout_seconds;
	tv.tv_usec = timeout_microsecs;
	err = recv_message(net_pack, &tv);
	return err;
}


int YAMRSlave::send_message(NetPackage &net_pack, const void *data, unsigned long data_len)
{
	int err;
	int n;

	for (err = -1, n = 0; err != 1 && n < max_send_fails; ++n)
	{
		err = net_pack.send(sockfd, data, data_len);
	}
	if (n >= max_send_fails)
	{
		cerr << "send to master failed " << max_send_fails << " times, giving..." << endl;
	}
	return err;
}


NetPackage::PackType YAMRSlave::run_model(Parameters &pars, Observations &obs, NetPackage &net_pack)
{
	NetPackage::PackType final_run_status = NetPackage::PackType::RUN_FAILED;
	bool done = false;
	int err = 0;
	
	thread_flag f_terminate(false);
	thread_flag f_finished(false);
	thread_exceptions shared_execptions;
	try 
	{
		vector<string> par_name_vec;
		vector<double> par_values;
		for (auto &i : pars)
		{
			par_name_vec.push_back(i.first);
			par_values.push_back(i.second);
		}
		
		vector<double> obs_vec;	
		thread run_thread(w_write_run_read, &f_terminate, &f_finished, &shared_execptions,
			&tplfile_vec, &outfile_vec, &insfile_vec,
			&inpfile_vec, &par_name_vec, &par_values,
			&obs_name_vec, &obs_vec, &comline_vec);
		pest_utils::thread_RAII raii(run_thread);
		while (true)
		{
			if (shared_execptions.size() > 0)
			{
				cout << "exception raised by run thread " << std::endl;
				//don't break here, need to check one last time for incoming messages
				done = true;
			}
			//check if the runner thread has finished
			if (f_finished.get())
			{
				cout << "received finished signal from run thread " << std::endl;
				//don't break here, need to check one last time for incoming messages
				done = true;
			}
			//this call includes a "sleep" for the timeout
			err = recv_message(net_pack, 0, 100000);
			if (err < 0)
			{
				f_terminate.set(true);
				exit(-1);
			}
			//timeout on recv
			else if (err == 2)
			{
			}
			else if (net_pack.get_type() == NetPackage::PackType::PING)
			{
				cout << "ping request recieved...";
				net_pack.reset(NetPackage::PackType::PING, 0, 0, "");
				const char* data = "\0";
				err = send_message(net_pack, &data, 0);
				if (err != 1)
				{
					f_terminate.set(true);
					exit(-1);
				}
				cout << "ping response sent" << endl;
			}
			else if (net_pack.get_type() == NetPackage::PackType::REQ_KILL)
			{
				cout << "received kill request signal from master" << endl;
				cout << "sending terminate signal to run thread" << endl;
				f_terminate.set(true);
				final_run_status = NetPackage::PackType::RUN_KILLED;
				break;
			}
			else if (net_pack.get_type() == NetPackage::PackType::TERMINATE)
			{
				cout << "received terminate signal from master" << endl;
				cout << "sending terminate signal to run thread" << endl;
				f_terminate.set(true);
				terminate = true;
				final_run_status = NetPackage::PackType::TERMINATE;
				break;
			}
			else
			{
				cerr << "Received unsupported message from master, only PING REQ_KILL or TERMINATE can be sent during model run" << endl;
				cerr << static_cast<int>(net_pack.get_type()) << endl;
				cerr << "something is wrong...exiting" << endl;
				f_terminate.set(true);
				final_run_status = NetPackage::PackType::TERMINATE;
				exit(-1);
			}
			if (done) break;
		}
		shared_execptions.rethrow();
		if (!f_terminate.get())
		{
			//update the parameter values
			pars.clear();
			for (int i = 0; i < par_name_vec.size(); ++i)
			{
				pars[par_name_vec[i]] = par_values[i];
			}
			// update observation values		
			obs.clear();
			for (int i = 0; i < obs_name_vec.size(); ++i)
			{
				obs[obs_name_vec[i]] = obs_vec[i];
			}
			final_run_status = NetPackage::PackType::RUN_FINISHED;
		}
	}
	catch(const std::exception& ex)
	{
		cerr << endl;
		cerr << "   " << ex.what() << endl;
		cerr << "   Aborting model run" << endl << endl;
		NetPackage::PackType::RUN_FAILED;
	}
	catch(...)
	{
 		cerr << "   Error running model" << endl;
		cerr << "   Aborting model run" << endl;
		NetPackage::PackType::RUN_FAILED;
	}
	return final_run_status;
}




void YAMRSlave::check_io()
{
	vector<string> inaccessible_files;
	for (auto &file : insfile_vec)
	if (!check_exist_in(file)) inaccessible_files.push_back(file);
	for (auto &file : outfile_vec)
	if (!check_exist_out(file)) inaccessible_files.push_back(file);
	for (auto &file : tplfile_vec)
	if (!check_exist_in(file)) inaccessible_files.push_back(file);
	for (auto &file : inpfile_vec)
	if (!check_exist_out(file)) inaccessible_files.push_back(file);

	if (inaccessible_files.size() != 0)
	{
		string missing;
		for (auto &file : inaccessible_files)
			missing += file + " , ";
		throw PestError("Could not access the following model interface files: " + missing);
	}
}

void YAMRSlave::check_par_obs()
{
	TemplateFiles templatefiles(false, false, tplfile_vec, inpfile_vec,par_name_vec);
	templatefiles.check_parameter_names();
	InstructionFiles instructionfiles(insfile_vec, outfile_vec);
	instructionfiles.check_obs_names(obs_name_vec);
}

void YAMRSlave::start(const string &host, const string &port)
{
	NetPackage net_pack;
	Observations obs;
	Parameters pars;
	vector<int8_t> serialized_data;
	int err;

	//class attribute - can be modified in run_model()
	terminate = false;	
	init_network(host, port);
	while (!terminate)
	{
		//get message from master
		err = recv_message(net_pack);
		if (err < 0)
		{
			terminate = true;
		}
		else if(net_pack.get_type() == NetPackage::PackType::REQ_RUNDIR)
		{
			// Send Master the local run directory.  This information is only used by the master
			// for reporting purposes
			net_pack.reset(NetPackage::PackType::RUNDIR, 0, 0,"");
			string cwd =  OperSys::getcwd();
			err = send_message(net_pack, cwd.c_str(), cwd.size());
			if (err != 1)
			{
				exit(-1);
			}
		}
		else if (net_pack.get_type() == NetPackage::PackType::PAR_NAMES)
		{
			//Don't check first8 bytes as these contain an interger which stores the size of the data.
			bool safe_data = NetPackage::check_string(net_pack.get_data(), 0, net_pack.get_data().size());
			if (!safe_data)
			{
				cerr << "recieved corrupt parameter name packet from master" << endl;
				cerr << "terminating execution ..." << endl << endl;
				net_pack.reset(NetPackage::PackType::CORRUPT_MESG, 0, 0, "");
				char data;
				int np_err = send_message(net_pack, &data, 0);
				exit(-1);
			}
			Serialization::unserialize(net_pack.get_data(), par_name_vec);
		}
		else if (net_pack.get_type() == NetPackage::PackType::OBS_NAMES)
		{
			//Don't check first8 bytes as these contain an interger which stores the size of the data.
			bool safe_data = NetPackage::check_string(net_pack.get_data(), 0, net_pack.get_data().size());
			if (!safe_data)
			{
				cerr << "recieved corrupt observation name packet from master" << endl;
				cerr << "terminating execution ..." << endl << endl;
				net_pack.reset(NetPackage::PackType::CORRUPT_MESG, 0, 0, "");
				char data;
				int np_err = send_message(net_pack, &data, 0);
				exit(-1);
			}
			Serialization::unserialize(net_pack.get_data(), obs_name_vec);
			cout << "checking model IO files...";
			try
			{
				check_io();
				//check_par_obs();
			}
			catch (exception &e)
			{
				cerr << e.what() << endl;
				net_pack.reset(NetPackage::PackType::IO_ERROR, 0, 0,"");
				string err(e.what());
				vector<char> data(err.begin(), err.end());
				data.push_back('\0');
				int np_err = send_message(net_pack, &data, data.size());
				exit(-1);
			}
			cout << "done" << endl;
		}
		else if(net_pack.get_type() == NetPackage::PackType::REQ_LINPACK)
		{
			linpack_wrap();
			net_pack.reset(NetPackage::PackType::LINPACK, 0, 0,"");
			char data;
			err = send_message(net_pack, &data, 0);
			if (err != 1)
			{
				exit(-1);
			}
		}
		else if(net_pack.get_type() == NetPackage::PackType::START_RUN)
		{
			Serialization::unserialize(net_pack.get_data(), pars, par_name_vec);
			// run model
			int group_id = net_pack.get_group_id();
			int run_id = net_pack.get_run_id();
			cout << "received parameters (group id = " << group_id << ", run id = " << run_id << ")" << endl;
			cout << "starting model run..." << endl;

			NetPackage::PackType final_run_status = run_model(pars, obs, net_pack);
			if (final_run_status == NetPackage::PackType::RUN_FINISHED)
			{
				//send model results back
				cout << "run complete" << endl;
				cout << "sending results to master (group id = " << group_id << ", run id = " << run_id << ")..." << endl;
				cout << "results sent" << endl << endl;
				serialized_data = Serialization::serialize(pars, par_name_vec, obs, obs_name_vec);
				net_pack.reset(NetPackage::PackType::RUN_FINISHED, group_id, run_id, "");
				err = send_message(net_pack, serialized_data.data(), serialized_data.size());
				if (err != 1)
				{
					exit(-1);
				}
			}
			else if (final_run_status == NetPackage::PackType::RUN_FAILED)
			{
				cout << "run failed" << endl;
				net_pack.reset(NetPackage::PackType::RUN_FAILED, group_id, run_id, "");
				char data;
				err = send_message(net_pack, &data, 0);
				if (err != 1)
				{
					exit(-1);
				}
			}
			else if (final_run_status == NetPackage::PackType::RUN_KILLED)
			{
				cout << "run killed" << endl;
				net_pack.reset(NetPackage::PackType::RUN_KILLED, group_id, run_id, "");
				char data;
				err = send_message(net_pack, &data, 0);
				if (err != 1)
				{
					exit(-1);
				}
			}
			else if (final_run_status == NetPackage::PackType::TERMINATE)
			{
				cout << "run preempted by termination requested" << endl;
				terminate = true;
			}

			if (!terminate)
			{
				// Send READY Message to master
				cout << "sending ready signal to master" << endl;
				net_pack.reset(NetPackage::PackType::READY, 0, 0, "");
				char data;
				err = send_message(net_pack, &data, 0);
				if (err != 1)
				{
					exit(-1);
				}
			}
		}
		else if (net_pack.get_type() == NetPackage::PackType::TERMINATE)
		{
			cout << "terminated requested" << endl;
			terminate = true;
		}
		else if (net_pack.get_type() == NetPackage::PackType::REQ_KILL)
		{
			cout << "received kill request from master. run already finished" << endl;
		}
		else if (net_pack.get_type() == NetPackage::PackType::PING)
		{
			cout << "ping request recieved...";
			net_pack.reset(NetPackage::PackType::PING, 0, 0, "");
			const char* data = "\0";
			err = send_message(net_pack, &data, 0);
			if (err != 1)
			{
					exit(-1);
			}
			cout << "ping response sent" << endl;
		}
		else 
		{
			cout << "received unsupported messaged type: " << int(net_pack.get_type()) << endl;
		}
		//w_sleep(100);
		this_thread::sleep_for(chrono::milliseconds(100));
	}
}

