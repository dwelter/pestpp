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
#include <winsock2.h>
#include <ws2tcpip.h>
#include "RunManagerYAM.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>
#include <cstring>
#include <map>
#include <direct.h>
#include <deque>
#include "network_wrapper.h"
#include "network_package.h"
#include "Transformable.h"
#include "ModelRunPP.h"
#include "Pest.h"
#include "utilities.h"
#include "Serialization.h"


using namespace std;
using namespace pest_utils;

YamModelRun::YamModelRun(int _run_id, int _sockfd) : sockfd(_sockfd), run_id(_run_id)
{
}

RunManagerYAM::RunManagerYAM(const ModelExecInfo &_model_exec_info, const std::vector<std::string>& _obs_name_vec, const string &_port, const string &stor_filename, ofstream &_f_rmr)
	: RunManagerAbstract(_model_exec_info), port(_port), file_stor(stor_filename), f_rmr(_f_rmr)
{
	w_init();
	int status;
	struct addrinfo hints;
	struct addrinfo *servinfo;
	// clean this up 
	obs_name_vec = _obs_name_vec;
	memset(&hints, 0, sizeof hints);
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;

	status = w_getaddrinfo(NULL, port.c_str(), &hints, &servinfo);
	cout << "Starting YAM - Yet Another run Manager" << endl;
	w_print_servinfo(servinfo, cout);
	//make socket, bind and listen
	listener = w_socket(servinfo->ai_family, servinfo->ai_socktype, servinfo->ai_protocol);
	w_bind(listener, servinfo->ai_addr, servinfo->ai_addrlen);
	w_listen(listener, BACKLOG);
	//free servinfo
	freeaddrinfo(servinfo);
	fdmax = listener;
	FD_ZERO(&master);
	FD_SET(listener, &master);
	return;
}


void RunManagerYAM::allocate_memory(const Parameters &model_pars, const Observations &obs, int _max_runs)
{
	file_stor.reset(model_pars, obs);
}

void  RunManagerYAM::free_memory()
{
	waiting_runs.clear();
	completed_runs.clear();
	zombie_runs.clear();
	file_stor.free_memory();
}

int RunManagerYAM::add_run(const Parameters &model_pars)
{
	int run_id = file_stor.add_run(model_pars);
	YamModelRun new_run(run_id);
	waiting_runs.push_back(new_run);
	nruns++;
	return run_id;
}


void RunManagerYAM::run()
{
	int ifail;
	int success_runs;
	stringstream message;
	NetPackage net_pack;
	int group_id = net_pack.get_new_group_id();

	success_runs = 0;
	if(slave_fd.empty())
	{
		cout << endl << "      waiting for slaves to appear..." << endl;
		f_rmr << endl << "Waiting for slaves to appear..." << endl << endl;
	}
	while (!active_runs.empty() || !waiting_runs.empty())
	{
		//schedule runs on available nodes
		schedule_runs(group_id);
		// get and process incomming messages
		listen();
	}
	total_runs += success_runs;
	std::cout << string(message.str().size(), '\b');
	message.str("");
	message << "(" << success_runs << "/" << nruns << " runs complete)";
	std::cout << message.str();
	//if (success_runs < i_run)
	//{
	//	cout << endl << endl;
	//	cout << "WARNING: " << i_run-success_runs << " out of " <<i_run << " runs failed" << endl << endl;
	//}
}


void RunManagerYAM::listen()
{
	struct sockaddr_storage remote_addr;
	fd_set read_fds; // temp file descriptor list for select()
	int addr_len;
	timeval tv;
	tv.tv_sec = 1;
	tv.tv_usec = 0;

	read_fds = master; // copy it
	if (w_select(fdmax+1, &read_fds, NULL, NULL, &tv) == -1) 
	{
		return;
	}
	// run through the existing connections looking for data to read
	for(int i = 0; i <= fdmax; i++) {
		if (FD_ISSET(i, &read_fds)) { // we got one!!
			if (i == listener)  // handle new connections
			{
				int newfd;
				addr_len = sizeof remote_addr;
				newfd = w_accept(listener,(struct sockaddr *)&remote_addr, &addr_len);
				if (newfd == -1) {}
				else 
				{
					FD_SET(newfd, &master); // add to master set
					if (newfd > fdmax) { // keep track of the max
						fdmax = newfd;
					}
					vector<string> sock_name = w_getnameinfo_vec(newfd);
					//f_rmr << "New Slave connection from: " << sock_name[0] <<":" <<sock_name[1] << endl;
					//slave_fd.push_back(newfd);
				}
			}
			else  // handle data from a client
			{
				process_message(i);				
			} // END handle data from client
		} // END got new incoming connection
	} // END looping through file descriptors
}


void RunManagerYAM::schedule_runs(int group_id)
{
	NetPackage net_pack;

	while(!slave_fd.empty() && !waiting_runs.empty())
	{
			// schedule a run on a slave
			int err;
			int sock_fd = slave_fd.front();
			YamModelRun tmp_run =  waiting_runs.front();
			tmp_run.set_socket(sock_fd);
			int run_id = tmp_run.get_id();
			vector<char> data = file_stor.get_serial_pars(run_id);
			vector<string> sock_name = w_getnameinfo_vec(sock_fd);
			f_rmr << "Sending run to: " << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
			net_pack.reset(NetPackage::START_RUN, group_id, run_id, "");
			err = net_pack.send(sock_fd, &data[0], data.size());
			active_runs.insert(pair<int, YamModelRun>(run_id, tmp_run));
			waiting_runs.pop_front();
			slave_fd.pop_front();
	}
}

void RunManagerYAM::process_message(int i)
{
	NetPackage net_pack;
	int err;
	vector<string> sock_name = w_getnameinfo_vec(i);
	if(( err=net_pack.recv(i)) <=0) // error or lost connection
	{
		if (err < 0) {
			f_rmr << "receive failed from slave: " << sock_name[0] <<":" <<sock_name[1] << endl;
			w_close(i); // bye!
			FD_CLR(i, &master); // remove from master set
		}
		else {
			f_rmr << "lost connection to slave: " << sock_name[0] <<":" <<sock_name[1] << endl;
			w_close(i); // bye!
			FD_CLR(i, &master); // remove from master set
		}
	}
	else if (net_pack.get_type() == net_pack.RUN_DIR)
	{
		string work_dir(net_pack.get_data().data(), net_pack.get_data().size());
		f_rmr << "New Slave connection from: " << sock_name[0] <<":" <<sock_name[1] << "   working dir: " << work_dir << endl;
	}
	else if (net_pack.get_type() == net_pack.REQ_CMD)
	{
		net_pack.reset(NetPackage::CMD, 0, 0, "");
		vector<char> data;
		vector<vector<string> *> tmp_vec;
		tmp_vec.push_back(&comline_vec);
		tmp_vec.push_back(&tplfile_vec);
		tmp_vec.push_back(&inpfile_vec);
		tmp_vec.push_back(&insfile_vec);
		tmp_vec.push_back(&outfile_vec);
		tmp_vec.push_back(&obs_name_vec);

		data = Serialization::serialize(tmp_vec);
		err = net_pack.send(i, &data[0], data.size());
	}
	else if (net_pack.get_type() == net_pack.READY)
	{
		// ready message received from slave and add slave to deque
		slave_fd.push_back(i);
	}
	else if (net_pack.get_type() == net_pack.RUN_FINISH)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_groud_id();
		f_rmr << "Run received from: " << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		if (process_model_run(YamModelRun(run_id, i)) == true)
		{
			Parameters pars;
			Observations obs;
			vector<Transformable *> tr_vec;
			tr_vec.push_back(&pars);
			tr_vec.push_back(&obs);
			Serialization::unserialize(net_pack.get_data(), tr_vec);
			file_stor.update_run(run_id, pars, obs);
		}
	}
	else if (net_pack.get_type() == net_pack.RUN_FAILED)
	{
		f_rmr << "Run failed on slave:" << endl;
		net_pack.print_header(f_rmr);
	}
	else
	{
		f_rmr << "Received unsupported message from slave:" << endl;
		net_pack.print_header(f_rmr);
		//save results from model run
	}
}

void RunManagerYAM::get_run(ModelRun &model_run, int run_id, PAR_UPDATE update_type)
{
	Parameters pars;
	Observations obs;
	file_stor.get_run(run_id, &pars, &obs);

	//Must set parameters before observations
	if(update_type == FORCE_PAR_UPDATE || model_run.get_par_tran().is_one_to_one())
	{
		// transform to numeric parameters
		model_run.get_par_tran().model2numeric_ip(pars);
		model_run.set_numeric_parameters(pars);
	}

	// Process Observations
	model_run.set_observations(obs);
}

 Parameters RunManagerYAM::get_model_parameters(int run_num) const
 {
	 return Parameters();
 }

 bool RunManagerYAM::process_model_run(YamModelRun model_run)
 {
	 bool use_run = false;
	 pair<unordered_multimap<int, YamModelRun>::iterator,unordered_multimap<int, YamModelRun>::iterator> range_pair;
	int run_id = model_run.get_id();
	//check if another instance of this model run has already completed 
	if (completed_runs.find(run_id) == completed_runs.end())
	{
		completed_runs.insert(pair<int, YamModelRun>(run_id,  model_run));
		use_run = true;
	}
	range_pair = active_runs.equal_range(run_id);
	//remaining runs with this id are not needed so mark them as zombies
	for ( unordered_multimap<int, YamModelRun>::iterator b=range_pair.first; b!=range_pair.second; ++b)
	{
		if ( (*b).second.get_socket() != model_run.get_socket())
		{
			zombie_runs.insert(*b);
		}
	}
	active_runs.erase(range_pair.first, range_pair.second);
	return use_run;
 }

RunManagerYAM::~RunManagerYAM(void)
{
	//close sockets and cleanup
	for(int i = 0; i <= fdmax; i++) {
		if (FD_ISSET(i, &master)) { // we got one!!
			NetPackage netpack(NetPackage::TERMINATE, 0, 0,"");
			char data;
			netpack.send(i, &data, 0);
			w_close(i);
		}
	}
	w_cleanup();
}
