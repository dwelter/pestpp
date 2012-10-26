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
#include <utility>
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

RunManagerYAM::RunManagerYAM(const ModelExecInfo &_model_exec_info, const string &_port, const string &stor_filename, ofstream &_f_rmr)
	: RunManagerAbstract(_model_exec_info, stor_filename), port(_port), f_rmr(_f_rmr)
{
	w_init();
	int status;
	struct addrinfo hints;
	struct addrinfo *servinfo;
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

unordered_multimap<int, YamModelRun>::iterator RunManagerYAM::get_active_run_id(int socket)
{
	auto i = active_runs.begin();
	auto end = active_runs.end();

	for(; i!=end && i->second.get_socket() != socket; ++i)
	{}
	return i;
}

void RunManagerYAM::allocate_memory(const Parameters &model_pars, const Observations &obs)
{
	file_stor.reset(model_pars, obs);
    obs_name_vec = obs.get_keys();
	cur_group_id = NetPackage::get_new_group_id();
}

void  RunManagerYAM::free_memory()
{
	waiting_runs.clear();
	completed_runs.clear();
	zombie_runs.clear();
	failed_runs.clear();
	file_stor.free_memory();
}

int RunManagerYAM::add_run(const Parameters &model_pars)
{
	int run_id = file_stor.add_run(model_pars);
	YamModelRun new_run(run_id);
	waiting_runs.push_back(new_run);
	return run_id;
}


void RunManagerYAM::run()
{
	stringstream message;
	NetPackage net_pack;

	cout << "    running model " << waiting_runs.size() << " times" << endl;
	f_rmr << "running model " << waiting_runs.size() << " times" << endl;
	if(master.fd_count < 2) // first entry is the listener, slave apper after this
	{
		cout << endl << "      waiting for slaves to appear..." << endl << endl;
		f_rmr << endl << "    waiting for slaves to appear..." << endl << endl;
	}
	while (!active_runs.empty() || !waiting_runs.empty())
	{
		//schedule runs on available nodes
		schedule_runs();
		// get and process incomming messages
		listen();
	}
	total_runs += completed_runs.size();
	message.str("");
	message << "    " << completed_runs.size() << " runs complete";
	std::cout << message.str() << endl;
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
				}
			}
			else  // handle data from a client
			{
				process_message(i);				
			} // END handle data from client
		} // END got new incoming connection
	} // END looping through file descriptors
}


void RunManagerYAM::schedule_runs()
{
	NetPackage net_pack;
	for (auto it_run=waiting_runs.begin(); !slave_fd.empty() &&  it_run!=waiting_runs.end();)
	{
		bool success = schedule_run(it_run->get_id());
		if (success)
		{
			it_run = waiting_runs.erase(it_run);
		}
		else
		{
			++it_run;
		}
	}
}

void RunManagerYAM::process_message(int i_sock)
{
	NetPackage net_pack;
	int err;
	vector<string> sock_name = w_getnameinfo_vec(i_sock);
	if(( err=net_pack.recv(i_sock)) <=0) // error or lost connection
	{
		if (err < 0) {
			f_rmr << "receive failed from slave: " << sock_name[0] <<":" <<sock_name[1] << " - terminating slave" << endl;
			cerr << "receive failed from slave: " << sock_name[0] <<":" <<sock_name[1] << " - terminating slave" << endl;
			w_close(i_sock); // bye!
			FD_CLR(i_sock, &master); // remove from master set
		}
		else {
			f_rmr << "lost connection to slave: " << sock_name[0] <<":" <<sock_name[1] << endl;
			cerr << "lost connection to slave: " << sock_name[0] <<":" <<sock_name[1] << endl;
			w_close(i_sock); // bye!
			FD_CLR(i_sock, &master); // remove from master set
		}
		// remove run from active queue and return it to the waiting queue
		auto it_active = get_active_run_id(i_sock);
		if (it_active != active_runs.end())
		{
			YamModelRun &cur_run = it_active->second;
			if (completed_runs.find(cur_run.get_id()) == completed_runs.end()) //check if run has already finish on another node
			{
				waiting_runs.push_front(it_active->second);
			}
			active_runs.erase(it_active);
		}
	}
	else if (net_pack.get_type() == net_pack.RUN_DIR)
	{
		string work_dir(net_pack.get_data().data(), net_pack.get_data().size());
		f_rmr << "New Slave connection from: " << sock_name[0] <<":" <<sock_name[1] << "   working dir: " << work_dir << endl;
		cout << endl << "New Slave connection from: " << sock_name[0] <<":" <<sock_name[1] << "   working dir: " << work_dir << endl << endl;
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
		err = net_pack.send(i_sock, &data[0], data.size());
	}
	else if (net_pack.get_type() == net_pack.READY)
	{
		// ready message received from slave and add slave to deque
		slave_fd.push_back(i_sock);
	}

	else if (net_pack.get_type() == net_pack.RUN_FINISH && net_pack.get_groud_id() != cur_group_id)
	{
		// this is an old run that did not finish on time
		// just ignore it
	}
	else if (net_pack.get_type() == net_pack.RUN_FINISH)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_groud_id();
		f_rmr << "Run received from: " << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		cout << "Run received from: " << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		process_model_run(i_sock, net_pack);

	}
	else if (net_pack.get_type() == net_pack.RUN_FAILED)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_groud_id();
		f_rmr << "Run failed on slave:" << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		cout << "Run failed on slave:" << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		failed_runs.insert(make_pair(run_id, i_sock));
		auto it = get_active_run_id(i_sock);
		waiting_runs.push_front(YamModelRun(run_id, i_sock));
		active_runs.erase(it);
	}
	else
	{
		f_rmr << "Received unsupported message from slave:" << endl;
		net_pack.print_header(f_rmr);
		//save results from model run
	}
}

 bool RunManagerYAM::process_model_run(int sock_id, NetPackage &net_pack)
 {
	bool use_run = false;
	int run_id =net_pack.get_run_id();
	YamModelRun model_run(run_id,  sock_id);
	//check if another instance of this model run has already completed 
	if (completed_runs.find(run_id) == completed_runs.end())
	{
		completed_runs.insert(pair<int, YamModelRun>(run_id,  model_run));
		Parameters pars;
		Observations obs;
		vector<Transformable *> tr_vec;
		tr_vec.push_back(&pars);
		tr_vec.push_back(&obs);
		Serialization::unserialize(net_pack.get_data(), tr_vec);
		file_stor.update_run(run_id, pars, obs);
		use_run = true;
	}
	auto range_pair = active_runs.equal_range(run_id);
	//remaining runs with this id are not needed so mark them as zombies
	for ( auto b=range_pair.first; b!=range_pair.second; ++b)
	{
		if ( (*b).second.get_socket() != model_run.get_socket())
		{
			zombie_runs.insert(*b);
		}
	}
	active_runs.erase(range_pair.first, range_pair.second);
	return use_run;
 }

 bool RunManagerYAM::schedule_run(int run_id)
 {
	bool scheduled = false;
	auto it_sock = slave_fd.end(); // iterator to current socket

	if (completed_runs.count(run_id) > 0)
	{
		// run already completed on different node.  Do nothing
	}
	else if (failed_runs.count(run_id) == 0)
	{
		 // schedule a run on a slave
		NetPackage net_pack(NetPackage::START_RUN, cur_group_id, run_id, "");
		it_sock = slave_fd.begin();
	}
	else if (failed_runs.count(run_id) > 0)
	{
		for(it_sock=slave_fd.begin(); it_sock!=slave_fd.end(); ++it_sock)
		{
			auto fail_iter_pair = failed_runs.equal_range(run_id);

			auto i = fail_iter_pair.first;
			for(i = fail_iter_pair.first;
				i!= fail_iter_pair.second && i->second != *it_sock;
				++i) {}
			if (i == fail_iter_pair.second)  // This is slave has not previously failed on this run
			{
				// This run has not previously failed on this slave
				// Schedule run on it_sock
				break;
			}
		}
	}
	if (it_sock != slave_fd.end())
	{
		YamModelRun tmp_run(run_id, *it_sock);
		vector<char> data = file_stor.get_serial_pars(run_id);
		vector<string> sock_name = w_getnameinfo_vec(*it_sock);
		f_rmr << "Sending run to: " << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << cur_group_id << ", run id = " << run_id << ")" << endl;
		NetPackage net_pack(NetPackage::START_RUN, cur_group_id, run_id, "");
		int err = net_pack.send(*it_sock, &data[0], data.size());
		active_runs.insert(pair<int, YamModelRun>(run_id, tmp_run));
		slave_fd.erase(it_sock);
		scheduled = true;
	}
	return scheduled;
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
