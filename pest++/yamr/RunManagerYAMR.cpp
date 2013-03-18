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

#include "network_wrapper.h"
#include "RunManagerYAMR.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>
#include <cstring>
#include <map>
#include <deque>
#include <utility>
#include <algorithm>
#include "network_wrapper.h"
#include "network_package.h"
#include "Transformable.h"
#include "utilities.h"
#include "Serialization.h"


using namespace std;
using namespace pest_utils;

YamrModelRun::YamrModelRun(int _run_id, int _sockfd) : sockfd(_sockfd), run_id(_run_id)
{
}

SlaveInfo::SlaveRec::SlaveRec()
{
	state = SlaveInfo::State::NEW;
	work_dir = "";
	linpack_time = std::chrono::hours(-500);
	run_time = std::chrono::hours(-500);
	start_time = std::chrono::system_clock::now();
}

bool SlaveInfo::CompareTimes::operator() (int a, int b)
{
	bool ret = false;
	auto &info_map = my_class_ptr->slave_info_map;
	auto ita = info_map.find(a);
	assert(ita != info_map.end());
	auto itb = info_map.find(b);
	assert(itb != info_map.end());
	if (itb->second.run_time > std::chrono::milliseconds(0) && ita->second.run_time > std::chrono::milliseconds(0))
	{
		bool ret = ! (itb->second.run_time > ita->second.run_time);
	}
	else if (itb->second.linpack_time > std::chrono::milliseconds(0) && ita->second.linpack_time > std::chrono::milliseconds(0))
	{
		bool ret = ! (itb->second.linpack_time > ita->second.linpack_time);
	}
	return ret;
}

SlaveInfo::SlaveInfo()
{
}

SlaveInfo::~SlaveInfo()
{
}

void SlaveInfo::add(int sock_id)
{
	slave_info_map[sock_id] = SlaveInfo::SlaveRec();
}

void SlaveInfo::erase(int sock_id)
{
	slave_info_map.erase(sock_id);
}

size_t SlaveInfo::size() const
{
	return slave_info_map.size();
}

SlaveInfo::State SlaveInfo::get_state(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	return it->second.state;
}

 void SlaveInfo::set_state(int sock_id, const State _state)
 {
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	it->second.state = _state;
 }

void SlaveInfo::set_work_dir(int sock_id, const std::string &work_dir)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	it->second.work_dir = work_dir;
}

string SlaveInfo::get_work_dir(int sock_id) const
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	string work_dir = it->second.work_dir;
	return work_dir;
}

void SlaveInfo::start_timer(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	it->second.start_time = std::chrono::system_clock::now();
}

void SlaveInfo::end_run(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	auto &run_time = it->second.run_time;
	auto dt = std::chrono::system_clock::now() - it->second.start_time;
	if (run_time > std::chrono::hours(0))
	{
		run_time = run_time + dt;
		run_time /= 2;
	}
	else
	{
		run_time = dt;
	}
}

void SlaveInfo::end_linpack(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	auto &linpack_time = it->second.linpack_time;
	linpack_time = std::chrono::system_clock::now() - it->second.start_time;
}


double SlaveInfo::get_runtime(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	auto &run_time = it->second.run_time;
	return double(run_time.count());
}

double SlaveInfo::get_linpack_time(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	auto &linpack_time = it->second.linpack_time;
	return double(linpack_time.count());
}

void SlaveInfo::sort_queue(deque<int> &slave_fd)
{
	CompareTimes cmp(this);
	sort(slave_fd.begin(), slave_fd.end(), cmp);
}


RunManagerYAMR::RunManagerYAMR(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, const string &_port, ofstream &_f_rmr)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
		_insfile_vec, _outfile_vec, stor_filename),
		 port(_port), f_rmr(_f_rmr)
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
	cout << "          starting YAMR (Yet Another run ManageR)..." << endl << endl;
	w_print_servinfo(servinfo, cout);
	cout << endl;
	//make socket, bind and listen
	addrinfo *connect_addr = w_bind_first_avl(servinfo, listener);
	if (connect_addr == nullptr)
	{
		stringstream err_str;
		err_str << "Error: port \"" << port << "\n is busy.  Can not bind port" << endl;
		throw(PestError(err_str.str()));
	}
	else {
		f_rmr << endl;
		cout << "YAMR Master listening on socket: " << w_get_addrinfo_string(connect_addr) << endl;
		f_rmr << "YAMR Master listening on socket:" << w_get_addrinfo_string(connect_addr) << endl;
	}
	w_listen(listener, BACKLOG);
	//free servinfo
	freeaddrinfo(servinfo);
	fdmax = listener;
	FD_ZERO(&master);
	FD_SET(listener, &master);
	cur_group_id = 0;
	return;
}

unordered_multimap<int, YamrModelRun>::iterator RunManagerYAMR::get_active_run_id(int socket)
{
	auto i = active_runs.begin();
	auto end = active_runs.end();

	for(; i!=end && i->second.get_socket() != socket; ++i)
	{}
	return i;
}

void RunManagerYAMR::initialize(const Parameters &model_pars, const Observations &obs)
{
	RunManagerAbstract::initialize(model_pars, obs);
	cur_group_id = NetPackage::get_new_group_id();
}

void RunManagerYAMR::reinitialize()
{
	RunManagerAbstract::reinitialize();
	cur_group_id = NetPackage::get_new_group_id();
}

void  RunManagerYAMR::free_memory()
{
	waiting_runs.clear();
	completed_runs.clear();
	zombie_runs.clear();
	failure_map.clear();
	file_stor.free_memory();
}

int RunManagerYAMR::add_run(const Parameters &model_pars)
{
	int run_id = file_stor.add_run(model_pars);
	YamrModelRun new_run(run_id);
	waiting_runs.push_back(new_run);
	return run_id;
}

int RunManagerYAMR::add_run(const std::vector<double> &model_pars)
{
	int run_id = file_stor.add_run(model_pars);
	YamrModelRun new_run(run_id);
	waiting_runs.push_back(new_run);
	return run_id;
}


void RunManagerYAMR::run()
{
	stringstream message;
	NetPackage net_pack;

	cout << "    running model " << waiting_runs.size() << " times" << endl;
	f_rmr << "running model " << waiting_runs.size() << " times" << endl;
	if(slave_info.size() == 0) // first entry is the listener, slave apper after this
	{
		cout << endl << "      waiting for slaves to appear..." << endl << endl;
		f_rmr << endl << "    waiting for slaves to appear..." << endl << endl;
	}
	while (!active_runs.empty() || !waiting_runs.empty())
	{
		init_slaves();
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


void RunManagerYAMR::listen()
{
	struct sockaddr_storage remote_addr;
	fd_set read_fds; // temp file descriptor list for select()
	socklen_t addr_len;
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
					slave_info.add(newfd);
				}
			}
			else  // handle data from a client
			{
				process_message(i);				
			} // END handle data from client
		} // END got new incoming connection
	} // END looping through file descriptors
}


void RunManagerYAMR::schedule_runs()
{
	NetPackage net_pack;
	slave_info.sort_queue(slave_fd);
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

void RunManagerYAMR::process_message(int i_sock)
{
	NetPackage net_pack;
	int err;
	vector<string> sock_name = w_getnameinfo_vec(i_sock);
	if(( err=net_pack.recv(i_sock)) <=0) // error or lost connection
	{
		if (err < 0) {
			f_rmr << "receive failed from slave: " << sock_name[0] <<":" <<sock_name[1] << " - terminating slave" << endl;
			cerr << "receive failed from slave: " << sock_name[0] <<":" <<sock_name[1] << " - terminating slave" << endl;
		}
		else {
			f_rmr << "lost connection to slave: " << sock_name[0] <<":" <<sock_name[1] << endl;
			cerr << "lost connection to slave: " << sock_name[0] <<":" <<sock_name[1] << endl;
		}
		w_close(i_sock); // bye!
		FD_CLR(i_sock, &master); // remove from master set
		slave_info.erase(i_sock); // remove information on this slave
		// remove run from active queue and return it to the waiting queue
		auto it_active = get_active_run_id(i_sock);
		if (it_active != active_runs.end())
		{
			YamrModelRun &cur_run = it_active->second;
			if (completed_runs.find(cur_run.get_id()) == completed_runs.end()) //check if run has already finish on another node
			{
				waiting_runs.push_front(it_active->second);
			}
			active_runs.erase(it_active);
		}
	}
	else if (net_pack.get_type() == NetPackage::PackType::RUNDIR)
	{
		string work_dir(net_pack.get_data().data(), net_pack.get_data().size());
		f_rmr << "New Slave connection from: " << sock_name[0] <<":" <<sock_name[1] << "   working dir: " << work_dir << endl;
		cout << endl << "New Slave connection from: " << sock_name[0] <<":" <<sock_name[1] << "   working dir: " << work_dir << endl << endl;
		slave_info.set_work_dir(i_sock, work_dir);
		slave_info.set_state(i_sock, SlaveInfo::State::CWD_RCV);
	}
	else if (net_pack.get_type() == NetPackage::PackType::LINPACK)
	{
		slave_info.end_linpack(i_sock);
		slave_info.set_state(i_sock, SlaveInfo::State::LINPACK_RCV);
	}
	else if (net_pack.get_type() == NetPackage::PackType::READY)
	{
		// ready message received from slave add slave to slave_fd
		slave_fd.push_back(i_sock);
	}

	else if (net_pack.get_type() == NetPackage::PackType::RUN_FINISH && net_pack.get_groud_id() != cur_group_id)
	{
		// this is an old run that did not finish on time
		// just ignore it
	}
	else if (net_pack.get_type() == NetPackage::PackType::RUN_FINISH)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_groud_id();
		// keep track of model run time
		slave_info.end_run(i_sock);
		f_rmr << "Run received from: " << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		cout << "Run received from: " << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		process_model_run(i_sock, net_pack);

	}
	else if (net_pack.get_type() == NetPackage::PackType::RUN_FAILED)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_groud_id();
		f_rmr << "Run failed on slave:" << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		cout << "Run failed on slave:" << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << group_id << ", run id = " << run_id << ")" << endl;
		file_stor.update_run_failed(run_id);
		failure_map.insert(make_pair(run_id, i_sock));
		auto it = get_active_run_id(i_sock);
		active_runs.erase(it);
		// TO DO add check for number of active nodes
		if (failure_map.count(run_id) < max_n_failure)
		{
			//put model run back into the waiting queue
			waiting_runs.push_front(YamrModelRun(run_id, i_sock));
		}
		else
		{
			//model run fail too many time.  Mark it as a bad run
			failed_runs.insert(run_id);
		}
	}
	else
	{
		f_rmr << "Received unsupported message from slave:" << endl;
		net_pack.print_header(f_rmr);
		//save results from model run
	}
}

 bool RunManagerYAMR::process_model_run(int sock_id, NetPackage &net_pack)
 {
	bool use_run = false;
	int run_id =net_pack.get_run_id();
	YamrModelRun model_run(run_id,  sock_id);
	//check if another instance of this model run has already completed 
	if (completed_runs.find(run_id) == completed_runs.end())
	{
		completed_runs.insert(pair<int, YamrModelRun>(run_id,  model_run));
		Parameters pars;
		Observations obs;
		Serialization::unserialize(net_pack.get_data(), pars, get_par_name_vec(), obs, get_obs_name_vec());
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

 bool RunManagerYAMR::schedule_run(int run_id)
 {
	bool scheduled = false;
	auto it_sock = slave_fd.end(); // iterator to current socket

	if (completed_runs.count(run_id) > 0)
	{
		// run already completed on different node.  Do nothing
	}
	else if (failure_map.count(run_id) == 0)
	{
		 // schedule a run on a slave
		NetPackage net_pack(NetPackage::PackType::START_RUN, cur_group_id, run_id, "");
		it_sock = slave_fd.begin();
	}
	else if (failure_map.count(run_id) > 0)
	{
		for(it_sock=slave_fd.begin(); it_sock!=slave_fd.end(); ++it_sock)
		{
			auto fail_iter_pair = failure_map.equal_range(run_id);

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
		YamrModelRun tmp_run(run_id, *it_sock);
		vector<char> data = file_stor.get_serial_pars(run_id);
		vector<string> sock_name = w_getnameinfo_vec(*it_sock);
		NetPackage net_pack(NetPackage::PackType::START_RUN, cur_group_id, run_id, "");
		int err = net_pack.send(*it_sock, &data[0], data.size());
		if (err != -1)
		{
			//start run timer
			slave_info.start_timer(*it_sock);
			f_rmr << "Sending run to: " << sock_name[0] <<":" <<sock_name[1] << "  (group id = " << cur_group_id << ", run id = " << run_id << ")" << endl;
			active_runs.insert(pair<int, YamrModelRun>(run_id, tmp_run));
			slave_fd.erase(it_sock);
			scheduled = true;
		}
	}
	return scheduled;
 }


 void RunManagerYAMR::init_slaves()
 {
	for (auto it_slv=slave_info.begin(); it_slv!=slave_info.end();++it_slv)
	{
		int i_sock = it_slv->first;
		SlaveInfo::State cur_state = slave_info.get_state(i_sock);
		if (cur_state == SlaveInfo::State::NEW)
		{
			NetPackage net_pack(NetPackage::PackType::REQ_RUNDIR, 0, 0, "");
			char data = '\0';
			int err = net_pack.send(i_sock, &data, sizeof(data));
			if (err != -1)
			{
				slave_info.set_state(i_sock, SlaveInfo::State::CWD_REQ);
			}
		}
		else if(cur_state == SlaveInfo::State::CWD_RCV)
		{
			// send Command line, tpl and ins information
			NetPackage net_pack(NetPackage::PackType::CMD, 0, 0, "");
			vector<char> data;
			vector<vector<string> const*> tmp_vec;
			tmp_vec.push_back(&comline_vec);
			tmp_vec.push_back(&tplfile_vec);
			tmp_vec.push_back(&inpfile_vec);
			tmp_vec.push_back(&insfile_vec);
			tmp_vec.push_back(&outfile_vec);
			tmp_vec.push_back(&file_stor.get_par_name_vec());
			tmp_vec.push_back(&file_stor.get_obs_name_vec());

			data = Serialization::serialize(tmp_vec);
			int err = net_pack.send(i_sock, &data[0], data.size());
			if (err != -1)
			{
				slave_info.set_state(i_sock, SlaveInfo::State::CMD_SENT);
			}
		}
		else if(cur_state == SlaveInfo::State::CMD_SENT)
		{
			NetPackage net_pack(NetPackage::PackType::REQ_LINPACK, 0, 0, "");
			char data = '\0';
			int err = net_pack.send(i_sock, &data, sizeof(data));
			if (err != -1)
			{
				slave_info.set_state(i_sock, SlaveInfo::State::LINPACK_REQ);
				slave_info.start_timer(i_sock);
			}
		}
		else if(cur_state == SlaveInfo::State::LINPACK_RCV)
		{
			slave_info.set_state(i_sock, SlaveInfo::State::ACTIVE);
			slave_fd.push_back(i_sock);
		}
	}
 }

RunManagerYAMR::~RunManagerYAMR(void)
{
	//close sockets and cleanup
	int err;
	err = w_close(listener);
	FD_CLR(listener, &master);
	// this is needed to ensure that the first slave closes properly
	w_sleep(2000);
	for(int i = 0; i <= fdmax; i++) {
		if (FD_ISSET(i, &master)) 
		{
			NetPackage netpack(NetPackage::PackType::TERMINATE, 0, 0,"");
			char data;
			netpack.send(i, &data, 0);
			err = w_close(i);
			FD_CLR(i, &master);
		}
	}
	w_cleanup();
}
