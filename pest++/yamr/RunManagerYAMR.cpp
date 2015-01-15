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
#include <chrono>
#include <ctime>
#include <iomanip>
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
	last_ping_time = std::chrono::system_clock::now();
	ping = false;
	failed_pings = 0;
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
		ret = ! (itb->second.run_time > ita->second.run_time);
	}
	else if (itb->second.linpack_time > std::chrono::milliseconds(0) && ita->second.linpack_time > std::chrono::milliseconds(0))
	{
		ret = ! (itb->second.linpack_time > ita->second.linpack_time);
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

double SlaveInfo::get_duration_sec(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());	
	chrono::system_clock::duration dt = chrono::system_clock::now() - it->second.start_time;
	return (double)std::chrono::duration_cast<std::chrono::milliseconds>(dt).count() / 1000.0;
}

double SlaveInfo::get_duration_minute(int sock_id)
{
	return this->get_duration_sec(sock_id) / 60.0;
}

double SlaveInfo::get_runtime_sec(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	auto &run_time = it->second.run_time;
	return(double)std::chrono::duration_cast<std::chrono::milliseconds>(run_time).count() / 1000.0;
}

double SlaveInfo::get_runtime_minute(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	auto &run_time = it->second.run_time;
	double run_minutes = std::chrono::duration_cast<std::chrono::milliseconds>(run_time).count() / 60000.0;
	return run_minutes;
}

double SlaveInfo::get_global_runtime_minute()
{
	double global_runtime = 0;
	double temp = 0;
	int count = 0;
	for (auto &si : slave_info_map)
	{
		temp = std::chrono::duration_cast<std::chrono::milliseconds>(si.second.run_time).count() / 60000.0;
		if (temp > 0)
		{
			count++;
			global_runtime += temp;
		}
	}
	return global_runtime / (double)count;
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

void SlaveInfo::reset_failed_pings(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	it->second.failed_pings = 0;
}

int SlaveInfo::add_failed_ping(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	it->second.failed_pings++;
	return it->second.failed_pings;
}

void SlaveInfo::set_ping(int sock_id, bool val)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	it->second.ping = val;
	//a success response
	if (!val) reset_failed_pings(sock_id);
	//sending a request
	else reset_last_ping_time(sock_id);
}

bool SlaveInfo::get_ping(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	return it->second.ping;
}

void SlaveInfo::reset_last_ping_time(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	it->second.last_ping_time = chrono::system_clock::now();
}

int SlaveInfo::seconds_since_last_ping_time(int sock_id)
{
	auto it = slave_info_map.find(sock_id);
	assert(it != slave_info_map.end());
	return chrono::duration_cast<std::chrono::seconds>
		(chrono::system_clock::now() - it->second.last_ping_time).count();
}


RunManagerYAMR::RunManagerYAMR(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, const string &_port, ofstream &_f_rmr, int _max_n_failure)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
	_insfile_vec, _outfile_vec, stor_filename, _max_n_failure),
	port(_port), f_rmr(_f_rmr)
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

unordered_multimap<int, YamrModelRun>::iterator RunManagerYAMR::get_zombie_run_id(int socket)
{
	auto i = zombie_runs.begin();
	auto end = zombie_runs.end();

	for (; i != end && i->second.get_socket() != socket; ++i)
	{
	}
	return i;
}


void RunManagerYAMR::initialize(const Parameters &model_pars, const Observations &obs, const string &_filename)
{
	RunManagerAbstract::initialize(model_pars, obs, _filename);
	cur_group_id = NetPackage::get_new_group_id();
}

void RunManagerYAMR::initialize_restart(const std::string &_filename)
{
	file_stor.init_restart(_filename);
	waiting_runs.clear();
	completed_runs.clear();
	zombie_runs.clear();
	failure_map.clear();
	vector<int> waiting_run_id_vec = get_outstanding_run_ids();
	for (int &id : waiting_run_id_vec)
	{
		YamrModelRun new_run(id);
		waiting_runs.push_back(new_run);
	}
}

void RunManagerYAMR::reinitialize(const std::string &_filename)
{
	waiting_runs.clear();
	completed_runs.clear();
	zombie_runs.clear();
	failure_map.clear();
	concurrent_map.clear();
	RunManagerAbstract::reinitialize(_filename);
	cur_group_id = NetPackage::get_new_group_id();
}

void  RunManagerYAMR::free_memory()
{
	waiting_runs.clear();
	completed_runs.clear();
	zombie_runs.clear();
	failure_map.clear();
}

int RunManagerYAMR::add_run(const Parameters &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	YamrModelRun new_run(run_id);
	waiting_runs.push_back(new_run);
	return run_id;
}

int RunManagerYAMR::add_run(const std::vector<double> &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	YamrModelRun new_run(run_id);
	waiting_runs.push_back(new_run);
	return run_id;
}

int RunManagerYAMR::add_run(const Eigen::VectorXd &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	YamrModelRun new_run(run_id);
	waiting_runs.push_back(new_run);
	return run_id;
}

void RunManagerYAMR::run()
{
	stringstream message;
	NetPackage net_pack;
	model_runs_done = 0;
	model_runs_failed = 0;
	cout << "    running model " << waiting_runs.size() << " times" << endl;
	f_rmr << "running model " << waiting_runs.size() << " times" << endl;
	if(slave_info.size() == 0) // first entry is the listener, slave apper after this
	{
		cout << endl << "      waiting for slaves to appear..." << endl << endl;
		f_rmr << endl << "    waiting for slaves to appear..." << endl << endl;
	}
	cout << endl;
	f_rmr << endl;
	while (!active_runs.empty() || !waiting_runs.empty() || !zombie_runs.empty())
	{
		init_slaves();
		//schedule runs on available nodes
		schedule_runs();
		// get and process incomming messages
		listen();		
	}
	total_runs += model_runs_done;
	echo();
	//kill any remaining active runs
	message.str("");
	message  << "    " << completed_runs.size() << " runs complete";
	cout << endl << "---------------------" << endl << message.str() << endl << endl;
	f_rmr << endl << "---------------------" << endl << message.str() << endl << endl;
	concurrent_map.clear();
	//if (init_run_obs.size() == 0)
	//	int status = file_stor.get_observations(0, init_run_obs);
	if (init_sim.size() == 0)
	{
		vector<double> pars;
		int status = file_stor.get_run(0, pars, init_sim);
	}
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
				//set the ping flag since the slave sent something back
				slave_info.set_ping(i, false);
				process_message(i);				
			} // END handle data from client
		} // END got new incoming connection
		else
		{						
			for (auto it_slv = slave_info.begin(); it_slv != slave_info.end(); ++it_slv)
			{
				if ((i == it_slv->first) && (slave_info.get_state(i) == SlaveInfo::State::ACTIVE))
				{					
					ping(i);
					break;
				}
			}			
		}
	} // END looping through file descriptors
}

void RunManagerYAMR::ping(int i_sock)
{				
	vector<string> sock_name = w_getnameinfo_vec(i_sock);
	fd_set read_fds = master;
	//if the slave hasn't communicated since the last ping request
	if ((!FD_ISSET(i_sock,&read_fds)) && (slave_info.get_ping(i_sock)))
	{
		int fails = slave_info.add_failed_ping(i_sock);		
		report("failed to receive ping response from slave: "+sock_name[0]+"$"+slave_info.get_work_dir(i_sock),false);		
		if (fails >= MAX_FAILED_PINGS)
		{
			report("max failed ping communications since last successful run for slave:"+sock_name[0]+"$"+slave_info.get_work_dir(i_sock)+"  -> terminating",false);
			close_slave(i_sock);
			return;
		}		
	}
	//check if it is time to ping again...
	double duration = (double)slave_info.seconds_since_last_ping_time(i_sock);
	double ping_time = max(double(PING_INTERVAL_SECS), slave_info.get_runtime_sec(i_sock));
	if (duration >= ping_time)
	{		
		const char* data = "\0";
		NetPackage net_pack(NetPackage::PackType::PING, 0, 0, "");
		int err = net_pack.send(i_sock, data, 0);
		if (err <= 0)
		{
			int fails = slave_info.add_failed_ping(i_sock);			
			report("failed to send ping request to slave:"+sock_name[0]+"$"+slave_info.get_work_dir(i_sock),false);			
			if (fails >= MAX_FAILED_PINGS)
			{
				report("max failed ping communications since last successful run for slave:" + sock_name[0] + "$" + slave_info.get_work_dir(i_sock) + "  -> terminating", true);
				close_slave(i_sock);
				return;
			}
		}
		else slave_info.set_ping(i_sock, true);
#ifdef _DEBUG
		//report("ping sent to slave:" + sock_name[0] + "$" + slave_info.get_work_dir(i_sock), false);
#endif
	}
}

void RunManagerYAMR::close_slave(int i_sock)
{	
	vector<string> sock_name = w_getnameinfo_vec(i_sock);
	w_close(i_sock); // bye!
	FD_CLR(i_sock, &master); // remove from master set
	slave_info.erase(i_sock); // remove information on this slave
	//remove slave from slave_fd 
	auto it_sfd = find(slave_fd.begin(), slave_fd.end(), i_sock);
	if (it_sfd != slave_fd.end())
	{
		slave_fd.erase(it_sfd);
	}
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
	stringstream ss;
	
	ss << "closed connection to slave: " << sock_name[0] << ":" << sock_name[1] << "; number of slaves: " << slave_info.size();
	report(ss.str(), false);
	//cout << "---closed slave connection to slave " << sock_name[0] << "$" << sock_name[1] << "---" << endl;
}

bool RunManagerYAMR::schedule_run(int run_id)
{
	bool scheduled = false;
	auto it_sock = slave_fd.end(); // iterator to current socket

	if (completed_runs.count(run_id) > 0)
	{
		// run already completed on different node.  Do nothing
	}
	else if (failure_map.count(run_id) == 0 || failure_map.count(run_id) >= slave_fd.size())
	{
		// schedule a run on a slave
		NetPackage net_pack(NetPackage::PackType::START_RUN, cur_group_id, run_id, "");
		it_sock = slave_fd.begin();
	}
	else if (failure_map.count(run_id) > 0)
	{
		for (it_sock = slave_fd.begin(); it_sock != slave_fd.end(); ++it_sock)
		{
			auto fail_iter_pair = failure_map.equal_range(run_id);

			auto i = fail_iter_pair.first;
			for (i = fail_iter_pair.first;
				i != fail_iter_pair.second && i->second != *it_sock;
				++i) {
			}
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
			int concur;
			auto it_concur = concurrent_map.find(run_id);
			if (it_concur != concurrent_map.end())
			{
				concurrent_map[it_concur->first]++;
				concur = it_concur->second;
			}
			else
			{
				concur = 1;
				concurrent_map.insert(pair<int, int>(run_id, concur));				
			}

			stringstream ss;
			ss << "Sending run " << run_id << " to: " << sock_name[0] << "$" << slave_info.get_work_dir(*it_sock) << "  (group id = " << cur_group_id << ", run id = " << run_id << ", concurrent runs = " << concur << ")";
			report(ss.str(), false);
			active_runs.insert(pair<int, YamrModelRun>(run_id, tmp_run));
			//reset the last ping time so we don't ping immediately after run is started
			slave_info.reset_last_ping_time(*it_sock);
			slave_fd.erase(it_sock);
			scheduled = true;			
		}
	}
	return scheduled;
}

void RunManagerYAMR::schedule_runs()
{
	NetPackage net_pack;
	slave_info.sort_queue(slave_fd);
	
	//first try to schedule waiting runs
	for (auto it_run = waiting_runs.begin(); !slave_fd.empty() && it_run != waiting_runs.end();)
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

	//check for overdue runs
	double duration, avg_runtime;
	double global_avg_runtime = slave_info.get_global_runtime_minute();
	bool should_schedule = false;
	for (auto it_active = active_runs.begin(); !slave_fd.empty() && it_active != active_runs.end();++it_active)
	{		
		should_schedule = false;
		int act_sock_id = it_active->second.get_socket();
		int run_id = it_active->second.get_id();
		unordered_map<int, int>::iterator it_concur;
		duration = slave_info.get_duration_minute(act_sock_id);
		avg_runtime = slave_info.get_runtime_minute(act_sock_id);
		if (avg_runtime <= 0) avg_runtime = global_avg_runtime;
		if (avg_runtime <= 0) avg_runtime = 1.0E+10;
		if (duration > avg_runtime*PERCENT_OVERDUE_GIVEUP)
		{
			vector<string> sock_name = w_getnameinfo_vec(act_sock_id);
			stringstream ss;
			ss << "killing overdue run " << run_id << " (" << duration << "|" << avg_runtime <<
				" minutes) on: " << sock_name[0] << "$" << slave_info.get_work_dir(act_sock_id);
 			report(ss.str(), false);
			NetPackage net_pack(NetPackage::PackType::REQ_KILL, 0, 0, "");
			char data = '\0';
			int err = net_pack.send(act_sock_id, &data, sizeof(data));
			if (err <= 0)
			{
				report("error sending kill request to slave:" + sock_name[0] + "$" +
					slave_info.get_work_dir(act_sock_id), true);
				close_slave(act_sock_id);
			}			
		}
		if (duration > avg_runtime*PERCENT_OVERDUE_RESCHED)
		{
			//check how many concurrent runs are going			
			it_concur = concurrent_map.find(it_active->first);
			if (it_concur == concurrent_map.end()) throw PestError("active run id not found in concurrent map");
			if (it_concur->second < MAX_CONCURRENT_RUNS) should_schedule = true;			
		}

		if (should_schedule)
		{
			vector<string> sock_name = w_getnameinfo_vec(act_sock_id);
			stringstream ss;
			int run_id = it_active->second.get_id();
			ss << "rescheduling overdue run " << run_id << " (" << duration << "|" <<
				avg_runtime << " minutes) on: " << sock_name[0] << "$" << 
				slave_info.get_work_dir(act_sock_id);
			report(ss.str(), false);
			bool success = schedule_run(run_id);
			if (success)
			{
				stringstream ss;
				ss << concurrent_map[it_concur->first] << " concurrent runs for run id = " << run_id;
				report(ss.str(), false);
			}
			else
			{
				stringstream ss;
				ss << "failed to schedule concurrent run for run id = " << run_id;
				report(ss.str(), false);
			}
		}
	}
}

void RunManagerYAMR::echo()
{
	cout << get_time_string() << "->" << setw(6) << model_runs_done << " runs complete " <<
		setw(6) << model_runs_failed << " runs failed " <<
		setw(4) << slave_info.size() << " slaves\r" << flush;

}

string RunManagerYAMR::get_time_string()
{
	std::time_t tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	string t_str = ctime(&tt);
	return t_str.substr(0, t_str.length() - 1);

}

void RunManagerYAMR::report(std::string message,bool to_cout)
{
	string t_str = get_time_string();
	f_rmr << t_str << "->" << message << endl;
	if (to_cout) cout << endl << t_str << "->" << message << endl;
}

void RunManagerYAMR::process_message(int i_sock)
{	
	echo();
	NetPackage net_pack;
	int err;
	vector<string> sock_name = w_getnameinfo_vec(i_sock);
	//std::time_t tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	if(( err=net_pack.recv(i_sock)) <=0) // error or lost connection
	{
		if (err < 0) {
			report("receive failed from slave: "+sock_name[0]+"$"+slave_info.get_work_dir(i_sock)+" - terminating slave",false);			
		}
		else {
			report("lost connection to slave: " + sock_name[0] + "$" + slave_info.get_work_dir(i_sock),false);
		}
		close_slave(i_sock);
	}
	else if (net_pack.get_type() == NetPackage::PackType::RUNDIR)
	{
		string work_dir(net_pack.get_data().data(), net_pack.get_data().size());
		stringstream ss;
		ss << "initializing new slave connection from: " << sock_name[0] << ":" << sock_name[1] << "; number of slaves: " << slave_info.size() << "; working dir: " << work_dir;
		report(ss.str(),false);		
		slave_info.set_work_dir(i_sock, work_dir);
		slave_info.set_state(i_sock, SlaveInfo::State::CWD_RCV);
	}
	else if (net_pack.get_type() == NetPackage::PackType::LINPACK)
	{
		slave_info.end_linpack(i_sock);
		slave_info.set_state(i_sock, SlaveInfo::State::LINPACK_RCV);
		stringstream ss;
		ss << "new slave ready: " << sock_name[0] << ":" << sock_name[1];
		report(ss.str(), false);
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
		//should never have this situation with the threaded slave
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_groud_id();
		stringstream ss;
		ss << "run " << run_id << " received from unexpected group id: " << group_id << ", should be group: " << cur_group_id;
		throw PestError(ss.str());
	}
	else if (net_pack.get_type() == NetPackage::PackType::RUN_FINISH)
	{		
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_groud_id();
		//decrement concurrent run counter
		auto it_concur = concurrent_map.find(run_id);
		
		if (it_concur == concurrent_map.end()) throw PestError("active run id not in concurrent map");		
		concurrent_map[it_concur->first]--;
		int concur = concurrent_map[it_concur->first];
		

		//check if this was a zombie run
		auto it_zombie = get_zombie_run_id(i_sock);
		if (it_zombie != zombie_runs.end())
		{
			stringstream ss;
			ss << "zombie run " << run_id << " finished on: " << sock_name[0] << "$" << slave_info.get_work_dir(i_sock) << 
				"  (run time = " << slave_info.get_runtime_minute(i_sock) << " min, group id = " << group_id <<
				", run id = " << run_id << " concurrent = " << concur << ")";
			report(ss.str(), false);
			zombie_runs.erase(it_zombie);
		}
		else
		{
			// keep track of model run time
			slave_info.end_run(i_sock);
			stringstream ss;
			ss << "run " << run_id << " received from: " << sock_name[0] << "$" << slave_info.get_work_dir(i_sock) << 
				"  (run time = " << slave_info.get_runtime_minute(i_sock) << " min, group id = " << group_id <<
				", run id = " << run_id << " concurrent = " << concur << ")";
			report(ss.str(), false);
			process_model_run(i_sock, net_pack);
		}

	}
	else if (net_pack.get_type() == NetPackage::PackType::RUN_FAILED)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_groud_id();
		//check if this was a concurrent run - decrement
		auto it_concur = concurrent_map.find(run_id);		
		if (it_concur == concurrent_map.end()) throw PestError("active run id not in concurrent map");
		concurrent_map[it_concur->first]--;
		int concur = concurrent_map[it_concur->first];
		// remove run from active queue and return it to the waiting queue
		auto it_zombie = get_zombie_run_id(i_sock);
		if (it_zombie != zombie_runs.end())
		{
			stringstream ss;
			ss << "zombie run " << run_id << " killed on slave: " << sock_name[0] << "$" << slave_info.get_work_dir(i_sock) << ", run id = " << run_id << " concurrent = " << concur;
			report(ss.str(),false);
			zombie_runs.erase(it_zombie);
		}
		else
		{
			stringstream ss;
			ss << "Run " << run_id << " failed on slave:" << sock_name[0] << "$" << slave_info.get_work_dir(i_sock) << "  (group id = " << group_id << ", run id = " << run_id << ", concurrent = " << concur << ") ";
			report(ss.str(), false);
			model_runs_failed++;
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
		}
	}
	else if (net_pack.get_type() == NetPackage::PackType::PING)
	{
#ifdef _DEBUG
		//report("ping received from slave" + sock_name[0] + "$" + slave_info.get_work_dir(i_sock),false);
#endif
	}
	else if (net_pack.get_type() == NetPackage::PackType::IO_ERROR)
	{
		//string err(net_pack.get_data().begin(),net_pack.get_data().end());		
		report("error in model IO files on slave: " + sock_name[0] + "$" + slave_info.get_work_dir(i_sock) + "-terminating slave. ",true);
		close_slave(i_sock);
	}
	else
	{
		report("received unsupported message from slave: ", false);
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
		model_runs_done++;
		//beopest-style screen output for run counting		
		//cout << setw(7) << model_runs_done;
		//if (model_runs_done % 9 == 0) cout << endl;				
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
	//kill all zombies
	for (auto it_zombie = zombie_runs.begin(); it_zombie != zombie_runs.end(); ++it_zombie)
	{
		if (it_zombie->second.get_id() == run_id)
		{

			int zombie_id = it_zombie->second.get_socket();
			vector<string> sock_name = w_getnameinfo_vec(zombie_id);
			stringstream ss;
			ss << "killing zombie run " << run_id << " on slave : " << sock_name[0] << "$" << slave_info.get_work_dir(zombie_id);
			report(ss.str(), false);
			NetPackage net_pack(NetPackage::PackType::REQ_KILL, 0, 0, "");
			char data = '\0';
			int err = net_pack.send(zombie_id, &data, sizeof(data));
			if (err <= 0)
			{
				report("error sending kill request to slave:" + sock_name[0] + "$" + slave_info.get_work_dir(zombie_id), true);
			}
		}
	}
	return use_run;
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
