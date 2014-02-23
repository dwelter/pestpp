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
#ifndef RUNMANAGERYAMR_H
#define RUNMANAGERYAMR_H
#include "network_wrapper.h"
#include <string>
#include <set>
#include <deque>
#include <unordered_map>
#include <chrono>
#include "network_wrapper.h"
#include "network_package.h"
#include "RunManagerAbstract.h"
#include "RunStorage.h"

class YamrModelRun
{
public:
	YamrModelRun(int _run_id, int _sockfd=0);
	int get_id(){return run_id;}
	void set_socket(int _sockfd) {sockfd = _sockfd;}
	int get_socket() const  {return sockfd;}
	~YamrModelRun() {}
private:
	int sockfd;
	int run_id;
};


class SlaveInfo
{
public:
	enum class State {NEW, CWD_REQ, CWD_RCV, CMD_SENT, LINPACK_REQ, LINPACK_RCV, ACTIVE};
	class SlaveRec {
		public:
			friend SlaveInfo;
			SlaveRec();
			~SlaveRec(){}
		private:
			State state;
			std::chrono::system_clock::duration linpack_time;
			std::chrono::system_clock::duration run_time;
			std::chrono::system_clock::time_point start_time;
			std::string work_dir;
		};
	typedef std::unordered_map<int, SlaveRec>::iterator iterator;
	typedef std::unordered_map<int, SlaveRec>::const_iterator const_iterator;
	SlaveInfo::iterator begin(){return slave_info_map.begin();}
	SlaveInfo::const_iterator begin() const {return slave_info_map.begin();}
	SlaveInfo::iterator end() {return slave_info_map.end();}
	SlaveInfo::const_iterator end() const {return slave_info_map.end();}
	SlaveInfo();
	size_t size() const;
	void add(int sock_id);
	void erase(int sock_id);
	State get_state(int sock_id);
	void set_state(int sock_id, const State);
	void set_work_dir(int sock_id, const std::string & wkd);
	std::string get_work_dir(int sock_id) const;
	void start_timer(int sock_id);
	void end_run(int sock_id);
	void end_linpack(int sock_id);
	double get_runtime(int sock_id);
	double get_linpack_time(int sock_id);
	void sort_queue(std::deque<int> &slave_fd);
	~SlaveInfo();
private:
	class CompareTimes
	{
	public:
		CompareTimes(SlaveInfo *_my_class_ptr) : my_class_ptr(_my_class_ptr) {}
		bool operator() (int a, int b);
		SlaveInfo *my_class_ptr;
	};

	std::unordered_map<int, SlaveRec> slave_info_map;
};

class RunManagerYAMR : public RunManagerAbstract
{
public:
	RunManagerYAMR(const std::vector<std::string> _comline_vec,
		const std::vector<std::string> _tplfile_vec, const std::vector<std::string> _inpfile_vec,
		const std::vector<std::string> _insfile_vec, const std::vector<std::string> _outfile_vec,
		const std::string &stor_filename, const std::string &port, std::ofstream &_f_rmr);
	virtual void initialize(const Parameters &model_pars, const Observations &obs, const std::string &_filename = std::string(""));
	virtual void initialize_restart(const std::string &_filename);
	virtual void reinitialize(const std::string &_filename = std::string(""));
	virtual void free_memory();
	virtual int add_run(const Parameters &model_pars, const std::string &info_txt="", double info_value=RunStorage::no_data);
	virtual int add_run(const std::vector<double> &model_pars, const std::string &info_txt="", double info_valuee=RunStorage::no_data);
	virtual int add_run(const Eigen::VectorXd &model_pars, const std::string &info_txt="", double info_valuee=RunStorage::no_data);
	virtual void run();
	~RunManagerYAMR(void);
private:
	std::string port;
	static const int BACKLOG = 10;
	int listener;
	int fdmax;
	static const int yamr_max_n_failure = 3; // maximium number of times to retry a failed model run
	std::deque<int> slave_fd; // list of slaves ready to accept a model run
	fd_set master; // master file descriptor list
	std::deque<YamrModelRun> waiting_runs;
	std::ofstream &f_rmr;
	std::unordered_multimap<int, YamrModelRun> active_runs;
	std::unordered_multimap<int, YamrModelRun> zombie_runs;
	std::unordered_map<int, YamrModelRun> completed_runs;
	SlaveInfo slave_info;
	std::unordered_multimap<int, int> failure_map;
	void listen();
	bool process_model_run(int sock_id, NetPackage &net_pack);
	void process_message(int i);
	bool schedule_run(int run_id);
	void schedule_runs();
	void init_slaves();
	std::unordered_multimap<int, YamrModelRun>::iterator get_active_run_id(int socket);
};

#endif /* RUNMANAGERYAMR_H */