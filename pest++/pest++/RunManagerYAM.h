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
#ifndef RUNMANAGERYAM_H
#define RUNMANAGERYAM_H
#include <winsock2.h>
#include <string>
#include <set>
#include <deque>
#include <unordered_map>
#include "network_wrapper.h"
#include "network_package.h"
#include "RunManagerAbstract.h"
#include "Transformable.h"
#include "RunStorage.h"

class YamModelRun
{
public:
	YamModelRun(int _run_id, int _sockfd=0);
	int get_id(){return run_id;}
	void set_socket(int _sockfd) {sockfd = _sockfd;}
	int get_socket() const  {return sockfd;}
	~YamModelRun() {}
private:
	int sockfd;
	int run_id;
};


class RunManagerYAM : public RunManagerAbstract
{
public:
	RunManagerYAM(const ModelExecInfo &_mode_exec_info, const std::string &port, const std::string &stor_filename, ofstream &_f_rmr);
	virtual void allocate_memory(const Parameters &pars, const Observations &obs, int _nruns);
	virtual void free_memory();
	virtual int add_run(const Parameters &model_pars);
	virtual void run();
	virtual void get_run(ModelRun &model_run, int run_id, PAR_UPDATE update_type=DEFAULT_PAR_UPDATE);
	virtual Parameters get_model_parameters(int run_id) const;
	virtual int get_total_runs(void) const {return total_runs;}
	virtual int get_nruns() {return nruns;}
	void RunManagerYAM::listen();
	~RunManagerYAM(void);
private:
	std::string port;
	static const int BACKLOG = 10;
	int listener;
	int fdmax;
	deque<int> slave_fd;
	fd_set master; // master file descriptor list
	std::deque<YamModelRun> waiting_runs;
	RunStorage file_stor;
	ofstream &f_rmr;
	unordered_multimap<int, YamModelRun> active_runs;
	unordered_multimap<int, YamModelRun> zombie_runs;
	unordered_map<int, YamModelRun> completed_runs;
	bool process_model_run(YamModelRun model_run);
	void process_message(int i);
	void schedule_runs(int group_id);
};

#endif /* RUNMANAGERYAM_H */