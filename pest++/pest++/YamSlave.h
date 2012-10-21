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

#ifndef YAMSLAVE_H_
#define YAMSLAVE_H_

#include <winsock2.h>
#include <ws2tcpip.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <fstream>
#include <memory>
#include "Pest.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include "network_wrapper.h"
#include "network_package.h"
#include "Transformable.h"
#include "FileManager.h"

class YAMSlave{
public:
	void init_network(const std::string &host, const std::string &port);
	void start(const std::string &host, const std::string &port);
	~YAMSlave();
	void run();
	void recv_message(NetPackage &net_pack);
	int send_message(NetPackage &net_pack, const void *data=NULL, unsigned long data_len=0);
	int run_model(Parameters &pars, Observations &obs);
	std::string tpl_err_msg(int i);
	std::string YAMSlave::ins_err_msg(int i);
private:
	int sockfd;
	int fdmax;
	fd_set master;
	std::vector<std::string> comline_vec;
	std::vector<std::string> tplfile_vec;
	std::vector<std::string> inpfile_vec;
	std::vector<std::string> insfile_vec;
	std::vector<std::string> outfile_vec;
	std::vector<std::string> obs_name_vec;
};

#endif /* YAMSLAVE_H_ */
