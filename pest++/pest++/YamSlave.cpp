#include "YamSlave.h"
#include "utilities.h"
#include "Serialization.h"
#include "system_variables.h"

using namespace pest_utils;

extern "C"
{

	void WRTTPL(int *,
		char *,
		char *,
		int *,
		char *,
		double *,
		int *);

	void READINS(int *,
		char *,
		char *,
		int *,
		char *,
		double *,
		int *);
}



void YAMSlave::init_network(const string &host, const string &port)
{
	w_init();


	int status;
	int err;
	struct addrinfo hints;
	struct addrinfo *servinfo;

	memset(&hints, 0, sizeof hints);
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;

	status = w_getaddrinfo(host.c_str(), port.c_str(), &hints, &servinfo);
	w_print_servinfo(servinfo, cout);
	// connect
	sockfd = w_socket(servinfo->ai_family, servinfo->ai_socktype, servinfo->ai_protocol);
	for (err = -1; err == -1;)
	{
		err = w_connect(sockfd, servinfo->ai_addr, servinfo->ai_addrlen);
	}
	freeaddrinfo(servinfo);

	fdmax = sockfd;
	FD_ZERO(&master);
	FD_SET(sockfd, &master);
	// send run directory to master
	cout << endl;
	cout << "connection to master succeeded..." << endl << endl;
	NetPackage net_pack(NetPackage::RUN_DIR, 0, 0,"");
	string cwd =  OperSys::getcwd();
	for (err = -1; err==-1;)
	{
		err = send_message(net_pack, cwd.c_str(), cwd.size());
	}
	 // Send READY Message to master
	//NetPackage net_pack(NetPackage::READY, 0, 0,"");
	//char data;
	//while (err == -1)
	//{
	//	err = send_message(net_pack, &data, 0);
	//}
}


YAMSlave::~YAMSlave()
{
	w_close(sockfd);
	w_cleanup();
}

void YAMSlave::recv_message(NetPackage &net_pack)
{
	fd_set read_fds;
	for(;;) {
	read_fds = master; // copy master
		if (w_select(fdmax+1, &read_fds, NULL, NULL, NULL) == -1) {
			exit(4);
		}
		for(int i = 0; i <= fdmax; i++) {
			if (FD_ISSET(i, &read_fds)) { // got message to read
				int err;
				if(( err=net_pack.recv(i)) <=0) // error or lost connection
				{
					vector<string> sock_name = w_getnameinfo_vec(i);
					if (err < 0) {
							cerr << "receive from master failed: " << sock_name[0] <<":" <<sock_name[1] << endl;
					}
					else {
						cerr << "lost connection to master: " << sock_name[0] <<":" <<sock_name[1] << endl;
						w_close(i); // bye!
						FD_CLR(i, &master); // remove from master set
					}
				}
				else  
				{
					// received data sored in net_pack return to calling routine to process it
				}
			return;
			}
		}
	}
}


int YAMSlave::send_message(NetPackage &net_pack, const void *data, unsigned long data_len)
{
	return net_pack.send(sockfd, data, data_len);
}

string YAMSlave::tpl_err_msg(int i)
{
	string err_msg;
	switch (i)
	{
	case 0:
		err_msg = "Routine completed successfully";
		break;
	case 1:
		err_msg = "TPL file does not exist";
		break;
	case 2:
		err_msg = "Not all parameters listed in TPL file";
		break;
	case 3:
		err_msg = "Illegal header specified in TPL file";
		break;
	case 4:
		err_msg = "Error getting parameter name from template";
		break;
	case 5:
		err_msg = "Error getting parameter name from template";
		break;
	case 10:
		err_msg = "Error writing to model input file";
	}
	return err_msg;
}


string YAMSlave::ins_err_msg(int i)
{
	string err_msg;
	switch (i)
	{
	case 0:
		err_msg = "Routine completed successfully";
		break;
	default:
		err_msg = "";
	}
	return err_msg;
}

int YAMSlave::run_model(Parameters &pars, Observations &obs)
{
	try 
	{
	//	message.str("");
		int ifail;
		int ntpl = tplfile_vec.size();
		int npar = pars.size();
		vector<string> par_name_vec;
		vector<double> par_values;
		for(auto &i : pars)
		{
			par_name_vec.push_back(i.first);
			par_values.push_back(i.second);
		}
		WRTTPL(&ntpl, StringvecFortranCharArray(tplfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
			StringvecFortranCharArray(inpfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
			&npar, StringvecFortranCharArray(par_name_vec, 50, pest_utils::TO_LOWER).get_prt(),
			par_values.data(), &ifail);
		if(ifail != 0)
		{
			throw PestError("Error processing template file");
		}
		// update parameter values
		pars.clear();
		for (int i=0; i<npar; ++i)
		{
			pars[par_name_vec[i]] = par_values[i];
		}

		// run model
		for (auto &i : comline_vec)
		{
			system(i.c_str());
		}
		// process instructio files
		int nins = insfile_vec.size();
		int nobs = obs_name_vec.size();
		std::vector<double> obs_vec;
		obs_vec.resize(nobs, -9999.00);
		READINS(&nins, StringvecFortranCharArray(insfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
			StringvecFortranCharArray(outfile_vec, 50, pest_utils::TO_LOWER).get_prt(),
			&nobs, StringvecFortranCharArray(obs_name_vec, 50, pest_utils::TO_LOWER).get_prt(),
			obs_vec.data(), &ifail);
		if(ifail != 0)
		{
			throw PestError("Error processing template file");
		}
		// update observation values
		obs.clear();
		for (int i=0; i<nobs; ++i)
		{
			obs[obs_name_vec[i]] = obs_vec[i];
		}

	}
	catch(...)
	{
		cerr << "   Aborting model run" << endl;
	}
	return 1;
}


void YAMSlave::init(const string &host, const string &port, const vector<string> _comline_vec,
		const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
		const vector<string> _insfile_vec, const vector<string> _outfile_vec,
		const vector<string> _obs_name_vec)
{
	comline_vec = _comline_vec;
	tplfile_vec = _tplfile_vec;
	inpfile_vec = _inpfile_vec;
	insfile_vec = _insfile_vec;
	outfile_vec = _outfile_vec;
	obs_name_vec = _obs_name_vec;

	bool terminate = false;
	NetPackage net_pack;
	Observations obs;
	Parameters pars;
	vector<char> serialized_data;

	init_network(host, port);

	while (!terminate)
	{
		//get message from master
		recv_message(net_pack);
		if (net_pack.get_type() == NetPackage::START_RUN)
		{
			Serialization::unserialize(net_pack.get_data(), pars);
			// run model
			int group_id = net_pack.get_groud_id();
			int run_id = net_pack.get_run_id();
			cout << "received parameters (group id = " << group_id << ", run id = " << run_id << ")" << endl;
			cout << "starting model run..." << endl;
			run_model(pars, obs);
			//send model results back
			cout << "run complete" << endl;
			cout << "sending results to master (group id = " << group_id << ", run id = " << run_id << ")" <<endl << endl;
			serialized_data = Serialization::serialize(pars, obs);
			net_pack.reset(NetPackage::RUN_FINISH, net_pack.get_groud_id(), net_pack.get_run_id(), "");
			send_message(net_pack, serialized_data.data(), serialized_data.size());
		}
		else if (net_pack.get_type() == NetPackage::TERMINATE)
		{
			terminate = true;
		}
		else 
		{
			cout << "received unsupported messaged type: " << net_pack.get_type() << endl;
		}
	}
	//cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
	//char buf[256];
    //gets_s(buf, sizeof(buf));
}