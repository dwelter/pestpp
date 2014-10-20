#include "network_wrapper.h"
#include "network_package.h"
#include "utilities.h"
#include "system_variables.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include <thread>

#ifdef OS_WIN
#include <Windows.h>
#include <conio.h>
#endif


#ifdef OS_LINUX
#include <arpa/inet.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/select.h>
#include<sys/wait.h>
#include <errno.h>
#endif

using namespace std;

void w_init()
{
#ifdef OS_WIN
	WSADATA wsaData;

	if (WSAStartup(MAKEWORD(2,0),  &wsaData)!=0)
	{
		cerr << "WSAStartup failed " << w_get_error_msg() << endl;
	}
#endif
}

int w_close(int sockfd) 
{
  int n;
	#ifdef OS_WIN
	shutdown(sockfd, SD_BOTH);
	if ((n = closesocket(sockfd)) != 0)
	{
		cerr << "error closing socket: "  << w_get_error_msg() << endl;  
	}
	return n;
	#endif
	#ifdef OS_LINUX
	n = shutdown(sockfd, 2);
	#endif
	return n;
}
void w_cleanup()
{
   #ifdef OS_WIN
	WSACleanup();
   #endif
}


int w_getaddrinfo(const char *node, const char *service,
			  const struct addrinfo *hints, struct addrinfo **res)
{
	int status;
	if ((status = getaddrinfo(node, service, hints, res)) !=0)
	{
		cerr << "getaddrinfo error: " << gai_strerror(status) << endl;
	}
	return status;
}

vector<string> w_getnameinfo_vec(int sockfd, int flags)
{
	int err;
	vector<string> name_info;
	char host[INET6_ADDRSTRLEN];
	char port[INET6_ADDRSTRLEN];
	struct sockaddr_storage addr;
	socklen_t addr_len = sizeof addr;
	err = getpeername(sockfd, (struct sockaddr*) &addr, &addr_len);
	err = getnameinfo((struct sockaddr*) &addr, addr_len, host, sizeof host, port, sizeof port, flags); 
	name_info.push_back(host);
	name_info.push_back(port);
	return name_info;
}

int w_socket(int domain, int type, int protocol)
{
	int sockfd = socket(domain, type, protocol);
	if (sockfd < 0) {
		cerr << "socket error: "  << w_get_error_msg() << endl;
	}
	return sockfd;
}

int w_connect(int sockfd, struct sockaddr *serv_addr, socklen_t addrlen)
{
	int n=0;
	if ((n=connect(sockfd, serv_addr, addrlen)) == -1 )
	{
		cerr << "connect error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_bind(int sockfd, struct sockaddr *my_addr, socklen_t addrlen)
{
	int n=0;
	if ((n=::bind(sockfd, my_addr, addrlen)) == -1 )
	{
		
		cerr << "bind error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_accept(int sockfd, struct sockaddr *addr, socklen_t *addrlen)
{
	int n=0;
	if ((n=accept(sockfd, addr, addrlen)) == -1) 
	{
		cerr << "bind error: " << w_get_error_msg() << endl;
	}
	return n;
}


int w_listen(int sockfd, int backlog)
{
	int n;
	if ((n = listen(sockfd, backlog)) == -1)
	{
		cerr << "listen error: "  << w_get_error_msg() << endl;
	}
	return n;
}

int w_recv(int sockfd, char *buf, size_t len, int flags)
{
	int n;
	n = recv(sockfd, buf, len, flags);
	if (n < 0){
		cerr << "recv error: "  << w_get_error_msg() << endl;
	}
	return n;
}
int w_send(int sockfd, char *buf, size_t len, int flags)
{
	int n;
	n = send(sockfd, buf, len, flags);
	if (n < 0){
		cerr << "send error: "  << w_get_error_msg() << endl;
	}
	return n;
}

int w_sendall(int sockfd, char *buf, unsigned long *len)
{
	unsigned long total = 0; // how many bytes we've sent
	unsigned long bytesleft = *len; // how many we have left to send
	int n;
	while(total < *len) {
		n = send(sockfd, buf+total, bytesleft, 0);
		if (n == -1) { break; }  //error
		if (n == 0) { break; } //connection closed
		total += n;
		bytesleft -= n;
	}
	*len = total; // return number actually sent here
		if (n < 0){
		cerr << "w_sendall error: " << n << endl;
	}
	if (n > 0) {n = 1;}
	return n; // return -1 on failure, 0 closed connection or 1 on success
}


int w_recvall(int sockfd, char *buf, unsigned long *len)
{
	unsigned long total = 0; // how many bytes we've received
	unsigned long bytesleft = *len; // how many we have left to receive
	int n;
	while(total < *len) {
		n = recv(sockfd, buf+total, bytesleft, 0);
		if (n == -1) { break; }  //error
		if (n == 0) { break; } //connection closed
		total += n;
		bytesleft -= n;
	}
	*len = total; // return number actually received here
	if (n < 0){
		//cerr << "w_recvall error: " << n << endl;
	}
	if (n > 0) {n = 1;}
	return n; // return -1 on failure, 0 closed connection or 1 on success
}

addrinfo* w_bind_first_avl(addrinfo *servinfo, int &sockfd)
{
	// loop through all the results and bind to the first we can
	struct addrinfo *p;
	char yes = '1';
	for(p = servinfo; p != nullptr; p = p->ai_next) 
	{
		if ((sockfd = w_socket(p->ai_family, p->ai_socktype,
		p->ai_protocol)) == -1) 
		{
			continue;
		}
		if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes,
		sizeof(int)) == -1)
		{
			return nullptr;
		}
		if (w_bind(sockfd, p->ai_addr, p->ai_addrlen) == -1) 
		{
			w_close(sockfd);
			continue;
		}
	break;
	}
	if (p == nullptr) {
		cerr << "server: failed to bind" << endl;
	}
	return p;
}

addrinfo* w_connect_first_avl(addrinfo *servinfo, int &sockfd)
{
	// loop through all the results and connect to the first we can
	struct addrinfo *p;
	for(p = servinfo; p != nullptr; p = p->ai_next) 
	{
		if ((sockfd = w_socket(p->ai_family, p->ai_socktype,
		p->ai_protocol)) == -1) 
		{
			continue;
		}
		if (w_connect(sockfd, p->ai_addr, p->ai_addrlen) == -1) 
		{
			w_close(sockfd);
			continue;
		}
	break;
	}
	if (p == nullptr) {
		// connection failed
	}
	return p;
}



void w_print_servinfo(addrinfo *res, ostream &fout)
{
	struct addrinfo *p;
	fout << "IP addresses:" << endl;
	for(p = res;p != NULL; p = p->ai_next) 
	{
		string socket_string = w_get_addrinfo_string(p);
	fout << "  " << socket_string << endl;
	}
}

string w_get_addrinfo_string(struct addrinfo *p)
{
	stringstream sstr;
	void *addr;
	char ipstr[INET6_ADDRSTRLEN];
	string ipver;
	unsigned short port = 0;
	// get the pointer to the address itself,
	// different fields in IPv4 and IPv6:
	if (p->ai_family == AF_INET) { // IPv4
		struct sockaddr_in *ipv4 = (struct sockaddr_in *)p->ai_addr;
		addr = &(ipv4->sin_addr);
		port = ntohs(ipv4->sin_port);
		ipver = "IPv4";
		inet_ntop(p->ai_family, addr, ipstr, sizeof ipstr);
		sstr << ipstr <<":" << port<< " (" << ipver << ")";
	}
	else { // IPv6
		struct sockaddr_in6 *ipv6 = (struct sockaddr_in6 *)p->ai_addr;
		addr = &(ipv6->sin6_addr);
		port = ntohs(ipv6->sin6_port);
		ipver = "IPv6";
		inet_ntop(p->ai_family, addr, ipstr, sizeof ipstr);
		sstr << "[" << ipstr <<"]:" << port << " (" << ipver << ")";
	}
	// convert the IP to a string and print it:
	return sstr.str();
}


int w_select(int numfds, fd_set *readfds, fd_set *writefds,
		   fd_set *exceptfds, struct timeval *timeout)
{
	int n;
	if ((n=select(numfds, readfds, writefds, exceptfds, timeout)) == -1)
	{
		cerr << "select error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_memcpy_s(void *dest, size_t numberOfElements, const void *src, size_t count)
{
	int err = 0;
		#ifdef OS_WIN
	err = memcpy_s(dest, numberOfElements, src, count);
		#endif
		#ifdef OS_LINUX
	memcpy(dest, src, count);
		#endif
	if (err) {
		cerr << "Error executing memcpy" << endl;
	}
	return err;
}

string w_get_error_msg()
{
	stringstream err_msg;
	#ifdef OS_WIN
	int err_no = WSAGetLastError();
	LPVOID lpMsgBuf;
	FormatMessage(
		FORMAT_MESSAGE_ALLOCATE_BUFFER |
		FORMAT_MESSAGE_FROM_SYSTEM |
		FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		err_no,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
		(LPTSTR) &lpMsgBuf,
		0,
		NULL);
	// Process any inserts in lpMsgBuf.
	// ..
	err_msg << "(tcp/ip error id = " << errno << ") " << (LPCTSTR)lpMsgBuf ;
	// Free the buffer.
	LocalFree( lpMsgBuf );
	#endif
	#ifdef OS_LINUX
	err_msg << "(tcp/ip error id = " << errno << ") " << strerror(errno);
	#endif
	return err_msg.str();
}

void w_sleep(int millisec)
{
   #ifdef OS_WIN
	Sleep(millisec);
   #endif
   #ifdef OS_LINUX
	sleep(millisec / 1000);
   #endif
}


#ifdef OS_WIN
PROCESS_INFORMATION start_command(string &cmd_string)
{
	char* cmd_line = _strdup(cmd_string.c_str());
	STARTUPINFO si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));
	if (!CreateProcess(NULL, cmd_line, NULL, NULL, false, 0, NULL, NULL, &si, &pi))
	{
		std::string cmd_string(cmd_line);
		throw std::runtime_error("CreateProcess() failed for command: " + cmd_string);
	}
	return pi;
}
#endif


#ifdef OS_LINUX
int start_command(string &cmd_string)
{
	//split cmd_string on whitespaces
	stringstream cmd_ss(cmd_string);
	string cmd;
	vector<string> cmds;
	while (getline(cmd_ss,cmd))
	{
		cmds.push_back(cmd);
	}
	//create the strurture for execv
	vector<char const*> arg_v;
	for (size_t icmd=0; icmd<cmds.size(); ++icmd)
	  {
	    arg_v.push_back(cmds[icmd].data());
	  }
	//char * const*argv = new char* [cmds.size()+1];
	//for (size_t icmd=0; icmd<cmds.size(); ++icmd)
	//{
	//  argv[icmd] = cmds[icmd].data();
	//}
	//argv[cmds.size() + 1] = NULL; //last arg must be NULL
		arg_v.push_back(NULL);
	      	
	pid_t pid = fork();
	if (pid == 0)
	{
	  int success = execv(arg_v[0], const_cast<char* const*>(&(arg_v[0])));
	     if (success == -1)
	     {
	       throw std::runtime_error("execv() failed for command: " + cmd_string);
	     }
	}
	return pid;
}

#endif


void w_run_commands(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished, vector<string> commands)
{
#ifdef OS_WIN
	//a flag to track if the run was terminated
	bool term_break = false;
	//create a job object to track child and grandchild process
	HANDLE job = CreateJobObject(NULL, NULL);
	if (job == NULL) throw PestError("could not create job object handle");
	JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = { 0 };
	jeli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
	if (0 == SetInformationJobObject(job, JobObjectExtendedLimitInformation, &jeli, sizeof(jeli)))
	{
		throw PestError("could not assign job limit flag to job object");
	}
	for (auto &cmd_string : commands)
	{		
		//start the command
		PROCESS_INFORMATION pi;
		try
		{
			pi = start_command(cmd_string);
		}
		catch (...)
		{
			finished->set(true);
			throw std::runtime_error("start_command() failed for command: " + cmd_string);
		}
		if (0 == AssignProcessToJobObject(job, pi.hProcess))
		{
			throw PestError("could not add process to job object: " + cmd_string);
		}
		DWORD exitcode;
		while (true)
		{
			//sleep			
			std::this_thread::sleep_for(std::chrono::milliseconds(OperSys::thread_sleep_milli_secs));
			//check if process is still active
			GetExitCodeProcess(pi.hProcess, &exitcode);
			//if the process ended, break
			if (exitcode != STILL_ACTIVE)
			{
				break;
			}
			//check for termination flag
			if (terminate->get())
			{
				std::cout << "recieved terminate signal" << std::endl;
				//try to kill the process
				bool success = CloseHandle(job);
				//bool success = TerminateProcess(pi.hProcess, 0);
				if (!success)
				{
					finished->set(true);
					throw std::runtime_error("unable to terminate process for command: " + cmd_string);
				}
				term_break = true;
				break;
			}
		}
		//jump out of the for loop if terminated
		if (term_break) break;
	}
	//set the finished flag for the listener thread
	finished->set(true);
	return;

#endif

#ifdef OS_LINUX
	//a flag to track if the run was terminated
	bool term_break = false;
	for (auto &cmd_string : commands)
	{		
		//start the command
		int command_pid = start_command(cmd_string);
		while (true)
		{
			//sleep
			std::this_thread::sleep_for(std::chrono::milliseconds(OperSys::thread_sleep_milli_secs));
			//check if process is still active
			int status;
			pid_t exit_code = waitpid(command_pid, &status, WNOHANG);
			//if the process ended, break
			if (exit_code == -1)
			{
				finished->set(true);
				throw std::runtime_error("waitpid() returned error status for command: " + cmd_string);
			}
			else if (exit_code != 0)
			{
				break;
			}
			//check for termination flag
			if (terminate->get())
			{
				std::cout << "recieved terminate signal" << std::endl;
				//try to kill the process
				errno = 0;
				int success = kill(command_pid, SIGKILL);
				if (success == -1)
				{
					finished->set(true);
					throw std::runtime_error("unable to terminate process for command: " + cmd_string);
				}
				term_break = true;
				break;
			}
		}
		//jump out of the for loop if terminated
		if (term_break) break;
	}
	//set the finished flag for the listener thread
	finished->set(true);
	return;
#endif

}

