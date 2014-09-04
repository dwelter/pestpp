// exe_tester.cpp : Defines the entry point for the console application.
//


#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <mutex>
#include <chrono>
#include <vector>
#include <exception>
#include "thread_structs.h"

int g_sleeptime = 1000; //milliseconds

//windows
#ifdef OS_WIN
#include <Windows.h>
#include <conio.h>

//starts a process and sends back handles
PROCESS_INFORMATION start_command(char* &cmd_line)
{
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

//runs the commands while checking for terminate signal
void runner(thread_flag* terminate, thread_flag* finished, std::vector<std::string> commands)
{
	//a flag to track if the run was terminated
	bool term_break = false;
	for (auto &cmd_string : commands)
	{
		char* cmd_line = _strdup(cmd_string.c_str());
		//start the command
		PROCESS_INFORMATION pi;
		try
		{
			pi = start_command(cmd_line);
		}
		catch (...)
		{
			finished->set(true);
			throw std::runtime_error("start_command() failed for command: " + cmd_string);
		}
		DWORD exitcode;
		while (true)
		{
			//sleep
			std::this_thread::sleep_for(std::chrono::milliseconds(g_sleeptime));
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
				bool success = TerminateProcess(pi.hProcess, 0);
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
}

//listens for keyboard input to send terminate signal
void listener(thread_flag* terminate, thread_flag* finished)
{
	while (true)
	{
		//sleep
		std::this_thread::sleep_for(std::chrono::milliseconds(g_sleeptime));
		//check for a key stroke
		if (_kbhit() != 0)
		{
			std::cout << "sending terminate signal" << std::endl;
			terminate->set(true);
			break;
		}
		//check if the runner thread has finished
		if (finished->get())
		{
			std::cout << "received finished signal " << std::endl;
			break;
		}
	}
	return;
}
#endif


//linux
#ifdef OS_LINIX
#include <sys\types.h>
#include <sys\select.h>
#include <errno.h>

int start_command(char* &cmd_line)
{
	pid_t pid = fork();
	if (pid == 0)
	{
		int success = execl(cmd_line);
		if (sucess == -1)
		{
			std::string cmd_string(cmd_line);
			throw std::runtime_error("execl() failed for command: " + cmd_string);
		}
	}
	else
	{
		return pid;
	}

}

void runner(thread_flag* terminate, thread_flag* finished, std::vector<std::string> commands)
{
	//a flag to track if the run was terminated
	bool term_break = false;
	for (auto &cmd_string : commands)
	{
		char* cmd_line = _strdup(cmd_string.c_str());
		//start the command
		int command_pid = start_command(cmd_line);
		int exit_code;
		while (true)
		{
			//sleep
			std::this_thread::sleep_for(std::chrono::milliseconds(g_sleeptime));
			//check if process is still active
			pid_t exit_code = waitpid(command_pid, &status, WNOHANG);
			//if the process ended, break
			if (exitcode == -1)
			{
				finished->set(true);
				throw std::runtime_error("waitpid() returned error status for command: " + cmd_string);
			}
			else if exitcode != 0)
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
}

int kbhit(void)
{
	struct timeval tv;
	fd_set read_fd;

	/* Do not wait at all, not even a microsecond */
	tv.tv_sec = 0;
	tv.tv_usec = 0;

	/* Must be done first to initialize read_fd */
	FD_ZERO(&read_fd);

	/* Makes select() ask if input is ready:
	* 0 is the file descriptor for stdin    */
	FD_SET(0, &read_fd);

	/* The first parameter is the number of the
	* largest file descriptor to check + 1. */
	if (select(1, &read_fd, NULL, /*No writes*/NULL, /*No exceptions*/&tv) == -1)
		return 0;  /* An error occured */

	/*  read_fd now holds a bit map of files that are
	* readable. We test the entry for the standard
	* input (file 0). */

	if (FD_ISSET(0, &read_fd))
		/* Character pending on stdin */
		return 1;

	/* no characters were pending */
	return 0;
}


void listener(thread_flag* terminate, thread_flag* finished)
{
	while (true)
	{
		//sleep
		std::this_thread::sleep_for(std::chrono::milliseconds(g_sleeptime));
		//check for a key stroke
		if (kbhit() != 0)
		{
			std::cout << "sending terminate signal" << std::endl;
			terminate->set(true);
			break;
		}
		//check if the runner thread has finished
		if (finished->get())
		{
			std::cout << "received finished signal " << std::endl;
			break;
		}
	}
	return;
}
#endif

int main(int argc, char* argv[])
{
	std::vector<std::string> commands;
	commands.push_back("python run.py");
	commands.push_back("python run.py");
	thread_flag terminate(false);
	thread_flag finished(false);
	std::thread runner_thread(runner, &terminate, &finished, commands);
	std::thread listener_thread(listener, &terminate, &finished);
	runner_thread.join();
	listener_thread.join();
	return 0;
}

