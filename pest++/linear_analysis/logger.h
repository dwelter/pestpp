#ifndef LOGGER_H_
#define LOGGER_H_

#include <iostream>
#include <fstream>
#include <chrono>
#include <map>

class Logger
{
public:
	Logger(bool _echo = false){ echo = _echo; }
	Logger(std::ofstream &_fout,bool _echo=false);
	void log(const std::string &message);
	void write(const std::string &message);
	void error(const std::string &message);
	void warning(const std::string &message);
	~Logger();
private:
	bool echo;
	std::ofstream* fout;
	std::chrono::system_clock::time_point prev_time;	
	std::map<std::string, std::chrono::system_clock::time_point> tagged_events;
	std::string time_to_string(const std::chrono::system_clock::time_point &tmp_time);
	std::string elapsed_time_to_string(std::chrono::system_clock::time_point &current_time, std::chrono::system_clock::time_point &prev_time);
	void writetime(std::stringstream &os, time_t tc);
};

#endif 
