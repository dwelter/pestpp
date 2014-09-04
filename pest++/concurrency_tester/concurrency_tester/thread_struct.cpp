#include <mutex>
#include "thread_structs.h"

thread_flag::thread_flag(bool _flag)
{
	flag = _flag;	
}

bool thread_flag::set(bool f)
{
	std::lock_guard<std::mutex> lock(m);
	flag = f;
	return flag;
}
bool thread_flag::get()
{
	std::lock_guard<std::mutex> lock(m);
	if (flag) return true;
	else return false;
}