#include <mutex>

class thread_flag
{
public:
	
	thread_flag(bool _flag);
	bool set(bool _flag);
	bool get();
	
private:
	bool flag;
	std::mutex m;
	
};
