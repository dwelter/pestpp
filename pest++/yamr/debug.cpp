#include "debug.h"

using namespace std;

std::ofstream fout_dbg;

std::ostream& operator<< (std::ostream &os, const std::set<std::string> val)
{
	for (const auto &i : val)
	{
		os << i << endl;
	}
	return os;
}

std::ostream& operator<< (std::ostream &os, const std::vector<std::string> val)
{
	for (const auto &i : val)
	{
		os << i << endl;
	}
	return os;
}