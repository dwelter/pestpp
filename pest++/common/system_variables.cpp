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

#include <string>
#include <sstream>
#include <cmath>
#include "system_variables.h"

#ifdef OS_WIN
 #include <direct.h>
#endif

#ifdef OS_LINUX
 #include <unistd.h>
#endif
using namespace std;

#ifdef OS_WIN
const std::string OperSys::DIR_SEP = "\\";
const std::string OperSys::COMMAND_LINE_APPEND = " ; ";
#endif
#ifdef OS_LINUX
const std::string OperSys::DIR_SEP = "/";
const std::string OperSys::COMMAND_LINE_APPEND = " & ";
#endif

void OperSys::string2pathname(string &s)
{
	size_t i;
	size_t len(s.size());
	string dir_chars("/\\");
	stringstream new_s;

	for (i=0; i<len; ++i) {
		if (s[i] == '/') {new_s << DIR_SEP;}
		else {new_s << s[i];}
		if (i!=len-1 && (s[i] == '\\' || s[i]=='/')) {
			new_s << DIR_SEP;
		}
	}
	s = new_s.str();
}

string OperSys::getcwd()
{
    #ifdef OS_WIN
	char *buffer;
	buffer = _getcwd( NULL, 0 );
	string cwd(buffer);
    free(buffer);
	return cwd;
    #endif
    #ifdef OS_LINUX
        char *buffer;
	buffer = ::getcwd( NULL, 0 );
	string cwd(buffer);
        free(buffer);
	return cwd;
    #endif
}

void OperSys::chdir(const char *str)
{
   #ifdef OS_WIN
      _chdir(str);
   #endif
   #ifdef OS_LINUX
      //chdir(str);
   #endif
}

char* OperSys::gets_s(char *str, size_t len)
{
 #ifdef OS_WIN
  return ::gets_s(str, len);
 #endif
 #ifdef OS_LINUX
  return gets(str);
 #endif

}

bool OperSys::double_is_invalid(double x)
{
#ifdef OS_WIN
   bool test = (_isnan(x)!=0 || _finite(x) == 0);
   return test;
#endif
#ifdef OS_LINUX
   bool test = (std::isnan(x) || std::isinf(x));
  return test;
#endif
}
