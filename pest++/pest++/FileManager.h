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
#ifndef FILEMANAGER_H_
#define FILEMANAGER_H_

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class FileManager
{
public:
	FileManager(const string &_base_filename, const string &_directory="");
	string ctl_filename() {return build_filename("pst");}
	string par_filename() {return build_filename("par");}
	string build_filename(const string &ext);
	string jacobian_filename() {return build_filename("jco");}
	string iteration_jacobian_filename() {return build_filename("jci");}
	const string &get_analytic_derivative_filename(){return analytic_derivative_filename;}
	void set_analytic_derivative_filename(const string &name) {analytic_derivative_filename = name;}
	ofstream &rec_ofstream () {return f_rec;}
	~FileManager(void);
private:
	string analytic_derivative_filename;
	string directory;
	string pest_base_filename;
	ofstream f_rec;
};

#endif /* FILEMANAGER_H_ */
