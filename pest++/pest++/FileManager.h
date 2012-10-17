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
#include <map>

using namespace std;

class FileManager
{
public:
	FileManager(const string &_base_filename, const std::string &_directory="");
	std::string ctl_filename() {return build_filename("pst");}
	std::string build_filename(const std::string &ext);
	std::string jacobian_filename() {return build_filename("jco");}
	std::string iteration_jacobian_filename() {return build_filename("jci");}
	void set_analytic_derivative_filename(const std::string &name) {analytic_derivative_filename = name;}
	std::ofstream &rec_ofstream();
	std::ofstream &sen_ofstream();
	std::ofstream &open_ofile_ext(const std::string &extension, ios_base::openmode mode = ofstream::out);
	std::ofstream &open_ofile_local(const std::string &tag, const std::string &filename, ofstream::openmode mode = ofstream::out);
	std::ofstream &open_ofile_absolute(const std::string &tag, const std::string &filename, ofstream::openmode mode = ofstream::out);
	std::ifstream &open_ifile_ext(const std::string &extension, ifstream::openmode mode = ifstream::in);
	std::ifstream &open_ifile_local(const std::string &tag, const std::string &filename, ifstream::openmode mode = ifstream::in);
	std::ifstream &open_ifile_absolute(const std::string &tag, const std::string &filename, ifstream::openmode mode = ifstream::in);
	std::fstream &open_iofile_ext(const std::string &extension, ios_base::openmode mode = fstream::in | fstream::out);
	std::fstream &open_iofile_local(const std::string &tag, const std::string &filename, fstream::openmode mode = fstream::in | fstream::out);
	std::fstream &open_iofile_absolute(const std::string &tag, const std::string &filename, fstream::openmode mode = fstream::in | fstream::out);
	void close_file(const std::string &extension);
	std::ofstream &get_ofstream(const std::string &tag);
	std::ifstream &get_ifstream(const std::string &tag);
	std::fstream &get_fstream(const std::string &tag);
	~FileManager(void);
private:
	std::string analytic_derivative_filename;
	std::string directory;
	std::string pest_base_filename;
	std::map<string, std::ofstream> ofile_map;
	std::map<string, std::ifstream> ifile_map;
	std::map<string, std::fstream> iofile_map;
	std::map<std::string, std::string> filename_map;
};

#endif /* FILEMANAGER_H_ */
