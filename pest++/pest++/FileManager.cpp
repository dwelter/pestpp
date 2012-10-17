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
#include "FileManager.h"
#include "utilities.h"
#include "OutputFileWriter.h"

using namespace pest_utils;

FileManager::FileManager(const string &_base_filename, const string &_directory)
	: pest_base_filename(strip_cp(_base_filename)), directory(strip_cp(_directory))
{
	ofstream &f_rec = open_file("rec");
	ofstream &f_sen = open_file("sen");
	OutputFileWriter::write_sen_header(f_sen, pest_base_filename);
}


FileManager::~FileManager(void)
{
	for (auto &i : file_map)
	{
		if(i.second.is_open()) i.second.close();
	}
}

string FileManager::build_filename(const string &ext)
{
	return directory + "\\" + pest_base_filename +"." + strip_cp(ext);
}

ofstream &FileManager::open_file(const string &extension)
{
	string filename = build_filename(extension);
	pair<map<string ,ofstream>::iterator, bool> ret;
	ret = file_map.insert(pair<string, ofstream>(extension, ofstream()));
	ofstream &f_new = ret.first->second;
	if (ret.second != false  && !f_new.is_open())
	{
		f_new.open(filename);
	}
	if (!f_new.good())
	{
		throw PestFileError(filename);
	}
	return f_new;
}

ofstream &FileManager::get_ofstream(const string &extension)
{
	map<string, ofstream>::iterator it;
	it = file_map.find(extension);
	// need to add error checking
	if (it == file_map.end())
	{
		string filename = build_filename(extension);
		throw PestFileErrorAccess(filename);
	}
	return it->second;
}

ofstream &FileManager::rec_ofstream()
{
	return get_ofstream("rec");
}

ofstream &FileManager::sen_ofstream()
{
	return get_ofstream("sen");
}