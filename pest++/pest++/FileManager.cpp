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
	string rec_filename = build_filename("rec");
	f_rec.open(rec_filename.c_str());
	string sen_filename = build_filename("sen");
	f_sen.open(sen_filename.c_str());
	OutputFileWriter::write_sen_header(f_sen, pest_base_filename);
}


FileManager::~FileManager(void)
{
	f_rec.close();
	f_sen.close();
}

string FileManager::build_filename(const string &ext)
{
	return directory + "\\" + pest_base_filename +"." + strip_cp(ext);
}