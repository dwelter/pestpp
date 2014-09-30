/*
© Copyright 2012, David Welter

This file is part of PEST++.

PEST++ is free software : you can redistribute it and / or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PEST++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PEST++.If not, see<http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include "RunStorage.h"
#include "utilities.h"
using namespace std;


void usage(ostream &fout)
{
	fout << "--------------------------------------------------------" << endl;
	fout << "usage:" << endl << endl;
	fout << "    ascii2pbin.exe bin_file  results_file" << endl << endl;
	fout << "    where:" << endl;
	fout << "          ext_file:   name of the pest++ .ext file containing the name binary run file" << endl;
	fout << "        results_file:   name of file which conatains the list of completed runs" << endl;
	fout << "--------------------------------------------------------" << endl;
}

void read_data_file(const string &filename, Transformable &data)
{
	ifstream fin;
	fin.open(filename);
	string line;
	string name;
	double value;
	data.clear();
	while (getline(fin, line))
	{
		vector<string> tokens;
		pest_utils::strip_ip(line);
		if (!line.empty() && line[0] != '#')
		{
			pest_utils::tokenize(line, tokens);
			name = pest_utils::upper_cp(tokens[0]);
			pest_utils::convert_ip(tokens[1], value);
			data.insert(name, value);
		}
	}
	fin.close();
}

//int run_name_2_run_id(const string &run_name)
//{
//	string str = pest_utils::strip_cp(run_name);
//	vector<string> tokens;
//	pest_utils::tokenize(str, tokens, "-");
//	int id = pest_utils::convert_cp<int>(tokens[1]);
//	return id;
//}

int main(int argc, char* argv[])
{
	string ext_filename = argv[1];
	string results_filename = argv[2];

	// Error checking
	if (argc != 3)
	{
		cerr << "Error: incorect number of command line arguements" << endl << endl;
		usage(cerr);
		cerr.flush();
		throw PestError("Error: incorect number of command line arguements");
	}

	string pbin_filename;
	int max_n_fail = 1;
	try
	{
		ifstream fin_ext(ext_filename);
		getline(fin_ext, pbin_filename);
		pest_utils::strip_ip(pbin_filename);
		string tmp_str;
		getline(fin_ext, tmp_str);
		pest_utils::strip_ip(tmp_str);
		max_n_fail = pest_utils::convert_cp<int>(tmp_str);

	}
	catch (exception &e)
	{
		cerr << "Error processing external run manager *.ext \"" << pbin_filename << "\"";
		cerr << e.what() << endl;
		usage(cerr);
		throw(e);
	}
	// Open PEST++ binary storage file
	RunStorage rs("");
	try
	{
		rs.init_restart(pbin_filename);
	}
	catch (exception &e)
	{
		cerr << "Error opening PEST++ binary file \"" << pbin_filename << "\"";
		cerr << e.what() << endl;
		usage(cerr);
		cerr.flush();
		throw(e);
	}

	//Open resutls file
	ifstream fin_results;
	try
	{
		fin_results.open(results_filename);
	}
	catch (exception &e)
	{
		cerr << "Error opening results file \"" << results_filename << "\"";
		cerr << e.what() << endl;
		usage(cerr);
		cerr.flush();
		throw(e);
	}

	// Finished with initial error checking.  Start the main program 
	int status = 0;
	string info_text;
	double info_value;
	vector<double> pars_vec;
	vector<double> obs_vec;
	vector<string> par_name_vec = rs.get_par_name_vec();
	vector<string> obs_name_vec = rs.get_obs_name_vec();
	int n_runs = rs.get_nruns();
	cout << "processing runs" << endl;

	//set the status of all outstanding runs to failed
	for (int id = 0; id < n_runs; ++id)
	{
		rs.get_info(id, status, info_text, info_value);
		if (status == 0)
		{
			rs.update_run_failed(id);
			rs.set_run_nfailed(id, max_n_fail);
		}
	}

	// Process each line in the results file
	string line;
	while (getline(fin_results, line))
	{
		try
		{
			vector<string> tokens;
			pest_utils::strip_ip(line);
			if (!line.empty() && line[0] != '#')
			{
				pest_utils::tokenize(line, tokens);
				int run_id = pest_utils::convert_cp<int>(tokens[0]);
				string out_filename = tokens[1];
				//read updated observations values
				Observations obs;
				read_data_file(out_filename, obs);
				if (tokens.size() > 2)
				{
					string par_filename = tokens[2];
					//read updated parameter values
					Parameters pars;
					read_data_file(par_filename, pars);
					//update PEST++ binary run file
					rs.update_run(run_id, pars, obs);
				}
				else
				{
					rs.update_run(run_id, obs);
				}
			}
		}
		catch (exception &e)
		{
			cerr << "Error processing line: \"" << line << endl;
			cerr << "The model run associated with this line will be marked as a failed run" << endl;
			cerr << e.what() << endl;
			cerr << "Error processing line: \"" << line << endl;
			cerr << endl;
			cerr.flush();
		}
		fin_results.close();
	}
}