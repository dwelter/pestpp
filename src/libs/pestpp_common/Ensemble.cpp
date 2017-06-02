#include <random>
#include <iomanip>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"

mt19937_64 Ensemble::rand_engine = mt19937_64(1);


Ensemble::Ensemble(Pest &_pest_scenario,
	FileManager &_file_manager,OutputFileWriter &_output_file_writer,
	PerformanceLog *_performance_log,  unsigned int seed)
	: pest_scenario(_pest_scenario), file_manager(_file_manager),
	output_file_writer(_output_file_writer), performance_log(_performance_log)
{
	// initialize random number generator
	rand_engine.seed(seed);
}


Ensemble::~Ensemble()
{
}

vector<string> Ensemble::prepare_csv(vector<string> names, ifstream &csv, bool forgive)
{
	if (!csv.good())
	{
		throw runtime_error("ifstream not good");
	}

	//process the header
	//any missing header labels will be marked to ignore those columns later
	string line;
	vector<string> header_tokens;
	if (!getline(csv, line))
		throw runtime_error("error reading header (first) line from csv file :");
	pest_utils::strip_ip(line);
	pest_utils::upper_ip(line);
	pest_utils::tokenize(line, header_tokens, ",", false);
	//cout << tokens << endl;
	//vector<string> header_tokens = tokens;

	// check for parameter names that in the pest control file but that are missing from the csv file
	vector<string> missing_names;
	string name;
	for (auto &name : names)
		if (find(header_tokens.begin(), header_tokens.end(), name) == header_tokens.end())
			missing_names.push_back(name);

	if (missing_names.size() > 0)
	{
		stringstream ss;
		ss << " the following names were not found in the csv file header:" << endl;
		for (auto &n : missing_names) ss << n << endl;
		if (!forgive)
			throw runtime_error(ss.str());
		else
			cout << ss.str() << endl << "continuing anyway..." << endl;
	}

	//build up a list of idxs to use
	//vector<int> header_idxs;
	//map<string, int> header_info;
	vector<string> header_names;
	for (int i = 0; i < header_tokens.size(); i++)
	{
		if (find(names.begin(), names.end(), header_tokens[i]) != names.end())
		{
			//header_idxs.push_back(i);
			//header_info[header_tokens[i]] = i;
			header_names.push_back(header_tokens[i]);
		}
	}
	return header_names;

}

void Ensemble::from_csv(ifstream &csv)
{
	int lcount = 1;
	vector<vector<double>> vectors;
	double val;
	string line;
	vector<string> tokens;
	vector<double> vals;
	string real_id;

	real_names.clear();

	while (getline(csv, line))
	{
		pest_utils::strip_ip(line);
		tokens.clear();
		vals.clear();
		pest_utils::tokenize(line, tokens, ",", false);
		if (tokens.size() != var_names.size() + 1) // +1 for run id in first column
		{
			stringstream ss;
			ss << "error parsing csv file on line " << lcount << ": wrong number of tokens, ";
			ss << "expecting " << var_names.size() + 1 << ", found " << tokens.size();
			throw runtime_error(ss.str());
		}
		try
		{
			pest_utils::convert_ip(tokens[0], real_id);
		}
		catch (exception &e)
		{
			stringstream ss;
			ss << "error converting token '" << tokens[0] << "' to <string> run_id on line " << lcount << ": " << line << endl << e.what();
			throw runtime_error(ss.str());
		}
		real_names.push_back(real_id);

		for (int i = 1; i<tokens.size(); ++i)
		{
			try
			{
				val = pest_utils::convert_cp<double>(tokens[i]);
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << "error converting token '" << tokens[i] << "' to double for " << var_names[i - 1] << " on line " << lcount << " : " << e.what();
				throw runtime_error(ss.str());
			}
			vals.push_back(val);

		}
		vectors.push_back(vals);
		lcount++;

	}
	//Eigen::MatrixXd reals(vectors.size(), vals.size());

	reals.resize(vectors.size(), vals.size());
	for (int i = 0; i < vectors.size(); i++)
	{
		double* ptr = &vectors[i][0];
		Eigen::Map<Eigen::VectorXd> vec_map(ptr, vectors[i].size());
		reals.row(i) = vec_map;
	}

}
void Ensemble::initialize_with_csv(string &file_name)
{
	
	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening parameter csv " + file_name + " for reading");
	var_names = prepare_csv(pest_scenario.get_ctl_ordered_par_names(),csv,false);
	
}
