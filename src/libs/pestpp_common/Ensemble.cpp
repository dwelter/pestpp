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

Mat Ensemble::to_matrix()
{
	return to_matrix(real_names, var_names);
}

Mat Ensemble::to_matrix(vector<string> &row_names, vector<string> &col_names)
{
	return Mat();
}


void Ensemble::to_csv(string &file_name)
{

}

Eigen::VectorXd Ensemble::get_real_vector(int ireal)
{
	if (ireal >= shape().first)
	{
		stringstream ss;
		ss << "Ensemble::get_real_vector() : ireal (" << ireal << ") >= reals.shape[0] (" << ireal << ")";
		throw_ensemble_error(ss.str());
	}
	return reals.row(ireal);
}

Eigen::VectorXd Ensemble::get_real_vector(const string &real_name)
{

	int idx = find(real_names.begin(), real_names.end(), real_name) - real_names.begin();
	if (idx >= real_names.size())
	{
		stringstream ss;
		ss << "Ensemble::get_real_vector() real_name '" << real_name << "' not found";
		throw_ensemble_error(ss.str());
	}
	return get_real_vector(idx);
}

void Ensemble::throw_ensemble_error(string &message)
{
	string full_message = "Ensemble Error: " + message;
	performance_log->log_event(full_message);
	performance_log->~PerformanceLog();
	file_manager.rec_ofstream() << endl << endl << "***" << full_message << endl;
	file_manager.~FileManager();
	throw runtime_error(full_message);
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

void Ensemble::from_csv(int num_reals,ifstream &csv)
{
	int lcount = 0;
	//vector<vector<double>> vectors;
	double val;
	string line;
	vector<string> tokens;
	//vector<double> vals;
	string real_id;

	real_names.clear();
	reals.resize(num_reals, var_names.size());
	int irow = 0;
	while (getline(csv, line))
	{
		pest_utils::strip_ip(line);
		tokens.clear();
		//vals.clear();
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

		for (int jcol = 1; jcol<tokens.size(); ++jcol)
		{
			try
			{
				val = pest_utils::convert_cp<double>(tokens[jcol]);
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << "error converting token '" << tokens[jcol] << "' to double for " << var_names[jcol - 1] << " on line " << lcount << " : " << e.what();
				throw runtime_error(ss.str());
			}
			//vals.push_back(val);
			reals(irow, jcol - 1) = val;
		}
		//vectors.push_back(vals);
		lcount++;

	}
	//Eigen::MatrixXd reals(vectors.size(), vals.size());
	if (lcount != num_reals)
		throw runtime_error("different number of reals found");

	/*reals.resize(vectors.size(), vals.size());
	for (int i = 0; i < vectors.size(); i++)
	{
		double* ptr = &vectors[i][0];
		Eigen::Map<Eigen::VectorXd> vec_map(ptr, vectors[i].size());
		reals.row(i) = vec_map;
	}*/

}


void Ensemble::initialize_with_csv(string &file_name)
{
	
	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening parameter csv " + file_name + " for reading");


	var_names = prepare_csv(pest_scenario.get_ctl_ordered_par_names(),csv,false);
	//blast through the file to get number of reals
	string line;
	int num_reals = 0;
	while (getline(csv, line))
		num_reals++;
	csv.close();
	
	csv.open(file_name);
	if (!csv.good())
		throw runtime_error("error re-opening parameter csv " + file_name + " for reading");
	getline(csv, line);
	from_csv(num_reals, csv);

}

ParameterEnsemble::ParameterEnsemble(const ParamTransformSeq &_par_transform, Pest &_pest_scenario,
	FileManager &_file_manager, OutputFileWriter &_output_file_writer,
	PerformanceLog *_performance_log, unsigned int seed):
	Ensemble(_pest_scenario,_file_manager,_output_file_writer,_performance_log,seed),par_transform(_par_transform)
{

}

void ParameterEnsemble::enforce_bounds()
{

}

void ParameterEnsemble::to_csv(string &file_name)
{
	Ensemble::to_csv(file_name);
}

ObservationEnsemble::ObservationEnsemble(ObjectiveFunc *_obj_func, Pest &_pest_scenario,
	FileManager &_file_manager, OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, unsigned int seed) :
	Ensemble(_pest_scenario, _file_manager, _output_file_writer, _performance_log, seed), obj_func_ptr(_obj_func)
{

}

EnsemblePair::EnsemblePair(ParameterEnsemble &_pe, ObservationEnsemble &_oe) : pe(_pe), oe(_oe)
{
	pair<int, int> par_shape = pe.shape();
	pair<int, int> obs_shape = oe.shape();
	if (par_shape.first != obs_shape.first)
	{
		stringstream ss;
		ss << "parameter ensemble shape[0] (" << par_shape.first << ") != observation ensemble shape[0] (" << obs_shape.first << ")";
		pe.throw_ensemble_error(ss.str());
	}
}

void EnsemblePair::run(RunManagerAbstract *run_mgr_ptr)
{
	
	Parameters pars = pe.get_pest_scenario_ptr()->get_ctl_parameters();
	//for (int ireal = 0; ireal < pe.shape().first; ireal++)
	Eigen::VectorXd evec;
	vector<double> svec;
	vector<string> var_names = pe.get_var_names();
	for (auto &real_name : pe.get_real_names())
	{
		evec = pe.get_real_vector(real_name);
		const vector<double> svec(evec.data(), evec.data() + evec.size());
		pars.update_without_clear(var_names,svec);

		//run_mgr_ptr->add_run()
	}
}