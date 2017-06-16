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

Eigen::MatrixXd Ensemble::get_eigen(vector<string> row_names, vector<string> col_names)
{
	if (row_names.size() == 0)
		throw_ensemble_error("Ensemble.get_eigen() row_names empty");
	if (col_names.size() == 0)
		throw_ensemble_error("Ensemble.get_eigen( col_names empty");
	vector<string> missing;
	vector<string>::iterator iter,start=real_names.begin(),end=real_names.end();
	vector<int> row_idxs,col_idxs;
	for (auto &rname : row_names)
	{
		iter = find(start, end, rname);
		if (iter == end)
			missing.push_back(rname);
		row_idxs.push_back(iter - start);
	}
	if (missing.size() > 0)
		throw_ensemble_error("Ensemble.get_eigen() the following row_names not found:", missing);
	start = var_names.begin();
	end = var_names.end();
	for (auto cname : col_names)
	{
		iter = find(start, end, cname);
		if (iter == end)
			missing.push_back(cname);
		col_idxs.push_back(iter - start);
	}
	if (missing.size() > 0)
		throw_ensemble_error("Ensemble.get_eigen() the following col_names not found:", missing);
	Eigen::MatrixXd mat;
	mat.resize(row_names.size(), col_names.size());
	int i=0, j=0;
	for (auto &ii : row_idxs)
	{
		j = 0;
		for (auto &jj : col_idxs)
		{
			//cout << i << ',' << j << ';' << ii << ',' << jj << endl;
			mat(i, j) = reals(ii, jj);
			j++;
		}
		i++;
	}
	return mat;
	
}

void Ensemble::to_csv(string &file_name)
{
	ofstream csv(file_name);
	if (!csv.good())
	{
		throw_ensemble_error("Ensemble.to_csv() error opening csv file " + file_name + " for writing");
	}
	for (auto &vname : var_names)
		csv << vname << ',';
	csv << endl;
	for (int ireal = 0; ireal < reals.rows(); ireal++)
	{
		csv << real_names[ireal] << ',';
		for (int ivar=0; ivar < reals.cols(); ivar++)
		{
			csv << reals.block(ireal, ivar,1,1) << ',';

		}
	}

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

void Ensemble::throw_ensemble_error(string message, vector<string> vec)
{
	stringstream ss;
	ss << ' ';
	for (auto &v : vec)
		ss << v << ',';
	throw_ensemble_error(message + ss.str());
}

void Ensemble::throw_ensemble_error(string message)
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
		enum transStatus { CTL, NUM, MODEL };
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
		irow++;

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




ParameterEnsemble::ParameterEnsemble(const ParamTransformSeq &_par_transform, Pest &_pest_scenario,
	FileManager &_file_manager, OutputFileWriter &_output_file_writer,
	PerformanceLog *_performance_log, unsigned int seed):
	Ensemble(_pest_scenario,_file_manager,_output_file_writer,_performance_log,seed),par_transform(_par_transform)
{

}

void ParameterEnsemble::initialize_with_csv(string &file_name)
{

	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening parameter csv " + file_name + " for reading");

	var_names = prepare_csv(pest_scenario.get_ctl_ordered_par_names(), csv, false);
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
	tstat = transStatus::CTL;
	//cout << reals << endl;
}

void ParameterEnsemble::enforce_bounds()
{
	throw_ensemble_error("pe.enforce_bounds() not implemented");
}

void ParameterEnsemble::to_csv(string &file_name)
{
	throw_ensemble_error("pe.to_csv() not implemented");
	if (tstat == ParameterEnsemble::transStatus::MODEL)
		

	Ensemble::to_csv(file_name);
}

ObservationEnsemble::ObservationEnsemble(ObjectiveFunc *_obj_func, Pest &_pest_scenario,
	FileManager &_file_manager, OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, unsigned int seed) :
	Ensemble(_pest_scenario, _file_manager, _output_file_writer, _performance_log, seed), obj_func_ptr(_obj_func)
{

}

void ObservationEnsemble::update_from_obs(int row_idx, Observations &obs)
{
	if (row_idx >= real_names.size())
		throw_ensemble_error("ObservtionEnsemble.update_from_obs() obs_idx out of range");
	reals.row(row_idx) = obs.get_data_eigen_vec(var_names);

}

void ObservationEnsemble::initialize_with_csv(string &file_name)
{

	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening observation csv " + file_name + " for reading");


	var_names = prepare_csv(pest_scenario.get_ctl_ordered_obs_names(), csv, false);
	//blast through the file to get number of reals
	string line;
	int num_reals = 0;
	while (getline(csv, line))
		num_reals++;
	csv.close();

	csv.open(file_name);
	if (!csv.good())
		throw runtime_error("error re-opening observation csv " + file_name + " for reading");
	getline(csv, line);
	from_csv(num_reals, csv);

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

	//initialize active real tracker
	for (int i = 0; i < par_shape.first; i++)
		active_real_indices.push_back(i);
}

void EnsemblePair::run(RunManagerAbstract *run_mgr_ptr)
{
	
	Parameters pars = pe.get_pest_scenario_ptr()->get_ctl_parameters();
	//for (int ireal = 0; ireal < pe.shape().first; ireal++)
	Eigen::VectorXd evec;
	vector<double> svec;
	vector<string> var_names = pe.get_var_names();
	//for (auto &real_name : pe.get_real_names())
	map<int, int> real_run_ids;
	int run_id;
	for (auto &act_idx : active_real_indices)
	{
		evec = pe.get_real_vector(act_idx);
		const vector<double> svec(evec.data(), evec.data() + evec.size());
		pars.update_without_clear(var_names,svec);
		if (pe.get_trans_status() == ParameterEnsemble::transStatus::CTL)
			pe.get_par_transform().ctl2model_ip(pars);
		else if (pe.get_trans_status() == ParameterEnsemble::transStatus::NUM)
			pe.get_par_transform().numeric2model_ip(pars);

		run_id = run_mgr_ptr->add_run(pars);
		real_run_ids[act_idx] =  run_id;
	}
	//pe.set_trans_status(ParameterEnsemble::transStatus::MODEL);
	//cout << pe.get_reals() << endl;
	
	run_mgr_ptr->run();

	set<int> failed_runs = run_mgr_ptr->get_failed_run_ids();
	vector<int> failed_real_idxs;
	Observations obs;
	cout << oe.get_reals() << endl;
	for (auto &real_run_id : real_run_ids)
	{
		if (failed_runs.find(real_run_id.second) != failed_runs.end())
			failed_real_idxs.push_back(real_run_id.first);
		else
		{
			run_mgr_ptr->get_run(real_run_id.second, pars, obs);
			oe.update_from_obs(real_run_id.first, obs);
		}
	}
	//cout << oe.get_reals() << endl;
	//cout << pe.get_reals() << endl;
}

Eigen::MatrixXd EnsemblePair::get_active_oe_eigen()
{
	return get_active_oe_eigen(oe.get_var_names());
}

Eigen::MatrixXd EnsemblePair::get_active_oe_eigen(const vector<string> &obs_names)
{
	vector<string> act_real_names, oe_real_names = oe.get_real_names();
	for (auto &i : active_real_indices)
		act_real_names.push_back(oe_real_names[i]);
	return pe.get_eigen(act_real_names, obs_names);
}

Eigen::MatrixXd EnsemblePair::get_active_pe_eigen()
{
	return get_active_pe_eigen(pe.get_var_names());
}

Eigen::MatrixXd EnsemblePair::get_active_pe_eigen(const vector<string> &par_names)
{
	vector<string> act_real_names, pe_real_names = pe.get_real_names();
	for (auto &i : active_real_indices)
		act_real_names.push_back(pe_real_names[i]);
	return pe.get_eigen(act_real_names, par_names);
}

