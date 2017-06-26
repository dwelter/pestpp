#include <random>
#include <iomanip>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "ParamTransformSeq.h"
#include "ObjectiveFunc.h"

mt19937_64 Ensemble::rand_engine = mt19937_64(1);

//Ensemble::Ensemble(Pest &_pest_scenario,
//	FileManager &_file_manager,OutputFileWriter &_output_file_writer,
//	PerformanceLog *_performance_log,  unsigned int seed)
//	: pest_scenario(_pest_scenario), file_manager(_file_manager),
//	output_file_writer(_output_file_writer), performance_log(_performance_log)
Ensemble::Ensemble(Pest *_pest_scenario_ptr): pest_scenario_ptr(_pest_scenario_ptr)
{
	// initialize random number generator
	rand_engine.seed(1123433458);
}

//Ensemble Ensemble::get(vector<string> &_real_names, vector<string> &_var_names)
//{
//	Ensemble new_en(pest_scenario_ptr);
//
//	if ((_real_names.size() == 0) && (_var_names.size() == 0))
//	{
//		new_en.from_eigen_mat(reals, real_names, var_names);
//	}
//	else
//	{
//		new_en.from
//	}
//
//}

Eigen::MatrixXd Ensemble::get_eigen_mean_diff()
{
	return get_eigen_mean_diff(vector<string>(),vector<string>());
}

Eigen::MatrixXd Ensemble::get_eigen_mean_diff(const vector<string> &_real_names, const vector<string> &_var_names)
{
	Eigen::MatrixXd _reals;
	if ((_real_names.size() == 0) && (_var_names.size() == 0))
		_reals = reals;
	else
		_reals = get_eigen(_real_names, _var_names);

	double mean;
	int s = _reals.rows();
	for (int j = 0; j < _reals.cols(); j++)
	{
		mean = _reals.col(j).mean();
		_reals.col(j) = _reals.col(j) - (Eigen::VectorXd::Ones(s) * mean);
	}
	return _reals;

}


void Ensemble::from_eigen_mat(Eigen::MatrixXd _reals, const vector<string> &_real_names, const vector<string> &_var_names)
{
	if (_reals.rows() != _real_names.size())
		throw_ensemble_error("Ensemble.from_eigen_mat() rows != real_names.size");
	if (_reals.cols() != _var_names.size())
		throw_ensemble_error("Ensemble.from_eigen_mat() cols != var_names.size");
	reals = _reals;
	var_names = _var_names;
	real_names = _real_names;
}

void Ensemble::set_eigen(Eigen::MatrixXd _reals)
{
	if (_reals.rows() != real_names.size())
		throw_ensemble_error("Ensemble.set_reals() rows != real_names.size");
	if (_reals.cols() != var_names.size())
		throw_ensemble_error("Ensemble.set_reals() cols != var_names.size");
	reals = _reals;
}

//Mat Ensemble::to_matrix()
//{
//	return to_matrix(real_names, var_names);
//}
//
//Mat Ensemble::to_matrix(vector<string> &row_names, vector<string> &col_names)
//{
//	return Mat();
//}

void Ensemble::reorder(vector<string> &_real_names, vector<string> &_var_names)
{
	reals = get_eigen(_real_names, _var_names);
	if (_var_names.size() != 0)
		var_names = _var_names;
	if (_real_names.size() != 0)
		real_names = _real_names;
}


void Ensemble::drop_rows(vector<int> &row_idxs)
{
	vector<int>::iterator start = row_idxs.begin(), end = row_idxs.end();
	vector<string> keep_names;
	for (int ireal = 0; ireal < reals.rows(); ireal++)
		if (find(start, end, ireal) == end)
			keep_names.push_back(real_names[ireal]);
	reals = get_eigen(keep_names, vector<string>());
	real_names = keep_names;
}

Eigen::MatrixXd Ensemble::get_eigen(vector<string> row_names, vector<string> col_names)
{
	vector<string> missing_rows,missing_cols;
	vector<string>::iterator iter, start = real_names.begin(), end = real_names.end();
	vector<int> row_idxs, col_idxs;
	Eigen::MatrixXd mat;
	//find row indices
	for (auto &rname : row_names)
	{
		iter = find(start, end, rname);
		if (iter == end)
			missing_rows.push_back(rname);
		row_idxs.push_back(iter - start);
	}
	//find col indices
	start = var_names.begin();
	end = var_names.end();
	for (auto cname : col_names)
	{
		iter = find(start, end, cname);
		if (iter == end)
			missing_cols.push_back(cname);
		col_idxs.push_back(iter - start);
	}

	// only mess with columns, keep rows the same
	if (row_names.size() == 0) 
	{		
		if (missing_cols.size() > 0)
			throw_ensemble_error("Ensemble.get_eigen() the following col_names not found:", missing_cols);
		mat.resize(real_names.size(), col_names.size());
		int j = 0;
		for (auto &jj : col_idxs)
		{
			mat.col(j) = reals.col(jj);
			j++;
		}
		return mat;
	}

	// only mess with rows, keep cols the same
	if (col_names.size() == 0) 
	{
		if (missing_rows.size() > 0)
			throw_ensemble_error("Ensemble.get_eigen() the following row_names not found:", missing_rows);
		mat.resize(row_names.size(), var_names.size());
		int i = 0;
		for (auto &ii : row_idxs)
		{
			mat.row(i) = reals.row(ii);
			i++;
		}
		return mat;
	}

	//the slow one where we are rearranging rows and cols
	if (missing_rows.size() > 0)
		throw_ensemble_error("Ensemble.get_eigen() the following row_names not found:", missing_rows);
	
	if (missing_cols.size() > 0)
		throw_ensemble_error("Ensemble.get_eigen() the following col_names not found:", missing_cols);
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
	csv << "real_name" << ',';
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
		csv << endl;
	}

}

const vector<string> Ensemble::get_real_names(vector<int> &indices)
{
	vector<string> names;
	for (auto &i : indices)
	{
		names.push_back(real_names[i]);
	}
	return names;
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
	//performance_log->log_event(full_message);
	//performance_log->~PerformanceLog();
	//file_manager.rec_ofstream() << endl << endl << "***" << full_message << endl;
	//file_manager.~FileManager();
	throw runtime_error(full_message);
}

Ensemble::~Ensemble()
{
}

vector<string> Ensemble::prepare_csv(const vector<string> &names, ifstream &csv, bool forgive)
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

void Ensemble::read_csv(int num_reals,ifstream &csv)
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




//ParameterEnsemble::ParameterEnsemble(const ParamTransformSeq &_par_transform, Pest &_pest_scenario,
//	FileManager &_file_manager, OutputFileWriter &_output_file_writer,
//	PerformanceLog *_performance_log, unsigned int seed):
//	Ensemble(_pest_scenario,_file_manager,_output_file_writer,_performance_log,seed),
//	par_transform(_par_transform)
ParameterEnsemble::ParameterEnsemble(Pest *_pest_scenario_ptr):Ensemble(_pest_scenario_ptr)
{
	par_transform = pest_scenario_ptr->get_base_par_tran_seq();
}

map<int,int> ParameterEnsemble::add_runs(RunManagerAbstract *run_mgr_ptr, vector<int> &real_idxs)
{
	map<int,int> real_run_ids;
	Parameters pars = get_pest_scenario_ptr()->get_ctl_parameters();
	//for (int ireal = 0; ireal < pe.shape().first; ireal++)
	Eigen::VectorXd evec;
	vector<double> svec;
	//for (auto &real_name : pe.get_real_names())
	int run_id;
	vector<string> run_real_names;
	if (real_idxs.size() > 0)
		for (auto i : real_idxs)
			run_real_names.push_back(real_names[i]);
	else
		run_real_names = real_names;

	for (auto &rname : run_real_names)
	{
		//evec = get_real_vector(rname);
		//const vector<double> svec(evec.data(), evec.data() + evec.size());
		pars.update_without_clear(var_names, get_real_vector(rname));
		if (tstat == ParameterEnsemble::transStatus::CTL)
			par_transform.ctl2model_ip(pars);
		else if (tstat == ParameterEnsemble::transStatus::NUM)
			par_transform.numeric2model_ip(pars);
		run_id = run_mgr_ptr->add_run(pars);
		real_run_ids[find(real_names.begin(), real_names.end(), rname) - real_names.begin()]  = run_id;
	}
	return real_run_ids;
}

void ParameterEnsemble::from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names)
{
	vector<string> missing;
	vector<string>::const_iterator start = _var_names.begin();
	vector<string>::const_iterator end = _var_names.end();

	for (auto &name : pest_scenario_ptr->get_ctl_ordered_par_names())
		if (find(start, end, name) == end)
			missing.push_back(name);
	if (missing.size() > 0)
		throw_ensemble_error("ParameterEnsemble.from_eigen_mat() the following par names no found: ", missing);
	Ensemble::from_eigen_mat(mat, _real_names, _var_names);
}


void ParameterEnsemble::from_csv(string &file_name)
{
	from_csv(file_name, pest_scenario_ptr->get_ctl_ordered_par_names());
}


void ParameterEnsemble::from_csv(string &file_name, const vector<string> &ordered_names)
{

	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening parameter csv " + file_name + " for reading"); 
	var_names = prepare_csv(ordered_names, csv, false);
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
	Ensemble::read_csv(num_reals, csv);

	//sort var_names according to ctl ordered par names
	vector<string> ordered_var_names;
	vector<string>::iterator start = var_names.begin();
	vector<string>::iterator end = var_names.end();

	for (auto &name : ordered_names)
		if (find(start, end, name) != end)
			ordered_var_names.push_back(name);
	if (var_names.size() != ordered_var_names.size())
	{
		stringstream ss;
		ss << "ordered names size " << ordered_var_names.size() << " != var names size " << var_names.size();
		throw_ensemble_error(ss.str());
	}
	if (var_names != ordered_var_names)
		reorder(vector<string>(), ordered_var_names);
	tstat = transStatus::CTL;
	//cout << reals << endl;
}

void ParameterEnsemble::enforce_bounds()
{
	throw_ensemble_error("pe.enforce_bounds() not implemented");
}

void ParameterEnsemble::to_csv(string &file_name)
{
	
	ofstream csv(file_name);
	if (!csv.good())
	{
		throw_ensemble_error("ParameterEnsemble.to_csv() error opening csv file " + file_name + " for writing");
	}
	csv << "real_name" << ',';
	for (auto &vname : var_names)
		csv << vname << ',';
	csv << endl;
	Parameters pars = pest_scenario_ptr->get_ctl_parameters();
	vector<double> svec;
	Eigen::VectorXd evec;
	svec.resize(reals.cols());
	for (int ireal = 0; ireal < reals.rows(); ireal++)
	{
		csv << real_names[ireal] << ',';
		//evec = reals.row(ireal);
		//svec.assign(evec.data(), evec.data() + evec.size());
		pars.update_without_clear(var_names, reals.row(ireal));
		if (tstat == transStatus::MODEL)
			par_transform.model2ctl_ip(pars);
		else if (tstat == transStatus::NUM)
			par_transform.numeric2ctl_ip(pars);

		for (auto &name : var_names)
			csv << pars[name] << ',';
		csv << endl;
	}
}

void ParameterEnsemble::transform_ip(transStatus to_tstat)
{
	if (to_tstat == tstat)
		return;
	if ((to_tstat == transStatus::NUM) && (tstat == transStatus::CTL))
	{
		Parameters pars = pest_scenario_ptr->get_ctl_parameters();
		for (int ireal = 0; ireal < reals.rows(); ireal++)
		{
			pars.update_without_clear(var_names, reals.row(ireal));
			par_transform.ctl2numeric_ip(pars);
			reals.row(ireal) = pars.get_data_eigen_vec(var_names);
		}

	}
	else
		throw_ensemble_error("ParameterEnsemble::transform_ip() only CTL to NUM implemented");
}


//ParameterEnsemble ParameterEnsemble::get_mean_diff()
//{
//	ParameterEnsemble new_pe = *this;
//	new_pe.set_reals(get_eigen_mean_diff());
//	return new_pe;
//}

//ObservationEnsemble::ObservationEnsemble(ObjectiveFunc *_obj_func, Pest &_pest_scenario,
//	FileManager &_file_manager, OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, unsigned int seed) :
//	Ensemble(_pest_scenario, _file_manager, _output_file_writer, _performance_log, seed), obj_func_ptr(_obj_func)
ObservationEnsemble::ObservationEnsemble(Pest *_pest_scenario_ptr): Ensemble(_pest_scenario_ptr)
{
}

void ObservationEnsemble::update_from_obs(int row_idx, Observations &obs)
{
	if (row_idx >= real_names.size())
		throw_ensemble_error("ObservtionEnsemble.update_from_obs() obs_idx out of range");
	reals.row(row_idx) = obs.get_data_eigen_vec(var_names);
}

vector<int> ObservationEnsemble::update_from_runs(map<int,int> &real_run_ids, RunManagerAbstract *run_mgr_ptr)
{
	set<int> failed_runs = run_mgr_ptr->get_failed_run_ids();
	vector<int> failed_real_idxs;
	Parameters pars = pest_scenario_ptr->get_ctl_parameters();
	Observations obs = pest_scenario_ptr->get_ctl_observations();
	//cout << oe.get_reals() << endl;
	for (auto &real_run_id : real_run_ids)
	{
		if (failed_runs.find(real_run_id.second) != failed_runs.end())
			failed_real_idxs.push_back(real_run_id.first);
		else
		{
			run_mgr_ptr->get_run(real_run_id.second, pars, obs);
			update_from_obs(real_run_id.first, obs);
		}
	}
	return failed_real_idxs;
}

void ObservationEnsemble::from_csv(string &file_name)
{
	from_csv(file_name, pest_scenario_ptr->get_ctl_ordered_obs_names());
}
void ObservationEnsemble::from_csv(string &file_name, const vector<string> &ordered_names)
{
	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening observation csv " + file_name + " for reading");
	var_names = prepare_csv(ordered_names, csv, false);
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
	Ensemble::read_csv(num_reals, csv);

	//sort var_names according to ctl ordered par names
	vector<string> ordered_var_names;
	vector<string>::iterator start = var_names.begin();
	vector<string>::iterator end = var_names.end();

	for (auto &name : ordered_names)
		if (find(start, end, name) != end)
			ordered_var_names.push_back(name);
	if (var_names.size() != ordered_var_names.size())
	{
		stringstream ss;
		ss << "ordered names size " << ordered_var_names.size() << " != var names size " << var_names.size();
		throw_ensemble_error(ss.str());
	}
	if (var_names != ordered_var_names)
		reorder(vector<string>(), ordered_var_names);
}

void ObservationEnsemble::from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names)
{
	vector<string> missing;
	vector<string>::const_iterator start = _var_names.begin();
	vector<string>::const_iterator end = _var_names.end();

	for (auto &name : pest_scenario_ptr->get_ctl_ordered_obs_names())
		if (find(start, end, name) == end)
			missing.push_back(name);
	if (missing.size() > 0)
		throw_ensemble_error("ObservationEnsemble.from_eigen_mat() the following obs names no found: ", missing);
	Ensemble::from_eigen_mat(mat, _real_names, _var_names);
}

//ObservationEnsemble ObservationEnsemble::get_mean_diff()
//{
//	ObservationEnsemble new_oe = *this;
//	new_oe.set_reals(get_eigen_mean_diff());
//	return new_oe;
//}



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

void EnsemblePair::queue_runs(RunManagerAbstract *run_mgr_ptr)
{

	Parameters pars = pe.get_pest_scenario_ptr()->get_ctl_parameters();
	//for (int ireal = 0; ireal < pe.shape().first; ireal++)
	Eigen::VectorXd evec;
	vector<double> svec;
	vector<string> var_names = pe.get_var_names();
	//for (auto &real_name : pe.get_real_names())
	real_run_ids.clear();
	int run_id;
	for (auto &act_idx : active_real_indices)
	{
		evec = pe.get_real_vector(act_idx);
		const vector<double> svec(evec.data(), evec.data() + evec.size());
		pars.update_without_clear(var_names, svec);
		if (pe.get_trans_status() == ParameterEnsemble::transStatus::CTL)
			pe.get_par_transform().ctl2model_ip(pars);
		else if (pe.get_trans_status() == ParameterEnsemble::transStatus::NUM)
			pe.get_par_transform().numeric2model_ip(pars);

		run_id = run_mgr_ptr->add_run(pars);
		real_run_ids[act_idx] = run_id;
	}
	//pe.set_trans_status(ParameterEnsemble::transStatus::MODEL);
	//cout << pe.get_reals() << endl;
}

void EnsemblePair::run(RunManagerAbstract *run_mgr_ptr)
{
	run_mgr_ptr->run();
}

vector<int> EnsemblePair::process_runs(RunManagerAbstract *run_mgr_ptr)
{
	
	set<int> failed_runs = run_mgr_ptr->get_failed_run_ids();
	vector<int> failed_real_idxs;
	Parameters pars = pe.get_pest_scenario_ptr()->get_ctl_parameters();
	Observations obs;
	//cout << oe.get_reals() << endl;
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
	//remove failed reals from active_idx
	vector<int>::iterator iter;
	for (auto &fi : failed_real_idxs)
	{
		iter = find(active_real_indices.begin(), active_real_indices.end(), fi);
		if (iter == active_real_indices.end())
		{
			stringstream ss;
			ss << "EnsemblePair.run() failed real idx " << fi << " not found in active_real_idxs";
			pe.throw_ensemble_error(ss.str());
		}
		active_real_indices.erase(iter);
	}
	//cout << oe.get_reals() << endl;
	//cout << pe.get_reals() << endl;
	return failed_real_idxs;
}

vector<string> EnsemblePair::get_oe_active_names()
{
	vector<string> act_real_names, real_names = oe.get_real_names();
	for (auto &i : active_real_indices)
		act_real_names.push_back(real_names[i]);
	return act_real_names;
}

vector<string> EnsemblePair::get_pe_active_names()
{
	vector<string> act_real_names, real_names = pe.get_real_names();
	for (auto &i : active_real_indices)
		act_real_names.push_back(real_names[i]);
	return act_real_names;
}


//these are too dangerous! dont use
//void EnsemblePair::set_pe(ParameterEnsemble &_pe)
//{
//	if (pe.shape() != _pe.shape())
//		pe.throw_ensemble_error("EnsemblePair::set_pe() pe.shape != _pe.shape");
//	pe = _pe;
//}
//
//void EnsemblePair::set_oe(ObservationEnsemble &_oe)
//{
//	if (oe.shape() != _oe.shape())
//		oe.throw_ensemble_error("EnsemblePair::set_oe() oe.shape != _oe.shape");
//	oe = _oe;
//}

//Eigen::MatrixXd EnsemblePair::get_active_oe_eigen()
//{
//	vector<string> act_real_names, oe_real_names = oe.get_real_names();
//	for (auto &i : active_real_indices)
//		act_real_names.push_back(oe_real_names[i]);
//	return oe.get_eigen(act_real_names, vector<string>());
//}
//
//Eigen::MatrixXd EnsemblePair::get_active_pe_eigen()
//{
//	vector<string> act_real_names, pe_real_names = pe.get_real_names();
//	for (auto &i : active_real_indices)
//		act_real_names.push_back(pe_real_names[i]);
//	return pe.get_eigen(act_real_names, vector<string>());
//}
//
//Eigen::MatrixXd EnsemblePair::get_active_pe_mean_diff()
//{
//	vector<string> act_real_names, pe_real_names = pe.get_real_names();
//	for (auto &i : active_real_indices)
//		act_real_names.push_back(pe_real_names[i]);
//	return pe.get_eigen_mean_diff(act_real_names);
//}
//
//Eigen::MatrixXd EnsemblePair::get_active_oe_mean_diff()
//{
//	vector<string> act_real_names, oe_real_names = oe.get_real_names();
//	for (auto &i : active_real_indices)
//		act_real_names.push_back(oe_real_names[i]);
//	return oe.get_eigen_mean_diff(act_real_names);
//
//}