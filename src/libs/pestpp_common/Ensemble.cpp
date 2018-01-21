#include <random>
#include <iomanip>
#include <unordered_set>
#include <iterator>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "ParamTransformSeq.h"
#include "ObjectiveFunc.h"
#include "RedSVD-h.h"
#include "covariance.h"

mt19937_64 Ensemble::rand_engine = mt19937_64(1);

Ensemble::Ensemble(Pest *_pest_scenario_ptr): pest_scenario_ptr(_pest_scenario_ptr)
{
	// initialize random number generator
	rand_engine.seed(1123433458);
}

void Ensemble::draw(int num_reals, Covariance &cov, Transformable &tran, const vector<string> &draw_names)
{
	if ((draw_names.size() > 50000) && (!cov.isdiagonal()))
		cout << "  ---  Ensemble::draw() warning: non-diagonal cov used to draw for lots of variables...this might run out of memory..." << endl << endl;
	
	Eigen::MatrixXd draws(num_reals, draw_names.size());
	draws.setZero();
	RedSVD::sample_gaussian(draws);
	if (cov.isdiagonal())
	{
		Eigen::VectorXd std = cov.e_ptr()->diagonal().cwiseSqrt();

		for (int j = 0; j < draw_names.size(); j++)
		{
			//cout << var_names[j] << " , " << std(j) << endl;
			draws.col(j) *= std(j);
		}
	}
	else
	{
		//totally arbitrary hack 
		int ncomps = std::min<int>(100,draw_names.size());
		RedSVD::RedSymEigen<Eigen::MatrixXd> eig;
		eig.compute(cov.e_ptr()->toDense(),ncomps);
		Eigen::MatrixXd evec(draw_names.size(), draw_names.size());
		evec.setZero();
		evec.leftCols(ncomps) = eig.eigenvectors();
		Eigen::VectorXd eval(draw_names.size());
		eval.setZero();
		eval.topRows(ncomps) = eig.eigenvalues().cwiseSqrt();

		//Eigen::MatrixXd proj = eig.eigenvectors() * eig.eigenvalues().cwiseSqrt().asDiagonal();
		Eigen::MatrixXd proj = evec * eval.asDiagonal();
		for (int i = 0; i < num_reals; i++)
		{
			/*Eigen::VectorXd v = reals.row(i);
			cout << v.rows() << " , " << v.cols() << endl;
			cout << proj.rows() << " , " << proj.cols() << endl;
			Eigen::MatrixXd t = proj * v;
			cout << t.rows() << " , " << t.cols() << endl;*/
			draws.row(i) = proj * draws.row(i).transpose();
		}
	}
	real_names.clear();
	stringstream ss;
	for (int i = 0; i < num_reals; i++)
	{
		ss.str("");
		ss << i;
		real_names.push_back(ss.str());
	}

	reals.resize(num_reals, var_names.size());
	reals.setZero();
	vector<string>::const_iterator start = draw_names.begin(), end=draw_names.end(), name;
	for (int j = 0; j < var_names.size(); j++)
	{
		int jj;
		name = find(start, end, var_names[j]);
		if (name != end)
		{
			jj = name - start;
			//cout << var_names[j] << " , " << tran.get_rec(var_names[j]) << endl;
			reals.col(j) = draws.col(jj).array() + tran.get_rec(var_names[j]);
		}
	}
}


Covariance Ensemble::get_diagonal_cov_matrix()
{
	
	Eigen::MatrixXd meandiff = get_eigen_mean_diff();
	double var;
	vector<Eigen::Triplet<double>> triplets;
	double num_reals = double(reals.rows());
	for (int j = 0; j < var_names.size(); j++)
	{
		var = (reals.col(j).cwiseProduct(reals.col(j))).sum() / num_reals;
		triplets.push_back(Eigen::Triplet<double>(j, j, var));
	}
	Eigen::SparseMatrix<double> mat;
	mat.conservativeResize(triplets.size(), triplets.size());
	mat.setFromTriplets(triplets.begin(), triplets.end());
	Covariance cov(var_names, mat, Covariance::MatType::DIAGONAL);
	return cov;
}

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
	if (keep_names.size() == 0)
		reals = Eigen::MatrixXd();
	else	
		reals = get_eigen(keep_names, vector<string>());
	real_names = keep_names;
}


void Ensemble::keep_rows(vector<int> &row_idxs)
{
	vector<int>::iterator start = row_idxs.begin(), end = row_idxs.end();
	vector<string> keep_names;
	for (int ireal = 0; ireal < reals.rows(); ireal++)
		if (find(start, end, ireal) != end)
			keep_names.push_back(real_names[ireal]);
	reals = get_eigen(keep_names, vector<string>());
	real_names = keep_names;
}


Eigen::MatrixXd Ensemble::get_eigen(vector<string> row_names, vector<string> col_names)
{
	vector<string> missing_rows,missing_cols;
	vector<string>::iterator iter, start = real_names.begin(), end = real_names.end();
	vector<int> row_idxs, col_idxs;

	//check for missing
	vector<string> missing;

	if (row_names.size() > 0)
	{
		map<string, int> real_map;
		for (int i = 0; i < real_names.size(); i++)
			real_map[real_names[i]] = i;
		set<string> real_set(real_names.begin(), end = real_names.end());
		set<string>::iterator end = real_set.end();
		for (auto &name : row_names)
		{
			if (real_set.find(name) == end)
				missing.push_back(name);
			row_idxs.push_back(real_map[name]);
		}
		if (missing.size() > 0)
			throw_ensemble_error("Ensemble.get_eigen() error: the following realization names were not found:", missing);
	}
	if (col_names.size() > 0)
	{
		map<string, int> var_map;
		for (int i = 0; i < var_names.size(); i++)
			var_map[var_names[i]] = i;
		set<string> var_set(var_names.begin(), var_names.end());
		set<string>::iterator end = var_set.end();
		for (auto &name : col_names)
		{
			if (var_set.find(name) == end)
				missing.push_back(name);
			col_idxs.push_back(var_map[name]);
		}
		if (missing.size() > 0)
			throw_ensemble_error("Ensemble.get_eigen() error: the following variable names were not found:", missing);
	}


	Eigen::MatrixXd mat;
	//find row indices
	//for (auto &rname : row_names)
	//{
	//	iter = find(start, end, rname);
	//	if (iter == end)
	//		missing_rows.push_back(rname);
	//	row_idxs.push_back(iter - start);
	//	row_idxs.push_back(row_map[rname]);
	//}
	//find col indices
	//start = var_names.begin();
	//end = var_names.end();
	/*for (auto cname : col_names)
	{
		iter = find(start, end, cname);
		if (iter == end)
			missing_cols.push_back(cname);
		col_idxs.push_back(iter - start);
		col_idxs.push_back(col_map[cname]);
	}*/

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

void Ensemble::set_real_names(vector<string>& _real_names)
{
	if (_real_names.size() != real_names.size())
	{
		stringstream ss;
		ss << " set_real_names() _real_names.size(): " << _real_names.size() << " != real_names.size(): " << real_names.size();
		throw_ensemble_error(ss.str());

	}
	real_names = _real_names;
}

Ensemble::~Ensemble()
{
}

map<string,int> Ensemble::prepare_csv(const vector<string> &names, ifstream &csv, bool forgive)
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
	unordered_set<string> hset(header_tokens.begin(), header_tokens.end());
	//cout << tokens << endl;
	//vector<string> header_tokens = tokens;

	// check for parameter names that in the pest control file but that are missing from the csv file
	vector<string> missing_names;
	unordered_set<string>::iterator end = hset.end();
	string name;
	for (auto &name : names)
		//if (find(header_tokens.begin(), header_tokens.end(), name) == header_tokens.end())
		if (hset.find(name) == end)
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
	map<string, int> header_info;
	hset.clear();
	hset = unordered_set<string>(names.begin(), names.end());
	end = hset.end();
	for (int i = 0; i < header_tokens.size(); i++)
	{
		//if (find(names.begin(), names.end(), header_tokens[i]) != names.end())
		if (hset.find(header_tokens[i]) != end)
		{
			//header_idxs.push_back(i);
			header_info[header_tokens[i]] = i;
			//var_names.push_back(header_tokens[i]);
			//header_names.push_back(header_tokens[i]);
		}
	}
	return header_info;

}

void Ensemble::add_to_cols(Eigen::MatrixXd &_reals, const vector<string> &_var_names)
{
	vector<string>::iterator start = var_names.begin(), end = var_names.end();
	vector<string> missing;
	for (auto &vname : _var_names)
		if (find(start, end, vname) == end)
			missing.push_back(vname);
	if (missing.size() > 0)
		throw_ensemble_error("add_to_cols(): the following var_names not found: ", missing);
	if (_var_names.size() != _reals.cols())
		throw_ensemble_error("add_to_cols(): _reals.cols() != _var_names.size()");
	if (_reals.rows() != reals.rows())
		throw_ensemble_error("add_to_cols(): _reals.rows() != reals.rows()");
	int i_this;
	string vname;
	for (int i = 0; i < _var_names.size(); i++)
	{
		vname = _var_names[i];
		i_this = find(start, end, vname) - start;
		reals.col(i_this) += _reals.col(i);
	}
}


void Ensemble::append_other_rows(Ensemble &other)
{
	if (other.shape().second != shape().second)
		throw_ensemble_error("append_other_rows(): different number of var_names in other");
	vector<string> probs;
	vector<string>::iterator start = var_names.begin(), end = var_names.end();
	for (auto &vname : other.get_var_names())
		if (find(start, end, vname) == end)
			probs.push_back(vname);
	if (probs.size() > 0)
		throw_ensemble_error("append_other_rows(): the following other::var_names not in this::var_names: ", probs);
	start = real_names.begin();
	end = real_names.end();
	for (auto &rname : other.get_real_names())

		if (find(start, end, rname) != end)
			probs.push_back(rname);
	if (probs.size() > 0)
		throw_ensemble_error("append_other_rows(): the following other::real_names are also in this::real_names: ", probs);
	vector<string> new_real_names = real_names;
	for (auto &rname : other.get_real_names())
		new_real_names.push_back(rname);
	//Eigen::MatrixXd org_reals = reals;
	//reals.resize(new_real_names.size(), var_names.size());
	reals.conservativeResize(new_real_names.size(), var_names.size());
		//for (int i = 0; i < org_reals.rows(), i++)
		//	reals.row(i) = 
	int iother = 0;
	for (int i = real_names.size(); i < new_real_names.size(); i++)
	{
		reals.row(i) = other.get_real_vector(iother);
		iother++;
	}
	real_names = new_real_names;
}


void Ensemble::from_binary(string &file_name, vector<string> &names, bool transposed)
{
	var_names.clear();
	real_names.clear();
	ifstream in;
	in.open(file_name.c_str(), ifstream::binary);

	int n_col;
	int n_nonzero;
	int n_row;
	int i, j, n;
	double data;
	char col_name[12];
	char row_name[20];

	// read header
	in.read((char*)&n_col, sizeof(n_col));
	in.read((char*)&n_row, sizeof(n_row));
	if (n_col > 0) throw runtime_error("Ensemble:::from_binary() error: binary matrix file " + file_name + " was produced by deprecated version of PEST");

	n_col = -n_col;
	n_row = -n_row;
	if ((n_col > 1e8) || (n_row > 1e8))
		throw_ensemble_error("Ensemble::from_binary() n_col or n_row > 1e8, something is prob wrong");
	////read number nonzero elements in jacobian (observations + prior information)
	in.read((char*)&n_nonzero, sizeof(n_nonzero));

	// record current position in file
	streampos begin_sen_pos = in.tellg();

	//advance to col names section
	in.seekg(n_nonzero*(sizeof(double) + sizeof(int)), ios_base::cur);

	//read col names	
	vector<string>* col_names = &var_names;
	vector<string>* row_names = &real_names;
	if (transposed)
	{
		col_names = &real_names;
		row_names = &var_names;
	}

	for (int i_rec = 0; i_rec<n_col; ++i_rec)
	{
		in.read(col_name, 12);
		string temp_col = string(col_name, 12);
		pest_utils::strip_ip(temp_col);
		pest_utils::upper_ip(temp_col);
		col_names->push_back(temp_col);
	}
	//read row names
	for (int i_rec = 0; i_rec<n_row; ++i_rec)
	{
		in.read(row_name, 20);
		string temp_row = pest_utils::strip_cp(string(row_name, 20));
		pest_utils::upper_ip(temp_row);
		row_names->push_back(temp_row);
	}

	//make sure that var_names is compatible with names
	if (var_names.size() != names.size())
	{
		set<string> vset(var_names.begin(), var_names.end());
		set<string> nset(names.begin(), names.end());
		vector<string> diff;
		set_symmetric_difference(vset.begin(), vset.end(), nset.begin(), nset.end(),std::back_inserter(diff));
		throw_ensemble_error("the following names are common between the var names in the binary file and the var names expected", diff);
	}

	//return to sensitivity section of file
	in.seekg(begin_sen_pos, ios_base::beg);

	// read matrix
	/*if (transposed)
		reals.resize(n_col, n_row);
	else*/
	reals.resize(n_row, n_col);
	reals.setZero();
	for (int i_rec = 0; i_rec<n_nonzero; ++i_rec)
	{
		in.read((char*)&(n), sizeof(n));
		--n;
		in.read((char*)&(data), sizeof(data));
		//j = int(n / (n_obs_and_pi)); // column index
		//i = (n - n_obs_and_pi*j) % n_obs_and_pi;  //row index
		j = int(n / (n_row)); // column index
		i = (n - n_row*j) % n_row;  //row index
	
		/*if (transposed)
			reals(j, i) = data;
		else*/
		reals(i, j) = data;
	}
	if (transposed)
		reals.transposeInPlace();
	in.close();
}



void Ensemble::read_csv(int num_reals,ifstream &csv, map<string,int> header_info)
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
		if (tokens[tokens.size() - 1].size() == 0)
			tokens.pop_back();
		//if (tokens.size() != var_names.size() + 1) // +1 for run id in first column
		//{
		//	stringstream ss;
		//	ss << "error parsing csv file on line " << lcount << ": wrong number of tokens, ";
		//	ss << "expecting " << var_names.size() + 1 << ", found " << tokens.size();
		//	throw runtime_error(ss.str());
		//}
		try
		{
			pest_utils::convert_ip(tokens[0], real_id);
		}
		catch (exception &e)
		{
			stringstream ss;
			ss << "error converting token '" << tokens[0] << "' to <int> run_id on line " << lcount << ": " << line << endl << e.what();
			throw runtime_error(ss.str());
		}
		real_names.push_back(real_id);
		//enum transStatus { CTL, NUM, MODEL };
		//for (int jcol = 1; jcol<tokens.size(); ++jcol)
		map<string, int> var_map;
		for (int i = 0; i < var_names.size(); i++)
			var_map[var_names[i]] = i;

		for (auto &hi : header_info)
		{
			try
			{
				val = pest_utils::convert_cp<double>(tokens[hi.second]);
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << "error converting token '" << tokens[hi.second] << "' to double for " << hi.first << " on line " << lcount << " : " << e.what();
				throw runtime_error(ss.str());
			}
			//vals.push_back(val);
			reals(irow, var_map[hi.first]) = val;
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




ParameterEnsemble::ParameterEnsemble(Pest *_pest_scenario_ptr):Ensemble(_pest_scenario_ptr)
{
	par_transform = pest_scenario_ptr->get_base_par_tran_seq();
}

void ParameterEnsemble::draw(int num_reals, Covariance &cov)
{
	var_names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
	Parameters par = pest_scenario_ptr->get_ctl_parameters();
	par_transform.active_ctl2numeric_ip(par);
	tstat = transStatus::NUM;
	Ensemble::draw(num_reals, cov, par, var_names);
	enforce_bounds();
	
}

void ParameterEnsemble::set_pest_scenario(Pest *_pest_scenario)
{
	pest_scenario_ptr = _pest_scenario;
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
			par_transform.active_ctl2model_ip(pars);
		else if (tstat == ParameterEnsemble::transStatus::NUM)
			par_transform.numeric2model_ip(pars);
		run_id = run_mgr_ptr->add_run(pars);
		real_run_ids[find(real_names.begin(), real_names.end(), rname) - real_names.begin()]  = run_id;
	}
	return real_run_ids;
}

void ParameterEnsemble::from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names, ParameterEnsemble::transStatus _tstat)
{
	vector<string> missing;
	vector<string>::const_iterator start = _var_names.begin();
	vector<string>::const_iterator end = _var_names.end();

	for (auto &name : pest_scenario_ptr->get_ctl_ordered_par_names())
		if (find(start, end, name) == end)
			missing.push_back(name);
	if (missing.size() > 0)
		throw_ensemble_error("ParameterEnsemble.from_eigen_mat() the following par names not found: ", missing);
	Ensemble::from_eigen_mat(mat, _real_names, _var_names);
	tstat = _tstat;
}


void ParameterEnsemble::from_binary(string &file_name)
{
	//Ensemble::from_binary(file_name, pest_scenario_ptr->get_ctl_ordered_par_names(), false);
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_par_names();
	Ensemble::from_binary(file_name, names, false);
	tstat = transStatus::CTL;
}

void ParameterEnsemble::from_csv(string &file_name)
{
	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening parameter csv " + file_name + " for reading"); 
	var_names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
	map<string,int>header_info = prepare_csv(var_names, csv, false);
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
	Ensemble::read_csv(num_reals, csv, header_info);

	//sort var_names according to ctl ordered par names
	/*vector<string> ordered_var_names;
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
		reorder(vector<string>(), ordered_var_names);*/
	tstat = transStatus::CTL;
	//cout << reals << endl;
}

void ParameterEnsemble::enforce_bounds()
{
	//throw_ensemble_error("pe.enforce_bounds() not implemented");
	if (tstat != ParameterEnsemble::transStatus::NUM)
	{
		throw_ensemble_error("pe.enforce_bounds() tstat != NUM not implemented");

	}
	ParameterInfo pinfo = pest_scenario_ptr->get_ctl_parameter_info();
	Parameters lower = pest_scenario_ptr->get_ctl_parameter_info().get_low_bnd(var_names);
	Parameters upper = pest_scenario_ptr->get_ctl_parameter_info().get_up_bnd(var_names);
	par_transform.ctl2numeric_ip(lower);
	par_transform.ctl2numeric_ip(upper);
	Eigen::VectorXd col;
	double l, u, v;
	for (int j = 0; j < reals.cols(); j++)
	{
		l = lower[var_names[j]];
		u = upper[var_names[j]];
		col = reals.col(j);
		//cout << var_names[j] <<',' << l <<','<< u << endl;
		//cout << col << endl << endl;
		for (int i = 0; i < reals.rows(); i++)
		{
			v = col(i);
			v = v < l ? l : v;
			col(i) = v > u ? u : v;
		}
		//cout << col << endl;
		reals.col(j) = col;
	}

}

void ParameterEnsemble::to_csv(string &file_name)
{
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_par_names();
	ofstream csv(file_name);
	if (!csv.good())
	{
		throw_ensemble_error("ParameterEnsemble.to_csv() error opening csv file " + file_name + " for writing");
	}
	csv << "real_name" << ',';
	for (auto &vname : names)
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

		for (auto &name : names)
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
		vector<string> adj_par_names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
		Eigen::MatrixXd new_reals = Eigen::MatrixXd(shape().first, adj_par_names.size());
		for (int ireal = 0; ireal < reals.rows(); ireal++)
		{
			pars.update_without_clear(var_names, reals.row(ireal));
			par_transform.ctl2numeric_ip(pars);
			assert(pars.size() == adj_par_names.size());
			new_reals.row(ireal) = pars.get_data_eigen_vec(adj_par_names);
		}
		reals = new_reals;
		var_names = adj_par_names;
		tstat = to_tstat;

	}
	else
		throw_ensemble_error("ParameterEnsemble::transform_ip() only CTL to NUM implemented");
	
}
Covariance ParameterEnsemble::get_diagonal_cov_matrix()
{
	transform_ip(transStatus::NUM);
	return Ensemble::get_diagonal_cov_matrix();
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

void ObservationEnsemble::draw(int num_reals, Covariance &cov)
{
	//var_names = pest_scenario_ptr->get_ctl_ordered_nz_obs_names();
	var_names = pest_scenario_ptr->get_ctl_ordered_obs_names();
	Observations obs = pest_scenario_ptr->get_ctl_observations();
	Ensemble::draw(num_reals, cov, obs, pest_scenario_ptr->get_ctl_ordered_nz_obs_names());
}

void ObservationEnsemble::update_from_obs(int row_idx, Observations &obs)
{
	if (row_idx >= real_names.size())
		throw_ensemble_error("ObservtionEnsemble.update_from_obs() obs_idx out of range");
	//Eigen::VectorXd temp = obs.get_data_eigen_vec(var_names);
	//cout << reals.row(row_idx) << endl << endl;
	reals.row(row_idx) = obs.get_data_eigen_vec(var_names);
	//cout << reals.row(row_idx) << endl;
	return;
}

void ObservationEnsemble::update_from_obs(string real_name, Observations &obs)
{
	vector<string>::iterator real,start = real_names.begin(), end = real_names.end();
	real = find(start, end, real_name);
	if (real == end)
	{
		throw_ensemble_error("real_name not in real_names:" + real_name);
	}
	update_from_obs(real - start, obs);
}

vector<int> ObservationEnsemble::update_from_runs(map<int,int> &real_run_ids, RunManagerAbstract *run_mgr_ptr)
{
	set<int> failed_runs = run_mgr_ptr->get_failed_run_ids();
	vector<int> failed_real_idxs;
	Parameters pars = pest_scenario_ptr->get_ctl_parameters();
	Observations obs = pest_scenario_ptr->get_ctl_observations();
	//cout << oe.get_reals() << endl;
	string real_name;
	for (auto &real_run_id : real_run_ids)
	{
		if (failed_runs.find(real_run_id.second) != failed_runs.end())
			failed_real_idxs.push_back(real_run_id.first);
		else
		{
			run_mgr_ptr->get_run(real_run_id.second, pars, obs);
			//real_name = real_names[real_run_id.first];
			update_from_obs(real_run_id.first, obs);
		}
	}
	return failed_real_idxs;
}

void ObservationEnsemble::from_binary(string &file_name)
{
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_obs_names();
	Ensemble::from_binary(file_name, names, true);
}

void ObservationEnsemble::from_csv(string &file_name)
{
	var_names = pest_scenario_ptr->get_ctl_ordered_obs_names();
	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening observation csv " + file_name + " for reading");
	map<string,int> header_info = prepare_csv(pest_scenario_ptr->get_ctl_ordered_nz_obs_names(), csv, false);
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
	Ensemble::read_csv(num_reals, csv,header_info);

	//sort var_names according to ctl ordered var names
	/*vector<string> ordered_var_names;
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
	}*/
	//if (var_names != ordered_var_names)
	//	reorder(vector<string>(), ordered_var_names);
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
