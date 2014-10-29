
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>

#include "Pest.h"
#include "utilities.h"
#include "covariance.h"

using namespace std;

//---------------------------------------
//Mat constructors
//---------------------------------------
Mat::Mat(vector<string> _row_names, vector<string> _col_names, Eigen::SparseMatrix<double> _matrix)
{
	row_names = _row_names;
	col_names = _col_names;
	assert(row_names.size() == _matrix.rows());
	assert(col_names.size() == _matrix.cols());
	matrix = _matrix;
}

Mat::Mat(vector<string> _row_names, vector<string> _col_names, Eigen::SparseMatrix<double> _matrix, MatType _mattype)
{
	row_names = _row_names;
	col_names = _col_names;
	assert(row_names.size() == _matrix.rows());
	assert(col_names.size() == _matrix.cols());
	matrix = _matrix;
	mattype = _mattype;
}


//---------------------------------------
//Mat operator
//--------------------------------------
Mat Mat::operator*(Mat &other_mat)
{
	//find common this.col_names and other_mat.row_names
	//if (row_names != other_mat.get_col_names())
	MatType new_mattype = MatType::DENSE;
	if ((mattype == MatType::DIAGONAL) && (other_mat.get_mattype() == MatType::DIAGONAL))
		new_mattype = MatType::DIAGONAL;
	vector<string> common;
	for (auto &row_name : other_mat.get_row_names())
	{
		if (find(col_names.begin(), col_names.end(), row_name) != col_names.end())
			common.push_back(row_name);
	}
	if (common.size() == 0)
	{
		cout << "no common this.col_names/other_mat.row_names:" << endl << "this.col_names:" << endl;
		for (auto &name : col_names) cout << name << ',';
		cout << endl << "other_mat.row_names:" << endl;
		for (auto &name : other_mat.get_row_names()) cout << name << ',';
		cout << endl;
		throw runtime_error("no common elements found in Mat::operator*");
	}
	//no realignment needed...
	if ((common == col_names) && (common == other_mat.get_row_names()))
	{
		Eigen::SparseMatrix<double> new_matrix = matrix * other_mat.get_matrix();
		Mat new_mat(row_names, other_mat.col_names, new_matrix, new_mattype);
		return new_mat;
	}

	//only other_mat needs realignment
	else if (common == col_names)
	{
		Mat new_other_mat = other_mat.get(common, other_mat.get_col_names());
		Eigen::SparseMatrix<double> new_matrix = matrix * new_other_mat.get_matrix();
		Mat new_mat(common, other_mat.col_names, new_matrix, new_mattype);
		return new_mat;
	}
	
	//only this needs realignment
	else if (common == other_mat.get_row_names())
	{
		Mat new_this = get(row_names, common);
		Eigen::SparseMatrix<double> new_matrix = new_this.get_matrix() * other_mat.get_matrix();
		Mat new_mat(row_names, other_mat.col_names, new_matrix, new_mattype);
		return new_mat;

	}
	//both need realignment
	else
	{
		Mat new_other_mat = other_mat.get(common, other_mat.get_col_names());
		Mat new_this = get(row_names, common);
		Eigen::SparseMatrix<double> new_matrix = new_this.get_matrix() * new_other_mat.get_matrix();
		Mat new_mat(row_names, other_mat.col_names, new_matrix, new_mattype);
		return new_mat;
	}


}


void Mat::transpose()
{
	if (mattype != MatType::DIAGONAL)
	{
		matrix = matrix.transpose();
		vector<string> temp = row_names;
		row_names = col_names;
		col_names = temp;
	}
}

//-----------------------------------------
//Mat IO
//-----------------------------------------
void Mat::to_ascii(const string &filename)
{
	ofstream out(filename);
	if (!out.good())
	{
		throw runtime_error("cannot open " + filename + " to write ASCII matrix");
	}
	out << setw(6) << nrow() << setw(6) << ncol() << setw(6) << icode << endl;
	out << matrix;
	if (icode == 1)
	{
		out<< "* row and column names" << endl;
		for (auto &name : row_names)
		{
			out << name << endl;
		}
	}
	else
	{
		out << "* row names" << endl;
		for (auto &name : row_names)
		{
			out << name << endl;
		}
		out << "* column names" << endl;
		for (auto &name : col_names)
		{
			out << name << endl;
		}
	}
	out.close();
}

void Mat::from_ascii(const string &filename)
{
	ifstream in(filename);
	if (!in.good())
	{
		throw runtime_error("cannot open " + filename + " to read ASCII matrix");
	}
	int nrow = -999, ncol = -999;
	if (in >> nrow >> ncol >> icode)
	{

	}
	else
	{
		throw runtime_error("error reading nrow ncol icode from first line of ASCII matrix file: " + filename);
	}
	//vector<double> vals;	
	vector<Eigen::Triplet<double>> triplet_list;
	double val;
	int irow = 0, jcol = 0;
	for (int inode = 0; inode < nrow*ncol;inode++)
	{
		if (in >> val)
		{			
			if (val != 0.0)
			{
				triplet_list.push_back(Eigen::Triplet<double>(irow,jcol,val));
			}	
			jcol++;
			if (jcol >= ncol)
			{
				irow++;
				jcol = 0;
			}
		}
		else
		{
			string i_str = to_string(inode);
			throw runtime_error("error reading entry number "+i_str+" from ASCII matrix file: "+filename);
		}
	}
	//read the newline char
	string header;
	getline(in, header);
	if (!getline(in,header))
	{
		throw runtime_error("error reading row/col description line from ASCII matrix file: " + filename);
	}
	pest_utils::upper_ip(header);
	string name;
	if (icode == 1)
	{
		if (nrow != ncol)
			throw runtime_error("nrow != ncol for icode type 1 ASCII matrix file:" + filename);
		if((header.find("ROW") == string::npos) || (header.find("COLUMN") == string::npos))
			throw runtime_error("expecting row and column names header instead of:" + header + " in ASCII matrix file: " + filename);
		try
		{
			row_names = read_namelist(in, nrow);
		}
		catch (exception &e)
		{
			throw runtime_error("error reading row/column names from ASCII matrix file: " + filename + "\n" + e.what());
		}
		if ((nrow != row_names.size()) || (ncol != row_names.size()))
			throw runtime_error("nuymber of row/col names does not match matrix dimensions");
		col_names = row_names;
	}
	else
	{
		if(header.find("ROW") == string::npos)
			throw runtime_error("expecting row names header instead of:" + header + " in ASCII matrix file: " + filename);
		try
		{
			row_names = read_namelist(in, nrow);
		}
		catch (exception &e)
		{
			throw runtime_error("error reading row names from ASCII matrix file: " + filename + "\n" + e.what());
		}
		if (!getline(in, header))
		{
			throw runtime_error("error reading column name descriptor from ASCII matrix file: " + filename);
		}
		pest_utils::upper_ip(header);
		if (header.find("COLUMN") == string::npos)
			throw runtime_error("expecting column names header instead of:" + header + " in ASCII matrix file: " + filename);
		try
		{
			col_names = read_namelist(in, ncol);
		}
		catch (exception &e)
		{
			throw runtime_error("error reading column names from ASCII matrix file: " + filename + "\n" + e.what());
		}
		if (nrow != row_names.size())
			throw runtime_error("nrow != row_names.size() in ASCII matrix file: " + filename);

		if(ncol != col_names.size())
			throw runtime_error("ncol != col_names.size() in ASCII matrix file: " + filename);

	}
	in.close();


	//
	Eigen::SparseMatrix<double> new_matrix(nrow, ncol);
	new_matrix.setZero();  // initialize all entries to 0
	new_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	matrix = new_matrix;
}

vector<string> Mat::read_namelist(ifstream &in, int &nitems)
{
	vector<string> names;
	string name;
	for (int i = 0; i < nitems; i++)
	{
		if (!getline(in, name))
		{
			string i_str = to_string(i);
			throw runtime_error("error name for entry " + i_str);
		}
		if (name.find("*") != string::npos)
		{
			string i_str = to_string(i);
			throw runtime_error("'*' found in item name: " + name+", item number: "+i_str);
		}
		pest_utils::upper_ip(name);
		if (find(names.begin(), names.end(), name) != names.end())
			throw runtime_error("duplicate name: " + name + " found in name list");
		names.push_back(name);
	}
	return names;
}

void Mat::to_binary(const string &filename)
{

}

void Mat::from_binary(const string &filename)
{

}


//-----------------------------------------
//Maninpulate the shape and ordering of Mats
//-----------------------------------------
void Mat::align(vector<string> &other_row_names, vector<string> &other_col_names)
{

}

Mat Mat::get(vector<string> &new_row_names, vector<string> &new_col_names)
{
	//check that every row and col name is listed
	vector<string> row_not_found;
	for (auto &n_row_name : new_row_names)
	{
		if (find(row_names.begin(), row_names.end(), n_row_name) == row_names.end())
			row_not_found.push_back(n_row_name);
	}
	vector<string> col_not_found;
	for (auto &n_col_name : new_col_names)
	{
		if (find(col_names.begin(), col_names.end(), n_col_name) == col_names.end())
			col_not_found.push_back(n_col_name);
	}

	if (row_not_found.size() != 0)
	{
		cout << "error in Mat::get(): the following row names were not found:" << endl;
		for (auto &name : row_not_found)
			cout << name << ",";
		cout << endl;
	}

	if (col_not_found.size() != 0)
	{
		cout << "error in Mat::get(): the following col names were not found:" << endl;
		for (auto &name : col_not_found)
			cout << name << ",";
		cout << endl;
	}

	if ((row_not_found.size() != 0) || (col_not_found.size() != 0))
	{
		throw runtime_error("atleast one row or col name not found in Mat::get()");
	}


	int nrow = new_row_names.size();
	int ncol = new_col_names.size();
	int irow_new;
	int icol_new;

	unordered_map<string, int> row_name2newindex_map;
	unordered_map<string, int> col_name2new_index_map;

	// Build mapping of parameter names to column number in new matrix to be returned
	icol_new = 0;
	for (vector<string>::const_iterator b = new_col_names.begin(), e = new_col_names.end();
		b != e; ++b, ++icol_new) {
		col_name2new_index_map[(*b)] = icol_new;
	}

	// Build mapping of observation names to row  number in new matrix to be returned
	irow_new = 0;
	for (vector<string>::const_iterator b = new_row_names.begin(), e = new_row_names.end();
		b != e; ++b, ++irow_new) {
		row_name2newindex_map[(*b)] = irow_new;
	}
	
	unordered_map<string, int>::const_iterator found_col;
	unordered_map<string, int>::const_iterator found_row;
	unordered_map<string, int>::const_iterator not_found_col_map = col_name2new_index_map.end();
	unordered_map<string, int>::const_iterator not_found_row_map = row_name2newindex_map.end();

	const string *row_name;
	const string *col_name;
	std::vector<Eigen::Triplet<double> > triplet_list;
	for (int icol = 0; icol<matrix.outerSize(); ++icol)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, icol); it; ++it)
		{			
			col_name = &col_names[it.col()];
			row_name = &row_names[it.row()];
			found_col = col_name2new_index_map.find(*row_name);
			found_row = row_name2newindex_map.find(*col_name);

			if (found_col != not_found_col_map && found_row != not_found_row_map)
			{
				triplet_list.push_back(Eigen::Triplet<double>(found_row->second, found_col->second, it.value()));
			}
		}
	}
	Eigen::SparseMatrix<double> new_matrix(nrow, ncol);
	new_matrix.setZero();
	new_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return Mat(new_row_names,new_col_names,new_matrix,mattype);
}

Mat Mat::extract(vector<string> &ext_row_names, vector<string> &ext_col_names)
{
	Mat new_mat = get(ext_row_names, ext_col_names);
	drop(ext_row_names, ext_col_names);
	return new_mat;
}

void Mat::drop(vector<string> &drop_row_names, vector<string> &drop_col_names)
{
	vector<string> new_row_names, new_col_names;
	vector<Eigen::Triplet<double>> triplet_list;

	//check that each of the drop_row_names and drop_col_names in row_names and col_names




	 

}



//-----------------------------------------
//covariance matrices
//-----------------------------------------
Covariance::Covariance(vector<string> &names)
{
	row_names = names;
	col_names = names;
	icode = 1;
}

Covariance::Covariance()
{
	icode = 1;
}

Covariance::Covariance(vector<string> _names, Eigen::SparseMatrix<double> _matrix)
{	
	if ((_names.size() != _matrix.rows()) || (_names.size() != _matrix.cols()))
		throw runtime_error("Covariance::Covariance(): names.size() does not match matrix dimensions");
	matrix = _matrix;
	row_names = _names;
	col_names = _names;
	icode = 1;
}

Covariance::Covariance(Mat _mat)
{
	if (_mat.get_row_names() != _mat.get_col_names())
		throw runtime_error("error instantiating Covariance from Mat: row_names != col_names");
	row_names = _mat.get_row_names();
	col_names = _mat.get_col_names();
	matrix = _mat.get_matrix();
	icode = 1;
	mattype = _mat.get_mattype();
}

Covariance Covariance::get(vector<string> &other_names)
{
	Covariance new_cov(Mat::get(other_names, other_names));
	return new_cov;

}

void Covariance::from_uncertainty_file(const string &filename)
{
	ifstream in(filename);
	if (!in.good())
	{
		throw runtime_error("cannot open " + filename + " to read uncertainty file: "+filename);
	}
	mattype = MatType::DIAGONAL;
	
	vector<Eigen::Triplet<double>> triplet_list;
	vector<string> names;
	string line,name;
	double val;
	vector<string> tokens;
	int irow=0, jcol=0;
	
	while (getline(in, line))
	{
		pest_utils::upper_ip(line);
		//if this is the start of some block
		if (line.find("START") != string::npos)
		{
			if (line.find("STANDARD_DEVIATION") != string::npos)
			{
				while (true)
				{
					if (!getline(in, line))
						throw runtime_error("EOF encountered while reading standard_deviation block\
							from uncertainty file:" + filename);
					pest_utils::upper_ip(line);
					if (line.find("END") != string::npos) break;

					tokens.clear();
					pest_utils::tokenize(line, tokens);
					pest_utils::convert_ip(tokens[1], val);					
					if (find(names.begin(), names.end(), name) != names.end())
						throw runtime_error(name + " listed more than once in uncertainty file:" + filename);
					names.push_back(tokens[0]);
					triplet_list.push_back(Eigen::Triplet<double>(irow, jcol, val));
					irow++, jcol++;
				}

			}
			else if (line.find("COVARIANCE_MATRIX") != string::npos)
			{
				string cov_filename = "none";
				double var_mult = 1.0;
				while (true)
				{
					if (!getline(in, line))
						throw runtime_error("EOF encountered while reading covariance_matrix block\
							from uncertainty file:" + filename);
					pest_utils::upper_ip(line);
					if (line.find("END") != string::npos) break;

					tokens.clear();
					pest_utils::tokenize(line, tokens);					
					if (tokens[0].find("FILE") != string::npos)
						cov_filename = tokens[1];
					else if (tokens[0].find("VARIANCE") != string::npos)
						pest_utils::convert_ip(tokens[1], var_mult);
					else
						throw runtime_error("unrecognized token:" + tokens[0] + " in covariance matrix block in uncertainty file:" + filename);
				}
				//read the covariance matrix
				Covariance cov;
				cov.from_ascii(cov_filename);

				//check that the names in the covariance matrix are not already listed
				vector<string> dup_names;
				for (auto &name : cov.get_row_names())
				{
					if (find(names.begin(), names.end(), name) != names.end())
						dup_names.push_back(name);
					else
						names.push_back(name);
				}
				if (dup_names.size() != 0)
				{
					cout << "the following names from covariance matrix file " << cov_filename << " have already be found in uncertainty file " << filename << endl;
					for (auto &name : dup_names)
						cout << name << ',';
					cout << endl;
					throw runtime_error("atleast one name in covariance matrix " + cov_filename + " is already listed in uncertainty file: " + filename);
				}

				//build triplets from the covariance matrix
				int start_irow = irow;
				Eigen::SparseMatrix<double> cov_matrix = cov.get_matrix();
				for (int icol = 0; icol < cov_matrix.outerSize(); ++icol)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(cov_matrix, icol); it; ++it)
					{
						triplet_list.push_back(Eigen::Triplet<double>(irow, jcol, it.value()));
						irow++;
					}
					jcol++;
					irow = start_irow;
				}
				mattype = MatType::BLOCK;
				irow = jcol;
			}
			else
				throw runtime_error("unrecognized block:" + line + " in uncertainty file:" + filename);
		}
	}

	Eigen::SparseMatrix<double> new_matrix(names.size(), names.size());
	new_matrix.setZero();  // initialize all entries to 0
	new_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	matrix = new_matrix;
	row_names = names;
	col_names = names;
}

void Covariance::from_parameter_bounds(Pest &pest_scenario)
{
	vector<Eigen::Triplet<double>> triplet_list;
	const ParameterRec* par_rec;
	int i = 0;
	double upper, lower;
	for (auto par_name : pest_scenario.get_ctl_ordered_par_names())
	{
		pest_utils::upper_ip(par_name);
		par_rec = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(par_name);
		upper = par_rec->ubnd;
		lower = par_rec->lbnd;
		if (par_rec->tranform_type == ParameterRec::TRAN_TYPE::LOG)
		{
			upper = log10(upper);
			lower = log10(lower);
		}
		if ((par_rec->tranform_type != ParameterRec::TRAN_TYPE::FIXED) && (par_rec->tranform_type != ParameterRec::TRAN_TYPE::TIED))
		{
			row_names.push_back(par_name);
			col_names.push_back(par_name);
			triplet_list.push_back(Eigen::Triplet<double>(i, i, pow((upper - lower) / 4.0, 2.0)));
			i++;
		}
	}
	if (triplet_list.size() > 0)
	{
		matrix.resize(row_names.size(), row_names.size());
		matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	}
	else
	{
		throw runtime_error("Error loading covariance from parameter bounds: no non-fixed/non-tied parameters found");
	}
	mattype = Mat::MatType::DIAGONAL;
}

void Covariance::from_parameter_bounds(const string &pst_filename)
{
	ifstream ipst(pst_filename);
	Pest pest_scenario;
	pest_scenario.process_ctl_file(ipst, pst_filename);
	from_parameter_bounds(pest_scenario);
}

void Covariance::from_observation_weights(const string &pst_filename)
{
	ifstream ipst(pst_filename);
	Pest pest_scenario;
	pest_scenario.process_ctl_file(ipst, pst_filename);
	from_observation_weights(pest_scenario);

}

void Covariance::from_observation_weights(Pest &pest_scenario)
{
	vector<Eigen::Triplet<double>> triplet_list;
	const ObservationRec* obs_rec;
	int i = 0;	
	for (auto obs_name : pest_scenario.get_ctl_ordered_obs_names())
	{
		pest_utils::upper_ip(obs_name);
		obs_rec = pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(obs_name);
		if (obs_rec->weight > 0.0)
		{
			row_names.push_back(obs_name);
			col_names.push_back(obs_name);
			triplet_list.push_back(Eigen::Triplet<double>(i, i, pow(1.0 / obs_rec->weight, 2.0)));
			i++;
		}
	}
	if (row_names.size() > 0)
	{
		matrix.resize(row_names.size(), row_names.size());
		matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	}
	else
	{
		throw runtime_error("Error loading covariance from obs weights: no non-zero weighted obs found");
	}
	mattype = Mat::MatType::DIAGONAL;


}

void Covariance::to_uncertainty_file(const string &filename)
{
	ofstream out(filename);
	if (!out.good())
	{
		throw runtime_error("Error opening file: " + filename + " to write an uncertainty file");
	}

	//check if diagonal, write stdevs
	if (mattype == Mat::MatType::DIAGONAL)
	{
		Eigen::VectorXd vec(matrix.diagonal());
		out << "START STANDARD_DEVIATION" << endl;
		int i=0;
		for (vector<string>::iterator name = row_names.begin(); name != row_names.end(); ++name, i++)
		{
			out << "  " << setw(20) << left << *name << "  " << setw(20) << left << vec(i) << endl;
		}
		out << "END STANDARD_DEVIATION" << endl;
		out.close();
	}
	else
	{
		out << "START COVARIANCE_MATRIX" << endl;
		out << "  file emu_cov.mat" << endl;
		out << " variance multiplier 1.0" << endl;
		out << "END COVARIANCE_MATRIX" << endl;
		out.close();
		to_ascii("emu_cov.mat");
	}

}
