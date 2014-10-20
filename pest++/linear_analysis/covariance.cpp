
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>

#include "Pest.h"
#include "covariance.h"

using namespace std;


Mat::Mat(vector<string> _row_names, vector<string> _col_names)
{
	row_names = _row_names;
	col_names = _col_names;
	matrix = Eigen::SparseMatrix<double>(row_names.size(),col_names.size());
}

Mat::Mat(vector<string> _row_names, vector<string> _col_names, Eigen::SparseMatrix<double> _matrix)
{
	row_names = _row_names;
	col_names = _col_names;
	assert(row_names.size() == _matrix.rows());
	assert(col_names.size() == _matrix.cols());
	matrix = _matrix;
}

Mat Mat::get_matrix(const vector<string> &row_names, const vector<string> &col_names)
{
	Mat new_matrix;
	return new_matrix;
}

void Mat::to_ascii(const string &filename)
{
	ofstream out(filename);
	if (!out.good())
	{
		throw runtime_error("cannot open " + filename + " to write ASCII matrix");
	}
	out << setw(6) << nrow() << setw(6) << ncol() << setw(6) << icode << endl;
	//Eigen::MatrixXd dmatrix(matrix)
	//int i, j;	
	//for (int irow = 0; irow < dmatrix.rows(); irow++)
	//{		
	//	//for (Eigen::SparseMatrix<double>::InnerIterator jcol(matrix, irow); jcol; ++jcol)
	//	//for (int jcol = 0; jcol < dmatrix.cols(); jcol++)
	//	for (auto &row : dmatrix.row(irow))
	//	{						
	//		out << left << setw(20) << dmatrix.value(irow,jcol);
	//	}
	//	out << endl;
	//}
	/*stringstream ss;
	ss << matrix;
	string smatrix(ss.str());
	string delim = "$";
	out << smatrix.substr(smatrix.find(delim)+3,smatrix.size());*/

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

	string name;
	if (icode == 1)
	{
		assert(nrow == ncol);
		assert((header.find("row") != string::npos) && (header.find("column") != string::npos));
		try
		{
			row_names = read_namelist(in, nrow);
		}
		catch (exception &e)
		{
			throw runtime_error("error reading row/column names from ASCII matrix file: " + filename + "\n" + e.what());
		}
		assert(nrow == row_names.size());
		assert(ncol == row_names.size());
	}
	else
	{
		assert(header.find("row") != string::npos);
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
		assert(header.find("column") != string::npos);
		try
		{
			col_names = read_namelist(in, ncol);
		}
		catch (exception &e)
		{
			throw runtime_error("error reading column names from ASCII matrix file: " + filename + "\n" + e.what());
		}
		assert(nrow == row_names.size());
		assert(ncol == col_names.size());
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

void Mat::align_rows(vector<string> &other_row_names)
{

}

void Mat::align_cols(vector<string> &other_col_names)
{

}

void Mat::align(vector<string> &other_row_names, vector<string> &other_col_names)
{

}

Covariance::Covariance(vector<string> names)
{

}

Covariance::Covariance()
{

}

void Covariance::from_uncertainty_file(string &filename)
{

}

void Covariance::from_parameter_bounds(Pest &pest_scenario)
{
	
}

void Covariance::from_parameter_bounds(string &pst_filename)
{
	ifstream ipst(pst_filename);
	pest_scenario.process_ctl_file(ipst, pst_filename);
	from_parameter_bounds(pest_scenario);
}

void Covariance::from_observation_weights(string &pst_filename)
{
	ifstream ipst(pst_filename);
	pest_scenario.process_ctl_file(ipst, pst_filename);
	from_observation_weights(pest_scenario);

}

void Covariance::from_observation_weights(Pest &pest_scenario)
{

}

void Covariance::to_uncertainty_file(string &filename)
{

}
