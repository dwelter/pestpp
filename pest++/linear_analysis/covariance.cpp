
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

//---------------------------------------
//Mat constructors
//---------------------------------------
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


//-----------------------------------------
//Maninpulate the shape and ordering of Mats
//-----------------------------------------
void Mat::align(vector<string> &other_row_names, vector<string> &other_col_names)
{

}

Mat Mat::get(vector<string> &other_row_names, vector<string> &other_col_names)
{
	return Mat();
}

Mat Mat::extract(vector<string> &other_row_names, vector<string> &other_col_names)
{
	return Mat();
}

void Mat::drop(vector<string> &other_row_names, vector<string> &other_col_names)
{

}



//-----------------------------------------
//covariance matrices
//-----------------------------------------
Covariance::Covariance(vector<string> &names)
{
	row_names = names;
	icode = 1;
}

Covariance::Covariance()
{
	icode = 1;
}

void Covariance::from_uncertainty_file(const string &filename)
{

}

void Covariance::from_parameter_bounds(Pest &pest_scenario)
{
	vector<Eigen::Triplet<double>> triplet_list;
	const ParameterRec* par_rec;
	int i = 0;
	double upper, lower;
	for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
	{
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
	for (auto &obs_name : pest_scenario.get_ctl_ordered_obs_names())
	{
		obs_rec = pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(obs_name);
		if (obs_rec->weight > 0.0)
		{
			row_names.push_back(obs_name);
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
