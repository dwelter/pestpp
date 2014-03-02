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
#include "eigen_tools.h"
#include <Eigen/Dense>
#include "Transformable.h"
#include <vector>
#include <algorithm>
#include <iomanip>


using namespace Eigen;
using namespace std;

MatrixXd diag_mat_mult(const VectorXd &diag, const MatrixXd &rhs)
{
	MatrixXd ret_val = diag.asDiagonal() * rhs;
	return ret_val;
}

Eigen::SparseMatrix<double> SVD_inv(const MatrixXd &U, const VectorXd &Sigma, 
					const MatrixXd &Vt, int max_sing, double eigthresh, int &num_sing)
{
	size_t s_size = Sigma.size();

	// Calculate V * S-1 * Ut
	// First Calculate S-1 
	num_sing = 0;
	for (int i=0; i < s_size; ++i) {
		if (Sigma(i) != 0 && i<max_sing && Sigma(i)/Sigma(0) > eigthresh) {
			++num_sing;
		}
	}
	Eigen::SparseMatrix<double> sigma_inv_trunc(num_sing, num_sing);
	sigma_inv_trunc.setZero();
	std::vector<Eigen::Triplet<double> > triplet_list;
	for (int i=0; i < num_sing; ++i) {
			triplet_list.push_back(Eigen::Triplet<double>(i,i,1.0 / Sigma(i)));
	}
	sigma_inv_trunc.setFromTriplets(triplet_list.begin(), triplet_list.end());

	// Calculate V * (S-1 * Ut) 
	Eigen::SparseMatrix<double> tmp_v = Vt.topRows(num_sing).transpose().sparseView();
	Eigen::SparseMatrix<double> tmp_ut =  U.leftCols(num_sing).transpose().sparseView();
	Eigen::SparseMatrix<double> ret_val = tmp_v * sigma_inv_trunc * tmp_ut;
	return ret_val;
}

void get_MatrixXd_row_abs_max(const MatrixXd &m, int row, int *max_col, double *max_val)
{
	size_t nrows = m.rows();
	size_t ncols = m.cols();
	assert(row <= nrows);

	*max_col = -999;
	*max_val = 0;
	for (size_t icol=0; icol<ncols; ++icol)
	{
		if(abs(m(row,icol)) > abs(*max_val)) 
		{
			*max_col = icol;
			*max_val = m(row,icol);
		}
	}
}

VectorXd stlvec_2_egienvec(const std::vector<double> &stl_vec)
{
	size_t len = stl_vec.size();
	VectorXd la_vec(len);
	for (size_t i=0; i<len; ++i)
	{
		la_vec(i) = stl_vec[i];
	}
	return la_vec;
}

void matrix_del_cols(MatrixXd &mat, const vector<int> &col_id_vec)
{
	size_t ncols = mat.cols();

	vector<int> del_id_vec = col_id_vec;
	std::sort(del_id_vec.begin(), del_id_vec.end(), [](int i, int j)->bool{return i>j;});

	// shift columns to over write deleted columns
	int n_shift = 0;
	for (int icol=0; icol<ncols; ++icol)
	{
		if (!del_id_vec.empty() && del_id_vec.back() == icol)
		{
			n_shift++;
			del_id_vec.pop_back();
		}
		else if (n_shift > 0)
		{
			mat.col(icol-n_shift) = mat.col(icol);
		}
	}
	//make sure all rows got removed;
	assert (del_id_vec.size() == 0);
	mat.conservativeResize(mat.rows(), ncols-n_shift);
}

void print(const MatrixXd &mat, ostream & fout)
{
	size_t nrows = mat.rows();
	size_t ncols = mat.cols();

	for (size_t i=0; i<nrows; ++i)
	{
		for (size_t j=0; j<ncols; ++j) 
		{
			fout << mat(i,j);
			if (j < ncols-1) {fout << ", ";}
		}
		fout << endl;
	}

}

void print(const MatrixXd &mat, ostream & fout, int n_per_line)
{
	size_t nrows = mat.rows();
	size_t ncols = mat.cols();
	size_t n = 0;

	for (size_t i=0; i<nrows; ++i)
	{
		for (size_t j=0; j<ncols; ++j) 
		{
			fout << setw(15) << setiosflags(ios::right) << mat(i,j);
			if ((j+1) % (n_per_line) == 0 || j+1==ncols)
			{
				fout << '\n';
			}
		}
		fout << endl;
	}

}

void print(const VectorXd &vec, ostream & fout, int n_per_line)
{
	size_t n = vec.size();

	for (size_t i=0; i<n; ++i)
	{
		fout << setw(15) << setiosflags(ios::right) << vec(i);
			if ((i+1) % (n_per_line) == 0 || i+1==n)
		{
			fout << '\n';
		}

	}
}


