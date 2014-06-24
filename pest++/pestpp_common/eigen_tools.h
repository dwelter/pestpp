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
#ifndef EIGEN_TOOLS_H_
#define EIGEN_TOOLS_H_

#include <vector>
#include <ostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Transformable;

void get_MatrixXd_row_abs_max(const Eigen::MatrixXd &m, int row, int *max_col, double *max_val);

Eigen::VectorXd stlvec_2_egienvec(const std::vector<double> &stl_vec);

Eigen::VectorXd transformable_2_egien_vec(const Transformable &data, std::vector<std::string> oredered_names);

void print(const Eigen::MatrixXd &mat, std::ostream &fout);

void matrix_del_cols(Eigen::SparseMatrix<double> &mat, const std::vector<int> &col_id_vec);
Eigen::SparseMatrix<double> get_diag_matrix(Eigen::SparseMatrix<double> &mat);

void print(const Eigen::MatrixXd &mat, std::ostream & fout, int n_per_line=7);
void print(const Eigen::VectorXd &vec, std::ostream & fout, int n_per_line=7);

#endif /* EIGEN_TOOLS_H_ */