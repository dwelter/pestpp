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
#ifndef LAPACK_TOOLS_H_
#define LAPACK_TOOLS_H_

#include <vector>
#include <ostream>
class Transformable;
class LaGenMatDouble;
class LaVectorDouble;

LaGenMatDouble diag_mat_mult(const LaVectorDouble &diag, const LaGenMatDouble &rhs);

LaGenMatDouble SVD_inv(const LaGenMatDouble &U, const LaVectorDouble &Sigma, 
					const LaGenMatDouble &Vt, int max_sing, double eigthresh, int &num_sing);

void get_LaGenMatDouble_row_abs_max(const LaGenMatDouble &m, int row, int *max_col, double *max_val);

LaVectorDouble stlvec2LaVec(const std::vector<double> &stl_vec);

void print(const LaGenMatDouble &mat, std::ostream &fout);

void add_LaVectorDouble_2_Transformable(Transformable &tr_data, const std::vector<std::string> &keys, 
										const LaVectorDouble &del_values);

#endif /* LAPACK_TOOLS_H_ */