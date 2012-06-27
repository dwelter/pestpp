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
#include <vector>
#include "lapack_tools.h"
#include <lapackpp.h>

LaGenMatDouble diag_mat_mult(const LaVectorDouble &diag, const LaGenMatDouble &rhs)
{
	LaGenMatDouble ret_val(rhs);
	int n_rows = rhs.rows();

	for (int i=0; i < n_rows; ++i) {
		ret_val.row(i).scale(diag(i));
	}
	return ret_val;
}

LaGenMatDouble SVD_inv(const LaGenMatDouble &U, const LaVectorDouble &Sigma, 
					const LaGenMatDouble &Vt, int max_sing, double eigthresh, int &num_sing)
{
	int u_size = U.rows();
	int vt_size = Vt.rows();
	int s_size = Sigma.size();
	double scale_fac;

	LaGenMatDouble SinvUt(s_size, u_size);
	LaGenMatDouble ret_val(vt_size, u_size);

	// Calculate V * S-1 * Ut
	// First Calculate S-1 * Ut
	num_sing = 0;
	for (int i=0; i < s_size; ++i) {
		if (Sigma(i) != 0 && i<max_sing && Sigma(i)/Sigma(0) > eigthresh) {
			scale_fac = 1.0 / Sigma(i);
			++num_sing;
		}
		else {
			scale_fac = 0.0;
		}
		for (int j=0; j < u_size; ++j)
		SinvUt(i,j) = U(j,i) * scale_fac;
	}
	// Calculate V * (S-1 * Ut) 
	Blas_Mat_Trans_Mat_Mult(Vt(LaIndex(0,s_size-1), LaIndex(0,vt_size-1)), SinvUt, ret_val, 1.0, 0.0);
	return ret_val;
}

void get_LaGenMatDouble_row_abs_max(const LaGenMatDouble &m, int row, int *max_col, double *max_val)
{
	int nrows = m.rows();
	int ncols = m.cols();
	assert(row <= nrows);

	*max_col = -999;
	*max_val = 0;
	for (int icol=0; icol<ncols; ++icol)
	{
		if(abs(m(row,icol)) > abs(*max_val)) 
		{
			*max_col = icol;
			*max_val = m(row,icol);
		}
	}
}
LaVectorDouble stlvec2LaVec(const std::vector<double> &stl_vec)
{
	int len = stl_vec.size();
	LaVectorDouble la_vec(len);
	for (int i=0; i<len; ++i)
	{
		la_vec(i) = stl_vec[i];
	}
	return la_vec;
}