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
#include <lapackpp.h>
#include <vector>
#include "SVD_PROPACK.h"

using namespace std;


extern "C" {
	double DLAMCH(char*);
	void DLANBPRO_SPARCE(int *m, int *n, int *k0, int *k, double *U, 
		int *ldu, double *V, int *ldv, double *B, int *ldb, 
		double *rnorm, double *doption, int *ioption, double *work,
		int *iwork, double *dparm, int *iparm, int *ierr);
}

SVD_PROPACK::SVD_PROPACK(void) : SVDPackage("PROPACK")
{
}

void SVD_PROPACK::solve_ip(LaGenMatDouble &A, LaVectorDouble &Sigma, LaGenMatDouble &U, LaGenMatDouble &Vt)
{
	class local_utils {
		public:
		static void init_array(double data[], int len, double value = 0.0)
		{
			for(int i=0; i<len; ++i)
			{
				data[i] = value;
			}
		}
	};

	int m_rows = A.rows();
	int n_cols = A.cols();
	int irow, icol;
	int n_nonzero;
	int k = 0;
	int kmax = min(m_rows, n_cols);
	kmax = min(n_max_sing, kmax);
	int ioption[] = {0, 1};
	char eps_char = 'e';
	double eps = DLAMCH(&eps_char);
	double d_option[] = {0.0, sqrt(eps), pow(eps, 3.0/4.0), 0.0};
	double rnorm = 0.0;
	int ierr = 0;

	Sigma = 0.0;
	U = 0.0;
	Vt = 0.0;
	//count number of nonzero entries
	for (n_nonzero=0, irow=0; irow<m_rows; ++irow)
	{
		for(icol=0; icol<n_cols; icol++)
		{
			if(A(irow, icol) != 0.0)
			{
				n_nonzero += 1;
			}
		}
	}
	// Allocate and initialize arrays
	double *dparm = new double[n_nonzero];
	int *iparm = new int[2*n_nonzero+1];
	double *tmp_u = new double[m_rows*(kmax+1)];
	double *tmp_v = new double[n_cols*kmax];
	double *tmp_b = new double[kmax*2];
	double *tmp_work = new double[2*(m_rows+n_cols+kmax+1)];
	int *tmp_iwork = new int[2*kmax+1];

	iparm[0] = n_nonzero;
	int i_nonzero;
	for (i_nonzero=0, irow=0; irow<m_rows; ++irow)
	{
		for(icol=0; icol<n_cols; icol++)
		{
			if(A(irow, icol) != 0.0)
			{
				dparm[i_nonzero] = A(irow, icol);
				iparm[i_nonzero+1] = (irow+1);
				iparm[n_nonzero+1+i_nonzero] = (icol+1);
				i_nonzero++;
			}
		}
	}
	local_utils::init_array(tmp_u, m_rows*(kmax+1), 0.0);
	local_utils::init_array(tmp_v, n_cols*kmax, 0.0);
	local_utils::init_array(tmp_b, kmax*2, 0.0);
	local_utils::init_array(tmp_work,2*(m_rows+n_cols+kmax+1), 0.0);

	int ld_tmpu = m_rows;
	int ld_tmpv = n_cols;
	int ld_tmpb = kmax;

	// Compute singluar values and vectors
	int k0 = -1;
    int k2 = 0;
    double eig_ratio = 1.0;
	while (k0<kmax-1 && eig_ratio > eign_thres)
	{
        k0 += 1;
        k2 = k0+1; // index of eignvlaue and eigenvector to be computed this time
		DLANBPRO_SPARCE(&m_rows, &n_cols, &k0, &k2, tmp_u, 
			&ld_tmpu, tmp_v, &ld_tmpv, tmp_b, &ld_tmpb, 
			&rnorm, d_option, ioption, tmp_work,
			tmp_iwork, dparm, iparm, &ierr);
        eig_ratio = tmp_b[k2-1] / tmp_b[0];
	}

	for (int i_sing=0; i_sing<kmax; ++i_sing)
	{
		Sigma(i_sing) = tmp_b[i_sing];
		for(int irow=0; irow<m_rows; ++ irow)
		{
			U(irow,i_sing) = tmp_u[i_sing*m_rows+irow];
		}
		for(int icol=0; icol<n_cols; ++ icol)
		{
			Vt(icol,i_sing) = tmp_v[i_sing*n_cols+icol];
		}
	}

		
	delete [] dparm;
	delete [] iparm;
	delete [] tmp_u;
	delete [] tmp_v;
	delete [] tmp_b;
	delete [] tmp_work;
    delete [] tmp_iwork;
}


void SVD_PROPACK::test()
{
	LaGenMatDouble A(3,3);
	LaVectorDouble Sigma(3);
	LaGenMatDouble U(3,3);
	LaGenMatDouble Vt(3,3);

	A = 0.0;
	A(0,0) = 20;
	A(1,1) = 30;
	A(2,2) = 50;

	set_max_sing(3);
	set_eign_thres(1e-7);


	solve_ip(A, Sigma, U, Vt);
	std::cout << "////////// PROPACK results /////////" << endl;
	std::cout << "A  = " << endl;
	std::cout << A << endl << endl;
	std::cout << "Sigma = " << endl;
	std::cout << Sigma << endl;
	std::cout << "U = " << endl;
	std::cout << U << endl << endl;
	std::cout << "Vt = " << endl;
	std::cout << Vt << endl << endl;

	LaSVD_IP(A, Sigma, U, Vt );
	std::cout << "////////// LAPACK results /////////" << endl;
	std::cout << "A  = " << endl;
	std::cout << A << endl << endl;
	std::cout << "Sigma = " << endl;
	std::cout << Sigma << endl;
	std::cout << "U = " << endl;
	std::cout << U << endl << endl;
	std::cout << "Vt = " << endl;
	std::cout << Vt << endl << endl;

}

SVD_PROPACK::~SVD_PROPACK(void)
{
}


