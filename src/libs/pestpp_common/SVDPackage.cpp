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

#include "SVDPackage.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "RedSVD-h.h"


using namespace Eigen;

SVDPackage::SVDPackage(std::string _descritpion, int _n_max_sing, double _eign_thres) : description(_descritpion), n_max_sing(_n_max_sing), eign_thres(_eign_thres) 
{
	performance_log = nullptr;
}


void SVDPackage::set_performance_log(PerformanceLog *_performance_log)
{
	performance_log = _performance_log;
}

void SVDPackage::set_max_sing(int _n_max_sing) {
	n_max_sing = _n_max_sing;
}

void SVDPackage::set_eign_thres(double _eign_thres)
{
	eign_thres = _eign_thres;
}

int SVDPackage::get_max_sing() {
	return n_max_sing;
}

double SVDPackage::get_eign_thres()
{
	return eign_thres;
}

void SVD_REDSVD::solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
	Eigen::SparseMatrix<double>& Vt, Eigen::VectorXd &Sigma_trunc)
{
	solve_ip(A, Sigma, U, Vt, Sigma_trunc, eign_thres);
}


void SVD_REDSVD::solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
	Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc, double _eigen_thres)
{
	RedSVD::RedSVD<MatrixXd> red_svd(A,n_max_sing);

	U = red_svd.matrixU().sparseView();
	VT = red_svd.matrixV().transpose().sparseView();
	VectorXd Sigma_full = red_svd.singularValues();

	int kmax = (Sigma_full.size() < n_max_sing) ? Sigma_full.size() : n_max_sing;
	int num_sing_used = 0;
	double eig_ratio;
	for (int i_sing = 0; i_sing < kmax; ++i_sing)
	{
		eig_ratio = Sigma_full[i_sing] / Sigma_full[0];
		if (eig_ratio > _eigen_thres)
		{
			++num_sing_used;
		}
		else
		{
			break;
		}
	}

	Sigma = Sigma_full.head(num_sing_used);
	Sigma_trunc = Sigma_full.tail(Sigma_full.size() - num_sing_used);

	VT = Eigen::SparseMatrix<double>(VT.topRows(num_sing_used));
	U = Eigen::SparseMatrix<double>(U.leftCols(num_sing_used));

	return;

}


void SVD_EIGEN::solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
	Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc)
{
	solve_ip(A, Sigma, U, VT, Sigma_trunc, eign_thres);
}


void SVD_EIGEN::solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
	Eigen::SparseMatrix<double>& Vt, Eigen::VectorXd &Sigma_trunc,  double _eigen_thres)
{

	JacobiSVD<MatrixXd> svd_fac(A,  ComputeThinU |  ComputeThinV);
	VectorXd Sigma_full = svd_fac.singularValues();
	U = svd_fac.matrixU().sparseView();
	Vt = svd_fac.matrixV().transpose().sparseView();

	//Compute number of singular values to be used in the solution
	int num_sing_used = 0;
	double eig_ratio;

	int kmax = (Sigma_full.size() < n_max_sing) ? Sigma_full.size() : n_max_sing;
	for (int i_sing = 0; i_sing < kmax; ++i_sing)
	{
		eig_ratio = Sigma_full[i_sing] / Sigma_full[0];
		if (eig_ratio > _eigen_thres)
		{
			++num_sing_used;
		}
		else
		{
			break;
		}
	}
	//Trim the Matricies based on the number of singular values to be used
	Sigma = Sigma_full.head(num_sing_used);
	Sigma_trunc = Sigma_full.tail(Sigma_full.size() - num_sing_used);

	Vt = Eigen::SparseMatrix<double>(Vt.topRows(num_sing_used));
	U = Eigen::SparseMatrix<double>(U.leftCols(num_sing_used));
}