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

#ifndef SVDPACKAGE_H_
#define SVDPACKAGE_H_
#include <string>
#include<Eigen/Dense>
#include<Eigen/Sparse>

class SVDPackage
{
public:
	SVDPackage(std::string _descritpion="undefined", int _n_max_sing=1000, double _eign_thres=1.0e-7);
		virtual void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::MatrixXd & U, Eigen::MatrixXd& VT ) = 0;
	virtual void set_max_sing(int _n_max_sing);
	virtual void set_eign_thres(double _eign_thres);
	virtual ~SVDPackage(void){};
	const std::string description;
protected:
	int n_max_sing;
	double eign_thres;
};

class SVD_EIGEN : public SVDPackage
{
public:
	SVD_EIGEN(void): SVDPackage("Eigen JacobiSVD")  {}
	virtual void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::MatrixXd& U, Eigen::MatrixXd& VT );
	virtual ~SVD_EIGEN(void) {}
};


#endif //SVDPACKAGE_H_
