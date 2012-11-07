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
#include <lapackpp.h>
#include <string>

SVDPackage::SVDPackage(std::string _descritpion, int _n_max_sing, double _eign_thres):n_max_sing(_n_max_sing), eign_thres(_eign_thres), description(_descritpion) {}


void SVDPackage::set_max_sing(int _n_max_sing) {
	n_max_sing = _n_max_sing;
}

void SVDPackage::set_eign_thres(double _eign_thres)
{
	eign_thres = _eign_thres;
}

void SVD_LAPACK::solve_ip(LaGenMatDouble& A, LaVectorDouble &Sigma, LaGenMatDouble& U, LaGenMatDouble& VT )
{
	LaSVD_IP(A, Sigma, U, VT );
}