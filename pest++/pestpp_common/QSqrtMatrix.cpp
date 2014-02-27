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
#include <utility>
#include <math.h>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include "QSqrtMatrix.h"
#include "Transformable.h"
#include "PriorInformation.h"
#include "pest_data_structs.h"


using namespace std;
using namespace Eigen;

QSqrtMatrix::QSqrtMatrix(const ObservationInfo &obs_info, const vector<string> &obs_names, 
	const PriorInformation *prior_info_ptr, double tikhonov_weight)
{
	weights.reserve(obs_names.size() + prior_info_ptr->size());
	unordered_map<string, ObservationRec>::const_iterator found_obsinfo_iter;
	unordered_map<string, ObservationRec>::const_iterator non_found_obsinfo_iter=obs_info.observations.end();
	PriorInformation::const_iterator found_prior_info;
	PriorInformation::const_iterator not_found_prior_info = prior_info_ptr->end();

	vector<string>::const_iterator b;
	vector<string>::const_iterator e;	

	int i;
	double weight = 0.0;

	for(i=0, b=obs_names.begin(), e=obs_names.end(); b != e; ++b, ++i)
	{
		found_obsinfo_iter = obs_info.observations.find(*b);
		found_prior_info = prior_info_ptr->find(*b);
		if (found_obsinfo_iter != non_found_obsinfo_iter)
		{
			weight = (*found_obsinfo_iter).second.weight;
			if( (*found_obsinfo_iter).second.is_regularization() ) {
				weight *= tikhonov_weight;
			}
			weights.push_back(sqrt(weight));
		}
		else if (found_prior_info != not_found_prior_info)
		{
			//calculate weights for prior information
			weight = (*found_prior_info).second.get_weight();
			if ((*found_prior_info).second.is_regularization()) {
				weight *= tikhonov_weight;
			}
			weights.push_back(sqrt(weight));
		}
		else {
			assert(true);  //observation not in standard observations or prior information 
		}
	}
}


bool ne_zero(double v) {return v!=0.0;}

//Eigen::SparseMatrix<double> QSqrtMatrix::operator*(const Eigen::SparseMatrix<double> &rhs) const
//{
//	Eigen::SparseMatrix<double> ret_val(rhs);
//	int n_rows = rhs.rows();
//
//	for (int i=0; i < n_rows; ++i) {
//		ret_val.row(i) *= weights[i];
//	}
//	return ret_val;
//}

int QSqrtMatrix::num_nonzero() const
{	
	int num = (int) count_if(weights.begin(), weights.end(), ne_zero);
	return num;
}

//Eigen::SparseMatrix<double> operator*(const Eigen::SparseMatrix<double> &lhs, const QSqrtMatrix &rhs)
//{
//	Eigen::SparseMatrix<double> ret_val(lhs);
//
//	int n_cols = lhs.cols();
//
//	for (int i=0; i < n_cols; ++i) {
//		ret_val.col(i) *= rhs.weights[i];
//	}
//	return ret_val;
//}


MatrixXd QSqrtMatrix::operator*(const MatrixXd &rhs) const
{
	MatrixXd ret_val(rhs);
	int n_rows = rhs.rows();

	for (int i=0; i < n_rows; ++i) {
		ret_val.row(i) *= weights[i];
	}
	return ret_val;
}

MatrixXd operator*(const MatrixXd &lhs, const QSqrtMatrix &rhs)
{
	MatrixXd ret_val(lhs);

	int n_cols = lhs.cols();

	for (int i=0; i < n_cols; ++i) {
		ret_val.col(i) *= rhs.weights[i];
	}
	return ret_val;
}



const VectorXd QSqrtMatrix::get_diag_vector() const
{
	Map<const VectorXd> vec(&weights[0], weights.size());
	return vec;
}


QSqrtMatrix::~QSqrtMatrix(void)
{
}
