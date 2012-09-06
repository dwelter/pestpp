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
#include <lapackpp.h>
#include <algorithm>
#include <vector>
#include "QSqrtMatrix.h"
#include "Transformable.h"
#include "PriorInformation.h"
#include "pest_data_structs.h"


using namespace std;

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

LaGenMatDouble QSqrtMatrix::operator*(const LaGenMatDouble &rhs) const
{
	LaGenMatDouble ret_val(rhs);
	LaGenMatDouble row;
	int n_rows = rhs.rows();

	for (int i=0; i < n_rows; ++i) {
		ret_val.row(i).scale(weights[i]);
	}
	return ret_val;
}

LaGenMatDouble  QSqrtMatrix::tran_q_mat_mult(const LaGenMatDouble &lhs) const
{
	LaGenMatDouble ret_val(lhs.cols(), lhs.rows());
	
	// Calculate Transpose of LHS
	for (int i=0, n_rows = lhs.rows(); i < n_rows; ++i) {
		for (int j=0, n_cols = lhs.cols(); j < n_cols; ++j) {
			ret_val(j,i) = lhs(i,j);
		}
	}
	// Calculate ret_val * Q
	for (int j=0, n_cols = ret_val.cols(); j < n_cols; ++j) {
		ret_val.col(j).scale(pow(weights[j], 2.0));
	}
	return ret_val;
}


bool ne_zero(double v) {return v!=0.0;}

int QSqrtMatrix::num_nonzero() const
{	
	int num = (int) count_if(weights.begin(), weights.end(), ne_zero);
	return num;
}

LaGenMatDouble operator*(const LaGenMatDouble &lhs, const QSqrtMatrix &rhs)
{
	LaGenMatDouble ret_val(lhs);
	LaGenMatDouble row;

	int n_cols = lhs.cols();

	for (int i=0; i < n_cols; ++i) {
		ret_val.col(i).scale(rhs.weights[i]);
	}
	return ret_val;
}



QSqrtMatrix::~QSqrtMatrix(void)
{
}
