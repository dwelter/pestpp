// linear_analysis_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include "Pest.h"
#include "covariance.h"


int main(int argc, char* argv[])
{
	//todos: schur: normalize weights by residuals; 
	//errvar: KL scaling, qhalf calc, R, I-R, G


	Pest pest_scenario;
	string pst_filename = "pest_1par.pst";

	Covariance parcov;
	parcov.from_parameter_bounds(pst_filename);

	Covariance obscov;
	obscov.from_observation_weights(pst_filename);

	Mat jacobian;
	jacobian.from_binary("pest_1par.jcb");
	jacobian.to_ascii("pest_1par_jcb.mat");

	//Eigen::SparseMatrix<double> prod = parcov_t.get_matrix() * jacobian.transpose().get_matrix() * parcov.get_matrix();
	//Mat result(parcov.get_row_names(), parcov.get_row_names(), prod, Mat::MatType::DENSE);
	vector<string> empty{}; 
	Mat pvec_t = jacobian.extract(vector<string>{"H02_08"},vector<string>() );
	Mat pvec = pvec_t.T();
 	Mat posterior = ((jacobian.T() * obscov.inv() * jacobian) + parcov.inv()).inv();
	Mat inv = posterior.inv();
	posterior.to_ascii("emu_test_result.mat");
	Mat prior_prd = pvec_t * parcov * pvec;
	cout << "prior" << prior_prd << endl;
	Mat post_prd = pvec_t * posterior * pvec;
	cout << "posterior" << post_prd << endl;



	/*
	Mat m;
	m.from_ascii("c_obs10_2");
	m.to_ascii("test.vec");*/
	//Covariance cov;
	//cov.from_parameter_bounds(pest_scenario);
	//vector<string> new_obs;
	//new_obs.push_back("H_OBS01_2");
	//new_obs.push_back("H_OBS01_1");
	//vector<string> new_par;
	//new_par.push_back("MULT1");
	//Covariance new_cov = cov.get(new_par);
	//cov.from_observation_weights(pest_scenario);
	//cov.to_ascii("test.cov");
	
	//Covariance cov;
	//cov.from_uncertainty_file("fake.unc");
	////cov.to_uncertainty_file("test.unc");

	//Mat mat1;
	//mat1.from_ascii("emu_mult_test1.mat");
	//
	//Mat mat2;
	//mat2.from_ascii("emu_mult_test2.mat");

	//Mat vec2 = mat2.get(vector<string>{"O1"}, mat2.get_col_names());
	//vec2.transpose();
	//Mat mat3 = mat1 * cov * mat2;
	//mat3.to_ascii("emu_mult_result.mat");





	return 0;
}

