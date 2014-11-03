// linear_analysis_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include "Pest.h"
#include "covariance.h"
#include "logger.h"


int main(int argc, char* argv[])
{
	//todos: schur: normalize weights by residuals; 
	//errvar: KL scaling, qhalf calc, R, I-R, G
	ofstream fout("emu.log");
	Logger log(fout,true);

	Pest pest_scenario;
	string pst_filename = "pest.pst";

	log.log("load parcov");
	Covariance parcov;
	parcov.from_parameter_bounds(pst_filename);
	log.log("load parcov");

	log.log("load obscov");
	Covariance obscov;
	obscov.from_observation_weights(pst_filename);
	log.log("load obscov");
	
	//vector<string> obs_names = obscov.get_row_names();
	//obs_names.erase(find(obs_names.begin(),obs_names.end(),string("O319_286")));
	//Covariance test = obscov.get(obs_names);

	log.log("load jco");
	Mat jacobian;
	jacobian.from_binary("pest.jco");
	//jacobian.to_ascii("pest_jco.mat");
	log.log("load jco");

	log.log("inv obscov");
	Mat obscov_inv = obscov.inv();
	log.log("inv obscov");

	log.log("jco_t");
	Mat jacobian_t = jacobian.T();
	log.log("jco_t");


	log.log("calc post1");
	Mat post1 = obscov_inv * jacobian;
	log.log("calc post1");
	log.log("calc post2");
	post1 = jacobian_t * post1;
	log.log("calc post2");
	log.log("calc schur");
	Mat posterior = (post1 + (parcov.inv())).inv();
	log.log("calc schur");
	log.log("write posterior");
	posterior.to_ascii("emu_test_result.mat");
	log.log("write posterior");


	log.log("get U");
	Mat U = jacobian.get_U();
	log.log("get U");
	log.log("get V");
	Mat V = jacobian.get_V();
	log.log("get V");
	log.log("get s");
	Mat s = jacobian.get_s();
	log.log("get s");
	fout.close();
	//cout << U;
	//cout << V;
	//cout << s;





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

