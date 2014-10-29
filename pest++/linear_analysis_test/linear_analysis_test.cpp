// linear_analysis_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include "Pest.h"
#include "covariance.h"


int main(int argc, char* argv[])
{
	Pest pest_scenario;
	string pst_filename = "pest.pst";
	ifstream ipst(pst_filename);
	pest_scenario.process_ctl_file(ipst, pst_filename);
	/*ipst.close();
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
	//cov.to_uncertainty_file("test.unc");

	Mat mat1;
	mat1.from_ascii("emu_mult_test1.mat");
	Mat mat2;
	mat2.from_ascii("emu_mult_test2.mat");
	Mat mat3 = mat1 * mat2;
	mat3.to_ascii("emu_mult_result.mat");


	return 0;
}

