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
	ipst.close();
	Mat m;
	m.from_ascii("c_obs10_2");
	return 0;
}

