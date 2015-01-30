// linear_analysis_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include "Pest.h"
#include "covariance.h"
#include "logger.h"
#include "utilities.h"
#include "linear_analysis.h"


int main(int argc, char* argv[])
{
	ofstream fout("linear_analysis.log");
	Logger log(fout);
	log.log("analysis");
	try
	{

		string jco("pest.jco");
		string pst("pest.pst");
		linear_analysis la(jco, pst, pst, &log);
		vector<string> preds;
		//preds.push_back("h02_08");
		//preds.push_back("c_obs01_1");
		preds.push_back("O10_1");
		la.set_predictions(preds);
		/*for (auto &pred : la.get_predictions())
		{
			pred.second.to_ascii(pred.first + ".vec");
		}*/
		vector<string> omitted;
		//omitted.push_back("mult1");
		omitted.push_back("k01_01_01");
		/*omitted.push_back("k01_01_02");
		omitted.push_back("k01_01_03");
		omitted.push_back("k01_01_04");*/
		la.extract_omitted(omitted);

		Covariance post = la.posterior_parameter_matrix();
		map<string, double> prpost = la.prior_prediction_variance();
		map<string, double> ptpost = la.posterior_prediction_variance();

		/*map<string, double> first = la.first_prediction(200);
		map<string, double> second = la.second_prediction(200);
		map<string, double> third = la.third_prediction(200);
		first = la.first_prediction(1);
		second = la.second_prediction(1);
		third = la.third_prediction(1);
		first = la.first_prediction(300);
		second = la.second_prediction(300);
		third = la.third_prediction(300);
		first = la.first_prediction(599);
		second = la.second_prediction(599);
		third = la.third_prediction(599);*/
	}
	catch (exception &e)
	{
		cout << e.what() << endl;
		log.error(e.what());
	}
	log.log("analysis");
	return 0;
}

