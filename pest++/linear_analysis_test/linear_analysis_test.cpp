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

map<string, double> get_obj_comps(string &filename)
{
	pest_utils::upper_ip(filename);
	ifstream ifile(filename);
	string line;
	
	vector<string> tokens;
	map<string, double> obj_comps;
	if (!ifile.good()) throw runtime_error("get_obj_comps error opening file: " + filename);
	if ((filename.find(".REI") != string::npos) || (filename.find(".RES") != string::npos))
	{
		double resid, weight;
		for (;;)
		{
			if (!getline(ifile, line)) throw runtime_error("get_obj_comps() error: EOF while looking for 'NAME'");
			pest_utils::upper_ip(line);
			if (line.find("NAME") != string::npos)
			{
				for (;;)
				{
					if (!getline(ifile, line)) break;
					pest_utils::upper_ip(line);
					tokens.clear();
					pest_utils::tokenize(line, tokens);
					pest_utils::convert_ip(tokens[4], resid);
					pest_utils::convert_ip(tokens[5], weight);
					if (obj_comps.find(tokens[1]) != obj_comps.end())obj_comps[tokens[1]] += pow(resid * weight, 2);
					else obj_comps[tokens[1]] = pow(resid * weight, 2);
				}
				break;
			}
		}

	}

	else if (filename.find(".REC") != string::npos)
	{
		throw runtime_error("get_obj_comps() .rec not implemented");
	}
	else if (filename.find(".IOBJ") != string::npos)
	{
		throw runtime_error("get_obj_comps() .iobj not implemented");
	}

	else throw runtime_error("get_obj_comps() error: unrecognized file type: " + filename + " must be .rei,.res,.rec,or .iobj");

	ifile.close();
	return obj_comps;

}

void normalize_weights_by_residual(Pest &pest_scenario, string &resid_filename)
{
	
	ObservationInfo obs_info = pest_scenario.get_ctl_observation_info();
	map<string, double> obj_comps = get_obj_comps(resid_filename);
	map<string, string> pst_grps = pest_scenario.get_observation_groups();
	vector<string> ogrp_names = pest_scenario.get_ctl_ordered_obs_group_names();
	const ObservationRec* obs_rec;
	double weight;

	map<string, int> grp_nnz;
	for (auto &ogrp : ogrp_names) grp_nnz[ogrp] = 0;

	for (auto &ogrp : pst_grps)
	{
		obs_rec = obs_info.get_observation_rec_ptr(ogrp.first);
		if (obs_rec->weight > 0.0) grp_nnz[ogrp.second]++;
	}

	for (auto &ogrp : pst_grps)
	{
		obs_rec = obs_info.get_observation_rec_ptr(ogrp.first);
		if (obs_rec->weight > 0.0)
		{
			weight = obs_info.get_observation_rec_ptr(ogrp.first)->weight * sqrt(((double)grp_nnz[ogrp.second]) / obj_comps[ogrp.second]);
			obs_info.set_weight(ogrp.first, weight);
		}
	}
}




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
		for (auto &pred : la.get_predictions())
		{
			pred.second.to_ascii(pred.first + ".vec");
		}
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

		map<string, double> first = la.first_prediction(200);
		map<string, double> second = la.second_prediction(200);
		map<string, double> third = la.third_prediction(200);
		first = la.first_prediction(200);
		second = la.second_prediction(200);
		third = la.third_prediction(200);
		first = la.first_prediction(300);
		second = la.second_prediction(300);
		third = la.third_prediction(300);
		first = la.first_prediction(599);
		second = la.second_prediction(599);
		third = la.third_prediction(599);
	}
	catch (exception &e)
	{
		cout << e.what() << endl;
		log.error(e.what());
	}
	log.log("analysis");
	return 0;
}

