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
		
	string jco("pest.jco");
	string pst("pest.pst");
	linear_analysis la(jco,pst,pst);
	vector<string> preds;
	preds.push_back("h02_08");
	preds.push_back("h01_08");
	la.set_predictions(preds);
	la.get_jacobian_ptr()->to_ascii("jco.mat");
	la.get_obscov_ptr()->to_ascii("obscov.mat");
	//Mat* normal = la.get_normal_ptr();
	//Mat* R = la.get_R_ptr(1);
	//Mat* G = la.get_G_ptr(1);
	Mat first = la.first_parameter(1);
	Mat second = la.second_parameter(1);

	//map<string, double> pvar = la.prior_prediction_variance();
	//for (auto pred : pvar)
	//	cout << pred.first << ":" << pred.second << endl;
	////pvar = la.posterior_prediction_variance(string("c_obs01_1"));
	//pvar = la.posterior_prediction_variance();
	//for (auto pred : pvar)
	//	cout << pred.first << ":" << pred.second << endl;
	//vector<string> cond_pars;
	//cond_pars.push_back("kr01c03");
	//cond_pars.push_back("kr01c28");
	//map<string, pair<double, double>> contr_result = la.contribution(cond_pars);
	//for (auto r : contr_result)
	//	cout << r.first << "prior reduction:" << r.second.first << " post reduction:" << r.second.second << endl;
	//pvar = la.prior_prediction_variance();
	//for (auto pred : pvar)
	//	cout << pred.first << ":" << pred.second << endl;
	////pvar = la.posterior_prediction_variance(string("c_obs01_1"));
	//pvar = la.posterior_prediction_variance();
	//for (auto pred : pvar)
	//	cout << pred.first << ":" << pred.second << endl;
	//vector<string> obs;
	//obs.push_back("h_obs18_1");
	//obs.push_back("h_obs20_1");
	//map<string,double> worth_result = la.worth(obs);
	//for (auto r : worth_result)
	//{
	//	cout << r.first << ' ' << r.second << endl;
	//}


	return 0;
}

