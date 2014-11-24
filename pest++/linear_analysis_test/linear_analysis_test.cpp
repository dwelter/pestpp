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
	la.get_parcov().to_ascii("emu_prior.cov");
	Covariance schur = la.posterior_covariance_matrix();
	schur.to_ascii("emu_post.cov");
	return 0;
}

