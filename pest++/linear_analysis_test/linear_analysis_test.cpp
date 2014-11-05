// linear_analysis_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include "Pest.h"
#include "covariance.h"
#include "logger.h"
#include "utilities.h"

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
	//todos: schur: normalize weights by residuals; 
	/*ifstream ipst(string("pest.pst"));
	if (!ipst.good()) throw runtime_error("no good");
	Pest pest_scenario;
	pest_scenario.process_ctl_file(ipst, string("pest.pst"));
	normalize_weights_by_residual(pest_scenario, string("pest.rei"));
*/

	//errvar: KL scaling, qhalf calc, R, I-R, G
	ofstream fout("emu.log");
	Logger log(fout,true);
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

	log.log("calc schur2");
	Eigen::SparseMatrix<double> obscov_inv_mat = obscov.inv().get_matrix();
	Eigen::SparseMatrix<double> jacobian_mat = jacobian.get_matrix();
	Eigen::SparseMatrix<double> parcov_inv_mat = parcov.inv().get_matrix();
	Eigen::SparseMatrix<double> post1(parcov.nrow(), parcov.ncol());
	post1 = jacobian_mat.transpose() * obscov_inv_mat * jacobian_mat;
	Eigen::SparseMatrix<double> post_mat(parcov.nrow(), parcov.ncol());
	post_mat = parcov.get_matrix() - (post1 + parcov_inv_mat);
	log.log("calc schur2");

	/*log.log("calc schur");
	log.log("post1");
	Mat post1 = jacobian.T() * obscov_inv * jacobian;
	log.log("post1");
	Mat posterior = parcov - (post1 + (parcov.inv()));
	log.log("calc schur");*/

	


	/*log.log("write posterior");
	posterior.to_ascii("emu_test_result.mat");
	log.log("write posterior");*/


	/*log.log("get U");
	Mat U = jacobian.get_U();
	log.log("get U");
	log.log("get V");
	Mat V = jacobian.get_V();
	log.log("get V");
	log.log("get s");
	Mat s = jacobian.get_s();
	log.log("get s");*/
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

