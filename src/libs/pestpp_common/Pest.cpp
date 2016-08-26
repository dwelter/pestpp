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

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include "Pest.h"
#include "utilities.h"
#include "pest_error.h"
#include <sstream>
#include "pest_data_structs.h"
#include "Transformation.h"
#include "FileManager.h"


using namespace::std;
using namespace::pest_utils;



Pest::Pest() : base_par_transform("PEST base_par_transform"), regul_scheme_ptr(0)
{
}

void Pest::set_defaults()
{
	svd_info.maxsing = 0;
	svd_info.eigthresh = 1.0e-6;
}

void Pest::check_inputs(ostream &f_rec)
{
	if ((control_info.noptmax == 0) && (pestpp_options.get_max_run_fail() > 1))
	{
		cout << "noptmax = 0, resetting max_run_fail = 1" << endl;
		f_rec << "noptmax = 0, resetting max_run_fail = 1" << endl;
		pestpp_options.set_max_run_fail(1);

	}

	Parameters numeric_pars = base_par_transform.ctl2model_cp(ctl_parameters);
	if (svd_info.maxsing == 0) {
		svd_info.maxsing = min(numeric_pars.size(), observation_values.size());
	}


	vector<string> par_warnings;
	vector<string> par_problems;
	bool unfixed_par = false;
	for (auto &pname : ctl_ordered_par_names)
	{
		//double pval = ctl_parameters[pname];
		//double lb = ctl_parameter_info.get_low_bnd(pname);
		const ParameterRec *prec = ctl_parameter_info.get_parameter_rec_ptr(pname);
		if (prec->tranform_type != ParameterRec::TRAN_TYPE::FIXED)
			unfixed_par = true;
		if (prec->init_value < prec->lbnd)
			par_problems.push_back(pname + " is less than lower bound");
		else if (prec->init_value == prec->lbnd)
			par_warnings.push_back(pname + " is at lower bound");
		if (prec->init_value > prec->ubnd)
			par_problems.push_back(pname + " is greater than upper bound");
		else if (prec->init_value == prec->ubnd)
			par_warnings.push_back(pname + " is at upper bound");
	}
	bool err = false;

	for (auto &str : par_warnings)
		cout << "parameter warning: " << str << endl;
	for (auto &str : par_problems)
	{
		cout << "parameter error: " << str << endl;
		err = true;
	}

	if (!unfixed_par)
	{
		cout << "parameter error: no adjustable parameters" << endl;
		err = true;
	}

	if (err)
		throw runtime_error("error in parameter data");

	int n_base = get_pestpp_options().get_n_iter_base();
	if (n_base == -1 || n_base > 0)
	{
	}
	else
	{
		stringstream ss;
		ss << "pest++ option 'n_iter_base' must either be -1 or greater than 0, not " << n_base;
		f_rec << "pest++ option 'n_iter_base' must either be -1 or greater than 0, not " << n_base;
		throw PestError(ss.str());
	}

	int n_super = get_pestpp_options().get_n_iter_super();
	if (n_super < 0)
	{
		stringstream ss;
		ss << "pest++ option 'n_iter_super' must >= 0, not " << n_super;
		f_rec << "pest++ option 'n_iter_super' must >= 0, not " << n_super;
		throw PestError(ss.str());
	}

	//check that prediction names are list in obs
	if ((pestpp_options.get_uncert_flag()) && (pestpp_options.get_prediction_names().size() > 0))
	{
		vector<string> missing;
		for (auto &pred_name : pestpp_options.get_prediction_names())
		{
			auto it_obs = find(ctl_ordered_obs_names.begin(), ctl_ordered_obs_names.end(), pred_name);
			if (it_obs == ctl_ordered_obs_names.end())
			{
				missing.push_back(pred_name);
			}
		}
		if (missing.size() > 0)
		{
			stringstream ss;
			ss << "Pest::check_inputs() the following predictions were not found in the observation names: ";
			for (auto &m : missing)
				ss << m << ',';
			f_rec << ss.str() << endl;
			throw PestError(ss.str());
		}
	}

	if (pestpp_options.get_auto_norm() > 0.0)
	{
		if (pestpp_options.get_parcov_scale_fac() > 0.0)
			throw PestError("Can't use 'autonorm' and 'parcov_scale_fac' > 0.0");
		f_rec << "pest++ option 'autonorm' is being deprecated in favor of 'use_parcov_scaling'" << endl;
		cout << "pest++ option 'autonorm' is being deprecated in favor of 'use_parcov_scaling'" << endl;
		//pestpp_options.set_auto_norm(-999.0);
	}
	if (pestpp_options.get_parcov_scale_fac() > 0.0)
	{
		if (pestpp_options.get_mat_inv() == PestppOptions::MAT_INV::Q12J)
		{
			throw PestError("pest++ mat_inv = q12j, but parcov_scale_fac > 0.0.");
		}
		if (pestpp_options.get_parcov_scale_fac() > 1.0)
		{
			cout << "'parcov_scale_fac' > 1.0, resetting to 1.0" << endl;
			pestpp_options.set_parcov_scale_fac(1.0);
		}

	}
}

void Pest::check_io()
{
	//make sure we can atleast access the model IO files
	vector<string> inaccessible_files;
	for (auto &file : model_exec_info.insfile_vec)
		if (!check_exist_in(file)) inaccessible_files.push_back(file);
	for (auto &file : model_exec_info.outfile_vec)
		if (!check_exist_out(file)) inaccessible_files.push_back(file);
	for (auto &file : model_exec_info.tplfile_vec)
		if (!check_exist_in(file)) inaccessible_files.push_back(file);
	for (auto &file : model_exec_info.inpfile_vec)
		if (!check_exist_out(file)) inaccessible_files.push_back(file);
	
	if (inaccessible_files.size() != 0)
	{
		string missing;
		for (auto &file : inaccessible_files)
			missing += file + " , ";

		throw PestError("Could not access the following model interface files: "+missing);
	}
}

const map<string, string> Pest::get_observation_groups() const
{
	map<string, string> obs_grp_map;
	for(auto &iobs : ctl_ordered_obs_names)
	{
		obs_grp_map[iobs] = (observation_info.get_observation_rec_ptr(iobs)->group);
	}
	return obs_grp_map;
}

vector<string> Pest::get_nonregul_obs() const
{
	vector<string>ret_val;
	for(Observations::const_iterator b=observation_values.begin(), e=observation_values.end();
		b!=e; ++b) {
		if (!observation_info.is_regularization((*b).first)) {
			ret_val.push_back((*b).first);
		}
	}
	return ret_val;
}

int Pest::process_ctl_file(ifstream &fin, string pst_filename)
{
	string line;
	string line_upper;
	string section("");
	vector<string> tokens;
	int sec_begin_lnum, sec_lnum;
	double value;
	string name;
	string *trans_type;
	string prior_info_string;
	pair<string, string> pi_name_group;
	int lnum;
	int num_par;
	int num_tpl_file;
	int i_tpl_ins = 0;
	double phimlim;
	double phimaccept;
	double fracphim;
	double wfinit;
	double wfmin;
	double wfmax;
	double wffac; 
	double wftol;
	bool use_dynamic_reg = false;
	bool reg_adj_grp_weights = false;
	vector<string> pestpp_input;

	regul_scheme_ptr = new DynamicRegularization(use_dynamic_reg);

	TranTied *t_tied = new TranTied("PEST to model tied transformation");
	TranOffset *t_offset = new TranOffset("PEST to model offset transformation");
	TranScale *t_scale = new TranScale("PEST to model scale transformation");
	TranLog10 *t_log = new TranLog10("PEST to model log transformation");
	TranFixed *t_fixed = new TranFixed("PEST to model fixed transformation");
	TranNormalize *t_auto_norm = new TranNormalize("PEST auto-normalization transformation");

	base_par_transform.push_back_ctl2model(t_scale);
	base_par_transform.push_back_ctl2model(t_offset);
	base_par_transform.push_back_ctl2active_ctl(t_tied);
	base_par_transform.push_back_ctl2active_ctl(t_fixed);
	base_par_transform.push_back_active_ctl2numeric(t_log);
	base_par_transform.push_back_active_ctl2numeric(t_auto_norm);

	try {
	prior_info_string = "";
	for(lnum=1, sec_begin_lnum=1; getline(fin, line); ++ lnum)
	{
		strip_ip(line);
		line_upper = upper_cp(line);
		tokens.clear();
		tokenize(line_upper, tokens);
		sec_lnum = lnum - sec_begin_lnum;
		if (tokens.empty())
		{
			//skip blank line
		}
		else if (line_upper.substr(0,2) == "++")
		{
			pestpp_input.push_back(line_upper);
		}
			
		else if (line_upper[0] == '*')
		{
			section = upper_cp(strip_cp(line_upper, "both", " *\t\n"));
			sec_begin_lnum = lnum;
		}
		else if (section == "CONTROL DATA")
		{
			if (sec_lnum == 1)
			{
				if (tokens[1] == "REGULARIZATION" || tokens[1] == "REGULARISATION")
				{
					use_dynamic_reg = true;
				}
			}

			if (sec_lnum == 2)
			{
				convert_ip(tokens[0], num_par);
			}
			if (sec_lnum == 3)
			{
				convert_ip(tokens[0], num_tpl_file);
				if(tokens.size() >= 5) {
					convert_ip(tokens[4], control_info.numcom);
				}
				else {
					control_info.numcom = 0;
				}
				if(tokens.size() >= 6) {
					convert_ip(tokens[5], control_info.jacfile);
				}
				else {
					control_info.jacfile = 0;
				}
			}
			else if (sec_lnum == 5)
			{
				convert_ip(tokens[0], control_info.relparmax);
				convert_ip(tokens[1], control_info.facparmax);
				convert_ip(tokens[2], control_info.facorig);
			}
			else if (sec_lnum == 6)
			{
				// remove text arguements from the line as these can be specified out of order
				// and PEST++ does not use them
				set<string> remove_tags = { "aui", "auid", "noaui", "senreuse", "nsenreuse", "boundscale", "noboundscale" };
				auto end_iter = std::remove_if(tokens.begin(), tokens.end(),
					[&remove_tags](string &str)->bool{return (remove_tags.find(upper_cp(str)) != remove_tags.end() 
					|| remove_tags.find(lower_cp(str)) != remove_tags.end()); });
				tokens.resize(std::distance(tokens.begin(), end_iter));

				convert_ip(tokens[0], control_info.phiredswh);
				if (tokens.size() >= 2) {
					convert_ip(tokens[1], control_info.noptswitch);
				}
				if (tokens.size() >= 3) {
					convert_ip(tokens[2], control_info.splitswh);
				}

			}
			else if (sec_lnum == 7)
			{
				convert_ip(tokens[0], control_info.noptmax);
				convert_ip(tokens[1], control_info.phiredstp);
				convert_ip(tokens[2], control_info.nphistp);
				convert_ip(tokens[3], control_info.nphinored);
				convert_ip(tokens[4], control_info.relparstp);
				convert_ip(tokens[5], control_info.nrelpar);
			}
			else {}
		}
		else if (section == "SINGULAR VALUE DECOMPOSITION")
		{
			if (sec_lnum == 2) {
				convert_ip(tokens[0], svd_info.maxsing);
				convert_ip(tokens[1], svd_info.eigthresh);
			}
			if (sec_lnum == 3) {
				convert_ip(tokens[0], svd_info.eigwrite);
			}
		}
		else if (section == "PARAMETER GROUPS")
		{
			ParameterGroupRec pgi;
			name = tokens[0];
			size_t n_tokens = tokens.size();
			pgi.name = name;
			ctl_ordered_par_group_names.push_back(name);
			convert_ip(tokens[1], pgi.inctyp);
			convert_ip(tokens[2], pgi.derinc);
			convert_ip(tokens[3], pgi.derinclb);
			convert_ip(tokens[4], pgi.forcen);
			convert_ip(tokens[5], pgi.derincmul);
			convert_ip(tokens[6], pgi.dermthd);
			if (n_tokens >= 8) convert_ip(tokens[7], pgi.splitthresh);
			if (n_tokens >= 9) convert_ip(tokens[8], pgi.splitreldiff);
			base_group_info.insert_group(name, pgi);
		}
		else if (section == "PARAMETER DATA")
		{
			if (sec_lnum <= num_par) {
				double scale;
				double offset;
				ParameterRec pi;
				name = tokens[0];
				trans_type = &tokens[1];
				convert_ip(tokens[2], pi.chglim);
				convert_ip(tokens[3], pi.init_value);
				convert_ip(tokens[4], pi.lbnd);
				convert_ip(tokens[5], pi.ubnd);
				convert_ip(tokens[6], pi.group);
				convert_ip(tokens[7], scale);
				convert_ip(tokens[8], offset);
				pi.scale = scale;
				pi.offset = offset;
				// add parameters to model parameter and paramter_info datasets
				ctl_ordered_par_names.push_back(name);
				if (*trans_type == "FIXED")
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::FIXED;
				}
				else if (*trans_type == "LOG")
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::LOG;
					n_adj_par++;
				}
				else if (*trans_type == "TIED")
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::TIED;
				}
				else if (*trans_type == "NONE")
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::NONE;
					n_adj_par++;
				}
				else
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::NONE;
					assert(true);
					n_adj_par++;
				}
				ctl_parameter_info.insert(name, pi);
				ctl_parameters.insert(name, pi.init_value);
				base_group_info.insert_parameter_link(name, pi.group);

				// build appropriate transformations
				if (*trans_type == "FIXED") {
					t_fixed->insert(name, pi.init_value);}
				else if (*trans_type == "LOG") {
					t_log->insert(name);
				}
				if (offset!=0) {
					t_offset->insert(name, offset);
				}
				if (scale !=1) {
					t_scale->insert(name, scale);
				}
			}
			// Get rest of information for tied paramters
			else {
				name = tokens[0];
				string name_tied =  tokens[1];
				double ratio =  ctl_parameters[name] / ctl_parameters[name_tied];
				t_tied->insert(name, pair<string, double>(name_tied, ratio));
			}
		}
		else if (section == "OBSERVATION GROUPS")
		{
			string name = tokens[0];
			ObservationGroupRec group_rec;
			observation_info.groups[name] = group_rec;
			vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), name);
			if (is == ctl_ordered_obs_group_names.end())
			{
				ctl_ordered_obs_group_names.push_back(name);
			}
		}
		else if (section == "OBSERVATION DATA")
		{
			ObservationRec obs_i;
			name = tokens[0];
			convert_ip(tokens[1], value);
			convert_ip(tokens[2], obs_i.weight);
			obs_i.group = tokens[3];
			ctl_ordered_obs_names.push_back(name);
			observation_info.observations[name] = obs_i;
			observation_values.insert(name, value);
		}

		else if (section == "PRIOR INFORMATION")
		{
			//This section processes the prior information.  It does not write out the
			//last prior infomration.  THis is because it must check for line continuations
			if (!prior_info_string.empty() && tokens[0] != "&"){
				pi_name_group = prior_info.AddRecord(prior_info_string);
				ctl_ordered_pi_names.push_back(pi_name_group.first);
				vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pi_name_group.second);
				if (is == ctl_ordered_obs_group_names.end())
				{
					ctl_ordered_obs_group_names.push_back(pi_name_group.second);
				}				
				prior_info_string.clear();
			}
			else if (tokens[0] == "&") {
				prior_info_string.append(" ");
			}
			prior_info_string.append(line_upper);
		}

		else if (section == "PRIOR INFORMATION" )
		{
			//This section processes the prior information.  It does not write out the
			//last prior infomration
			if (!prior_info_string.empty() && tokens[0] != "&") {
				pi_name_group = prior_info.AddRecord(prior_info_string);
				ctl_ordered_pi_names.push_back(pi_name_group.first);
				vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pi_name_group.second);
				if (is == ctl_ordered_obs_group_names.end())
				{
					ctl_ordered_obs_group_names.push_back(pi_name_group.second);
				}
				prior_info_string.clear();
			}
			else if (tokens[0] != "&") {
				prior_info_string.append(" ");
			}
			prior_info_string.append(line);
		}
		else if (section == "MODEL COMMAND LINE" )
		{
			model_exec_info.comline_vec.push_back(line);
		}
		else if (section == "MODEL INPUT/OUTPUT" )
		{
			vector<string> tokens_case_sen;
			tokenize(line, tokens_case_sen);
			if(i_tpl_ins < num_tpl_file)
			{
				model_exec_info.tplfile_vec.push_back(tokens_case_sen[0]);
				model_exec_info.inpfile_vec.push_back(tokens_case_sen[1]);
			}
			else
			{
				model_exec_info.insfile_vec.push_back(tokens_case_sen[0]);
				model_exec_info.outfile_vec.push_back(tokens_case_sen[1]);
			}
			++i_tpl_ins;
		}
		else if (section == "REGULARISATION" || section=="REGULARIZATION" )
		{
			if (sec_lnum == 1) {
				convert_ip(tokens[0], phimlim);
				convert_ip(tokens[1], phimaccept);
				fracphim = 0.0;
				if(tokens.size() >=3) convert_ip(tokens[2], fracphim);
			}
			else if (sec_lnum == 2) {
				convert_ip(tokens[0], wfinit);
				convert_ip(tokens[1], wfmin);
				convert_ip(tokens[2], wfmax);
			}
			else if (sec_lnum == 3) {
				int iregadj;
				convert_ip(tokens[0], wffac);
				convert_ip(tokens[1], wftol);
				if (tokens.size() > 2)
				{
					convert_ip(tokens[2], iregadj);
					if (iregadj == 1) reg_adj_grp_weights = true;
				}
				delete regul_scheme_ptr;
				regul_scheme_ptr = new DynamicRegularization(use_dynamic_reg, reg_adj_grp_weights, phimlim,
					phimaccept, fracphim, wfmin, wfmax, wffac, wftol, wfinit);
			}
		}
	}

	// write out last prior information record
	if (!prior_info_string.empty())
	{
		pi_name_group = prior_info.AddRecord(prior_info_string);
		ctl_ordered_pi_names.push_back(pi_name_group.first);
		vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pi_name_group.second);
		if (is == ctl_ordered_obs_group_names.end())
		{
			ctl_ordered_obs_group_names.push_back(pi_name_group.second);
		}
		prior_info_string.clear();
	}
	}
	catch (PestConversionError &e) {
		std::stringstream out;
		out << "Error parsing \"" << pst_filename << "\" on line number " << lnum << endl;
		out << e.what() << endl;
		e.add_front(out.str());
		e.raise();
	}
	fin.close();
	// process pest++ options last
	pestpp_options.set_n_iter_super(0);
	pestpp_options.set_n_iter_base(max(1, control_info.noptmax));
	pestpp_options.set_super_eigthres(svd_info.eigthresh);
	pestpp_options.set_max_n_super(ctl_parameters.size());
	pestpp_options.set_max_super_frz_iter(5);
	pestpp_options.set_max_reg_iter(20);
	pestpp_options.set_uncert_flag(true);
	pestpp_options.set_prediction_names(vector<string>());
	pestpp_options.set_parcov_filename(string());
	pestpp_options.set_basejac_filename(string());
	pestpp_options.set_sweep_parameter_csv_file(string());
	pestpp_options.set_sweep_output_csv_file("sweep_out.csv");
	pestpp_options.set_sweep_base_run(false);
	pestpp_options.set_sweep_forgive(false);
	pestpp_options.set_sweep_chunk(500);
	//pestpp_options.set_use_parcov_scaling(false);
	pestpp_options.set_parcov_scale_fac(-999.0);
	for(vector<string>::const_iterator b=pestpp_input.begin(),e=pestpp_input.end();
		b!=e; ++b) {
			pestpp_options.parce_line(*b);
	}

	if (pestpp_options.get_auto_norm() > 0.0)
	{
		cout << "WARNING 'autonorm' option is being deprecated in favor of use_parcov_scaling. " << endl;
		double u_bnd;
		double l_bnd;
		double avg;
		double spread;
		double auto_norm = pestpp_options.get_auto_norm();
		const string *par_name;
		Parameters upper_bnd;
		Parameters lower_bnd;
		for(vector<string>::const_iterator b=ctl_ordered_par_names.begin(), e=ctl_ordered_par_names.end();
			b!=e; ++b) {
				upper_bnd.insert(*b, ctl_parameter_info.get_parameter_rec_ptr(*b)->ubnd);
				lower_bnd.insert(*b, ctl_parameter_info.get_parameter_rec_ptr(*b)->lbnd);
		}
		base_par_transform.ctl2numeric_ip(upper_bnd);
		base_par_transform.ctl2numeric_ip(lower_bnd);
		for(Parameters::const_iterator b=upper_bnd.begin(), e=upper_bnd.end();
			b!=e; ++b) {
				par_name = &((*b).first);
				u_bnd = (*b).second;
				l_bnd = lower_bnd.get_rec(*par_name);
				spread = u_bnd - l_bnd;
				avg = (u_bnd + l_bnd) / 2.0;
				t_auto_norm->insert(*par_name, -avg, auto_norm/spread);
		}
	}

	regul_scheme_ptr->set_max_reg_iter(pestpp_options.get_max_reg_iter());
//	//Make sure we use Q1/2J is PROPACK is chosen
//	if (pestpp_options.get_svd_pack() == PestppOptions::SVD_PACK::PROPACK)
//	{
//		pestpp_options.set_mat_inv(PestppOptions::MAT_INV::Q12J);
//	}
	return 0;
}


const vector<string> &Pest::get_comline_vec()
{
	return model_exec_info.comline_vec;
}
const vector<string> &Pest::get_tplfile_vec()
{
	return model_exec_info.tplfile_vec;
}
const  vector<string> &Pest::get_inpfile_vec()
{
	return model_exec_info.inpfile_vec;
}
const  vector<string> &Pest::get_insfile_vec()
{
	return model_exec_info.insfile_vec;
}
const vector<string> &Pest::get_outfile_vec()
{
	return model_exec_info.outfile_vec;
}

Pest::~Pest() {
	//if (regul_scheme_ptr !=0) delete regul_scheme_ptr;
}

ostream& operator<< (ostream &os, const Pest& val)
{
	os << val.control_info;
	os << val.svd_info;
	return os;
}
