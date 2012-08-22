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

void Pest::check_inputs()
{
	Parameters numeric_pars = base_par_transform.ctl2model_cp(ctl_parameters);
	if (svd_info.maxsing == 0) {
		svd_info.maxsing = min(numeric_pars.size(), observation_values.size());
	}
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


int Pest::process_ctl_file(const string &filename, FileManager &file_manager)
{
	ifstream fin;
	string line;
	string section("");
	vector<string> tokens;
	int sec_begin_lnum, sec_lnum;
	double value;
	string name;
	string *trans_type;
	string prior_info_string;
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
	vector<string> pestpp_input;

	regul_scheme_ptr = new Regularization();

	TranTied *t_tied = new TranTied("PEST to model tied transformation");
	TranOffset *t_offset = new TranOffset("PEST to model offset transformation");
	TranScale *t_scale = new TranScale("PEST to model scale transformation");
	TranLog10 *t_log = new TranLog10("PEST to model log transformation");
	TranFixed *t_fixed = new TranFixed("PEST to model fixed transformation");
	TranFrozen *t_frozen = new TranFrozen("PEST frozen parameter transformation");
	TranNormalize *t_auto_norm = new TranNormalize("PEST auto-normalization transformation");

	base_par_transform.push_back_ctl2model(t_scale);
	base_par_transform.set_scale_ptr(t_scale);
	base_par_transform.push_back_ctl2model(t_offset);
	base_par_transform.set_offset_ptr(t_offset);
	base_par_transform.push_back_ctl2numeric(t_tied);
	base_par_transform.push_back_ctl2numeric(t_fixed);
	base_par_transform.push_back_ctl2numeric(t_frozen);
	base_par_transform.set_frozen_ptr(t_frozen);
	base_par_transform.add_default_deep_copy(t_frozen);
	base_par_transform.push_back_ctl2numeric(t_log);
	base_par_transform.set_log10_ptr(t_log);
	base_par_transform.push_back_ctl2numeric(t_auto_norm);
	 
	try {
		fin.open(filename.c_str());
		if (! fin) throw(PestFileError(filename));
	}
	catch(...)
	{
		cerr << "Error: can not open control file: \"" << filename << "\"" << endl;
		throw(PestFileError(filename));
	}

	try {
	prior_info_string = "";
	for(lnum=1, sec_begin_lnum=1; getline(fin, line); ++ lnum)
	{
		strip_ip(line);
		upper_ip(line);
		tokens.clear();
		tokenize(line, tokens);
		sec_lnum = lnum - sec_begin_lnum;
		if (tokens.empty())
		{
			//skip blank line
		}
		else if (line.substr(0,2) == "++")
		{
			pestpp_input.push_back(line);
		}
			
		else if (line[0] == '*')
		{
			section = upper_cp(strip_cp(line, "both", " *\t\n"));
			sec_begin_lnum = lnum;
		}
		else if (section == "CONTROL DATA")
		{
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
				convert_ip(tokens[0], control_info.phiredswh);
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
		}
		else if (section == "PARAMETER GROUPS")
		{
			ParameterGroupRec pgi;
			name = tokens[0];
			convert_ip(tokens[1], pgi.inctyp);
			convert_ip(tokens[2], pgi.derinc);
			convert_ip(tokens[3], pgi.derinclb);
			convert_ip(tokens[4], pgi.forcen);
			convert_ip(tokens[5], pgi.derincmul);
			convert_ip(tokens[6], pgi.dermthd);
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
				// add parameters to model parameter and paramter_info datasets
				ctl_ordered_par_names.push_back(name);
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
		}
		else if (section == "OBSERVATION DATA")
		{
			ObservationRec obs_i;
			name = tokens[0];
			convert_ip(tokens[1], value);
			convert_ip(tokens[2], obs_i.weight);
			obs_i.weight *= obs_i.weight;
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
				prior_info.AddRecord(prior_info_string);
				prior_info_string.clear();
			}
			else if (tokens[0] == "&") {
				prior_info_string.append(" ");
			}
			prior_info_string.append(line);
		}

		else if (section == "DERIVATIVES COMMAND LINE")
		{
			if (sec_lnum == 1) {
			//command_line_derivates = line;
			}
			else if (sec_lnum == 2) {
				file_manager.set_analytic_derivative_filename(line);
			}
		}

		else if (section == "PRIOR INFORMATION" )
		{
			//This section processes the prior information.  It does not write out the
			//last prior infomration
			if (!prior_info_string.empty() && tokens[0] != "&") {
				prior_info.AddRecord(prior_info_string);
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
			if(i_tpl_ins < num_tpl_file)
			{
				model_exec_info.tplfile_vec.push_back(tokens[0]);
				model_exec_info.inpfile_vec.push_back(tokens[1]);
			}
			else
			{
				model_exec_info.insfile_vec.push_back(tokens[0]);
				model_exec_info.outfile_vec.push_back(tokens[1]);
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
				convert_ip(tokens[0], wffac);
				convert_ip(tokens[1], wftol);
				delete regul_scheme_ptr;
				regul_scheme_ptr = new RegularizationPest(phimlim, phimaccept, fracphim, wfmin, wfmax,
								   wffac, wftol);
			}
		}
	}

	// write out last prior information record
	if (!prior_info_string.empty())
	{
		prior_info.AddRecord(prior_info_string);
		prior_info_string.clear();
	}
	}
	catch (PestConversionError &e) {
		std::stringstream out;
		out << "Error parsing \"" << filename << "\" on line number " << lnum << endl;
		e.add_front(out.str());
		e.raise();
	}
	fin.close();
	// process pest++ options last
	pestpp_options.set_n_iter_super(0);
	pestpp_options.set_n_iter_base(control_info.noptmax);
	pestpp_options.set_super_eigthres(svd_info.eigthresh);
	pestpp_options.set_svd_pack(PestppOptions::LAPACK);
	pestpp_options.set_auto_norm(-999.0);
	for(vector<string>::const_iterator b=pestpp_input.begin(),e=pestpp_input.end();
		b!=e; ++b) {
			pestpp_options.parce_line(*b);
	}

	if (pestpp_options.get_auto_norm() > 0.0)
	{
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
	if (regul_scheme_ptr !=0) delete regul_scheme_ptr;
}

ostream& operator<< (ostream &os, const Pest& val)
{
	os << val.control_info;
	os << val.svd_info;
	return os;
}
