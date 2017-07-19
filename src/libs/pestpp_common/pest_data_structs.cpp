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
#include <iomanip>
#include <vector>
#include <map>
#include <sstream>
#include <list>
#include <regex>
#include "pest_data_structs.h"
#include "utilities.h"
#include "pest_error.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "pest_error.h"
#include "Transformable.h"

using namespace::std;
using namespace::pest_utils;


ostream& operator<< (ostream &os, const ControlInfo& val)
{
	os << "PEST Control Information" << endl;
	os << "    relparmax = " << val.relparmax << endl;
	os << "    facparmax = " << val.facparmax << endl;
	os << "    facorig = " << val.facorig << endl;
	os << "    phiredswh = " << val.phiredswh << endl;
	os << "    noptmax = " << val.noptmax << endl;
	os << "    phiredstp = " << val.phiredstp << endl;
	os << "    nphistp = " << val.nphistp << endl;
	os << "    nphinored = " << val.nphinored << endl;
	os << "    relparstp = " << val.relparstp << endl;
	os << "    nrelpar = " << val.nrelpar << endl;
	return os;
}

//////////////////  ParameterGroupRec Methods ////////////////////////////////

ParameterGroupRec& ParameterGroupRec::operator=(const ParameterGroupRec &rhs)
{
	name = rhs.name;
	inctyp = rhs.inctyp;
	derinc = rhs.derinc;
	derinclb = rhs.derinclb;
	forcen = rhs.forcen;
	derincmul = rhs.derincmul;
	dermthd = rhs.dermthd;
	splitthresh = rhs.splitthresh;
	splitreldiff = rhs.splitreldiff;
	return *this;
}

ostream& operator<< (ostream &os, const ParameterGroupRec& val)
{
	os << "PEST Parameter Group Information" << endl;
	os << "    name   = " << val.name << endl;
	os << "    inctyp = " << val.inctyp << endl;
	os << "    derinc = " << val.derinc << endl;
	os << "    derinclb = " << val.derinclb << endl;
	os << "    forcen = " << val.forcen << endl;
	os << "    derincmul = " << val.derincmul << endl;
	return os;
}

/////////////////////////////////////  ParameterGroupInfo Methods ///////////////////////
const ParameterGroupRec* ParameterGroupInfo::get_group_rec_ptr(const string &name) const
{
	const ParameterGroupRec *ret_val = 0;
	unordered_map<string, ParameterGroupRec*>::const_iterator g_iter;

	g_iter = parameter2group.find(name);
	if(g_iter != parameter2group.end()) {
		ret_val = (*g_iter).second;
	}
	return ret_val;
}

ParameterGroupRec* ParameterGroupInfo::get_group_rec_ptr_4_mod(const string &name)
{
	ParameterGroupRec *ret_val = 0;
	unordered_map<string, ParameterGroupRec*>::const_iterator g_iter;

	g_iter = parameter2group.find(name);
	if (g_iter != parameter2group.end()) {
		ret_val = (*g_iter).second;
	}
	return ret_val;
}

string ParameterGroupInfo::get_group_name(const string &par_name) const
{
	return get_group_rec_ptr(par_name)->name;
}

void ParameterGroupInfo::insert_group(const string &group_name, ParameterGroupRec &rec) 
{
	groups[group_name] = new ParameterGroupRec(rec);
}

void ParameterGroupInfo::insert_parameter_link(const string &parameter_name, const string & group_name)
{
	unordered_map<string, ParameterGroupRec*>::const_iterator g_iter;

	g_iter = groups.find(group_name);
	if(g_iter == groups.end()) {
		throw PestIndexError(group_name, "Invalid parameter group name");
	}
	else {
		parameter2group[parameter_name] = (*g_iter).second;
	}

}

const ParameterGroupInfo& ParameterGroupInfo::operator=(const ParameterGroupInfo &rhs)
{
	unordered_map<ParameterGroupRec*, ParameterGroupRec*> old2new;
	unordered_map<string, ParameterGroupRec*>::const_iterator it(rhs.groups.begin());
	unordered_map<string, ParameterGroupRec*>::const_iterator end(rhs.groups.end());
	for (; it != end; ++it) {
		ParameterGroupRec* new_ptr = new ParameterGroupRec(*(*it).second);
		groups[(*it).first] = new_ptr;
		old2new[(*it).second] = new_ptr;
	}
	unordered_map<ParameterGroupRec*, ParameterGroupRec*>::iterator it_find;
	it = rhs.parameter2group.begin();
	end =  rhs.parameter2group.end();
	for (; it != end; ++it) {
		it_find = old2new.find((*it).second);
		parameter2group[(*it).first] = (*it_find).second;
	}
	return *this;
}

vector<string> ParameterGroupInfo::get_group_names() const
{
	vector<string> group_names;
	for (auto &g : groups)
		group_names.push_back(g.first);
	return group_names;
}

bool ParameterGroupInfo::have_switch_derivative() const
{
	bool switch_der = false;
	for (const auto irec : groups)
	{
		if (lower_cp(irec.second->forcen) == "switch")
		{
			switch_der = true;
			break;
		}
	}
	return switch_der;
}

ParameterGroupInfo::~ParameterGroupInfo()
{
	unordered_map<string, ParameterGroupRec*>::iterator it(groups.begin());
	unordered_map<string, ParameterGroupRec*>::iterator end(groups.end());
	for (; it != end; ++it) {
		delete (*it).second;
	}
}

ostream& operator<< (ostream &os, const ParameterGroupInfo &val)
{
	unordered_map<string, ParameterGroupRec*>::const_iterator it(val.parameter2group.begin());
	unordered_map<string, ParameterGroupRec*>::const_iterator end(val.parameter2group.end());
	for (; it != end; ++it) {
		os << *(it->second) << endl;
	}
	return os;
}


/////////////////////////////////////  ParameterRec Methods ///////////////////////

ostream& operator<< (ostream &os, const ParameterRec& val)
{
	os << "PEST Parameter Information" << endl;
	os << "    chglim = " << val.chglim << endl;
	os << "    lbnd = " << val.lbnd << endl;
	os << "    ubnd = " << val.ubnd << endl;
	os << "    group = " << val.group << endl;
	os << "    dercom = " << val.dercom << endl;
	return os;
}

/////////////////////////////////////  CtlParameterInfo Methods ///////////////////////



const ParameterRec* ParameterInfo::get_parameter_rec_ptr(const string &name) const
{
	const ParameterRec *ret_val = 0;
	unordered_map<string, ParameterRec>::const_iterator p_iter;

	p_iter = parameter_info.find(name);
	if(p_iter != parameter_info.end()) {
		ret_val = &((*p_iter).second);
	}
	return ret_val;
}

Parameters ParameterInfo::get_low_bnd(const vector<string> &keys) const
{
	Parameters l_bnd;
	vector<string>::const_iterator iend=keys.end();
	ParameterRec const *v_ptr;
	for (vector<string>::const_iterator i=keys.begin(); i!=iend; ++i)
	{
		v_ptr = get_parameter_rec_ptr(*i);
		if (v_ptr) {
			l_bnd.insert(*i, v_ptr->lbnd);
		}
		else {
			l_bnd.insert(*i, Parameters::no_data);
		}
	}
	return l_bnd;
}

Parameters ParameterInfo::get_up_bnd(const vector<string> &keys) const
{
	Parameters u_bnd;
	vector<string>::const_iterator iend=keys.end();
	ParameterRec const *v_ptr;
	for (vector<string>::const_iterator i=keys.begin(); i!=iend; ++i)
	{
		v_ptr = get_parameter_rec_ptr(*i);
		if (v_ptr) {
			u_bnd.insert(*i, v_ptr->ubnd);
		}
		else {
			u_bnd.insert(*i, Parameters::no_data);
		}
	}
	return u_bnd;
}

Parameters ParameterInfo::get_init_value(const vector<string> &keys) const
{
	Parameters init_value;
	vector<string>::const_iterator iend=keys.end();
	ParameterRec const *v_ptr;
	for (vector<string>::const_iterator i=keys.begin(); i!=iend; ++i)
	{
		v_ptr = get_parameter_rec_ptr(*i);
		if (v_ptr) {
			init_value.insert(*i, v_ptr->init_value);
		}
		else {
			init_value.insert(*i, Parameters::no_data);
		}
	}
	return init_value;
}


ostream& operator<< (ostream &os, const PestppOptions& val)
{
	os << "PEST++ Options" << endl;
	os << "    n_iter_base = " << left << setw(20) << val.get_n_iter_base() << endl;
	os << "    n_iter_super = " << left << setw(20) << val.get_n_iter_super() << endl;
	os << "    max_n_super = " << left << setw(20) << val.get_max_n_super() << endl;
	os << "    super eigthres = " << left << setw(20) << val.get_super_eigthres() << endl;
	os << "    svd pack = " << left << setw(20) << val.get_svd_pack() << endl;
	os << "    auto norm = " << left << setw(20) << val.get_auto_norm() << endl;
	os << "    super relparmax = " << left << setw(20) << val.get_super_relparmax() << endl;
	os << "    max super frz iter = " << left << setw(20) << val.get_max_super_frz_iter() << endl;
	os << "    mat inv = " << left << setw(20) << val.get_mat_inv() << endl;
	os << "    max run fail = " << left << setw(20) << val.get_max_run_fail() << endl;
	os << "    max reg iter = " << left << setw(20) << val.get_max_reg_iter() << endl;	
	os << "    use jacobian scaling a la PEST? = ";
	if (val.get_jac_scale())
		os << " yes" << endl;
	else
		os << " no" << endl;
	if (val.get_reg_frac() > 0.0)
		os << "    regularization fraction of total phi = " << left << setw(10) << val.get_reg_frac() << endl;
	os << "    lambdas = " << endl;
	for (auto &lam : val.get_base_lambda_vec())
	{
		os << right << setw(15) << lam << endl;
	}
	os << "    lambda scaling factors = " << endl;
	for (auto &ls : val.get_lambda_scale_vec())
		os << right << setw(15) << ls << endl;

	if (!val.get_basejac_filename().empty())
	{
		os << "   restarting with existing jacobian matrix file: " << val.get_basejac_filename() << endl;
		if (!val.get_hotstart_resfile().empty())
		{
			os << "   and using existing residual file " << val.get_hotstart_resfile() << " to forego initial model run" << endl;
		}
	}
	if (val.get_uncert_flag())
	{
		os << "    using FOSM-based uncertainty estimation for parameters" << endl;
		os << "    parameter covariance file = " << left << setw(20) << val.get_parcov_filename() << endl;
		if (val.get_prediction_names().size() > 0)
		{
			os << "    using FOSM-based uncertainty for forecasts" << endl;
			os << "    forecast names = " << endl;
		}
		
	}
	for (auto &pname : val.get_prediction_names())
		os << right << setw(15) << pname << endl;
	os << "    derivative run failure forgive = " << left << setw(15) << val.get_der_forgive() << endl;
	os << "    run overdue reschedule factor = " << left << setw(20) << val.get_overdue_reched_fac() << endl;
	os << "    run overdue giveup factor = " << left << setw(20) << val.get_overdue_giveup_fac() << endl;
	os << "    base parameter jacobian filename = " << left << setw(20) << val.get_basejac_filename() << endl;
	os << "    prior parameter covariance upgrade scaling factor = " << left << setw(10) << val.get_parcov_scale_fac() << endl;
	if (val.get_global_opt() == PestppOptions::GLOBAL_OPT::OPT_DE)
	{
		os << "    global optimizer = differential evolution (DE)" << endl;
		os << "    DE CR = " << left << setw(10) << val.get_de_cr() << endl;
		os << "    DE F = " << left << setw(10) << val.get_de_f() << endl;
		os << "    DE population size = " << setw(10) << val.get_de_npopulation() << endl;
		os << "    DE max generations = " << setw(10) << val.get_de_max_gen() << endl;
		os << "    DE F dither = " << left << setw(10) << val.get_de_dither_f() << endl;
	}
	os << endl;
	return os;
}

PestppOptions::PestppOptions(int _n_iter_base, int _n_iter_super, int _max_n_super, double _super_eigthres,
	SVD_PACK _svd_pack, MAT_INV _mat_inv, double _auto_norm, double _super_relparmax, int _max_run_fail,
	bool _iter_summary_flag, bool _der_forgive, double _overdue_reched_fac, double _overdue_giveup_fac, double _reg_frac,
	GLOBAL_OPT _global_opt, double _de_f, double _de_cr, int _de_npopulation, int _de_max_gen, bool _de_dither_f)
	: n_iter_base(_n_iter_base), n_iter_super(_n_iter_super), max_n_super(_max_n_super), super_eigthres(_super_eigthres), 
	svd_pack(_svd_pack), mat_inv(_mat_inv), auto_norm(_auto_norm), super_relparmax(_super_relparmax),
	max_run_fail(_max_run_fail), max_super_frz_iter(50), max_reg_iter(50), base_lambda_vec({ 0.1, 1.0, 10.0, 100.0, 1000.0 }),
	lambda_scale_vec({0.9, 0.8, 0.7, 0.5}),
	iter_summary_flag(_iter_summary_flag), der_forgive(_der_forgive), overdue_reched_fac(_overdue_reched_fac),
	overdue_giveup_fac(_overdue_giveup_fac), reg_frac(_reg_frac), global_opt(_global_opt),
	de_f(_de_f), de_cr(_de_cr), de_npopulation(_de_npopulation), de_max_gen(_de_max_gen), de_dither_f(_de_dither_f)
{
}

void PestppOptions::parce_line(const string &line)
{
	string key;
	string value;
	regex lambda_reg("(\\w+)(?:\\s*\\()([^\\)]+)(?:\\))");
	const std::sregex_iterator end_reg;
	//regex lambda_reg("((\\w+\\s*\\([^\\)]+\\))+?)");
	cmatch mr;

	size_t found = line.find_first_of("#");
	if (found == string::npos) {
		found = line.length();
	}
	string tmp_line = line.substr(0, found);
	strip_ip(tmp_line, "both", "\t\n\r+ ");
	//upper_ip(tmp_line);


	for (std::sregex_iterator i(tmp_line.begin(), tmp_line.end(), lambda_reg); i != end_reg; ++i)
	{
		string key = (*i)[1];
		string org_value = (*i)[2];
		upper_ip(key);
		string value = upper_cp(org_value);

		if (key=="MAX_N_SUPER"){
			convert_ip(value, max_n_super); 
		}
		else if (key=="SUPER_EIGTHRES"){
			convert_ip(value, super_eigthres); 
		}
		else if (key=="N_ITER_BASE"){
			convert_ip(value, n_iter_base); 
		}
		else if (key=="N_ITER_SUPER"){
			convert_ip(value, n_iter_super); 
		}
		else if (key=="SVD_PACK"){

			if (value == "PROPACK")
				svd_pack = PROPACK;
			else if (value == "REDSVD")
				svd_pack = REDSVD;
			else if ((value == "EIGEN") || (value == "JACOBI"))
				svd_pack == EIGEN;
			else
				throw PestParsingError(line, "Invalid ++svd_pack: \"" + value + "\"");


		}
		else if (key=="AUTO_NORM"){
			convert_ip(value, auto_norm); 
		}
		else if (key == "SUPER_RELPARMAX"){
			convert_ip(value, super_relparmax);
		}
		else if (key == "MAX_SUPER_FRZ_ITER"){
			convert_ip(value, max_super_frz_iter);
		}
		else if (key == "MAT_INV"){
			if (value == "Q1/2J") mat_inv = Q12J;
		}
		else if (key == "MAX_RUN_FAIL"){
			convert_ip(value, max_run_fail);
		}
		else if (key == "MAX_REG_ITER"){
			convert_ip(value, max_reg_iter);
		}
		else if (key == "LAMBDAS")
		{
			base_lambda_vec.clear();
			vector<string> lambda_tok;
			tokenize(value, lambda_tok, ",");
			for (const auto &ilambda : lambda_tok)
			{
				base_lambda_vec.push_back(convert_cp<double>(ilambda));
			}
		}	
		else if (key == "LAMBDA_SCALE_FAC")
		{
			lambda_scale_vec.clear();
			vector<string> scale_tok;
			tokenize(value, scale_tok, ",");
			for (const auto &iscale : scale_tok)
			{
				lambda_scale_vec.push_back(convert_cp<double>(iscale));
			}
		}
		else if (key == "ITERATION_SUMMARY")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> iter_summary_flag;
		}
		else if (key == "DER_FORGIVE")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> der_forgive;
		}
		else if (key == "UNCERTAINTY")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> uncert;
		}
		else if (key == "PREDICTIONS" || key == "FORECASTS")
		{
			prediction_names.clear();
			vector<string> prediction_tok;
			tokenize(value, prediction_tok, ",");
			for (const auto &pname : prediction_tok)
			{
				prediction_names.push_back(pname);
			}
		}
		else if ((key == "PARCOV") || (key == "PARAMETER_COVARIANCE") 
			|| (key == "PARCOV_FILENAME"))
		{
			convert_ip(org_value, parcov_filename);
		}

		else if ((key == "BASE_JACOBIAN") || (key == "BASE_JACOBIAN_FILENAME"))
		{
			convert_ip(org_value, basejac_filename);
		}

		else if (key == "HOTSTART_RESFILE")
		{
			convert_ip(org_value, hotstart_resfile);
		}

		else if (key == "OVERDUE_RESCHED_FAC"){
			convert_ip(value, overdue_reched_fac);
		}
		else if (key == "OVERDUE_GIVEUP_FAC"){
			convert_ip(value, overdue_giveup_fac);
		}
		else if (key == "CONDOR_SUBMIT_FILE")
		{
			convert_ip(value, condor_submit_file);
		}
		else if (key == "SWEEP_PARAMETER_CSV_FILE")
			convert_ip(org_value, sweep_parameter_csv_file);
		
		else if (key == "SWEEP_OUTPUT_CSV_FILE")
			convert_ip(org_value, sweep_output_csv_file);
		else if (key == "SWEEP_CHUNK")
			convert_ip(value, sweep_chunk);
		else if (key == "SWEEP_FORGIVE")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> sweep_forgive;
		}
		else if (key == "SWEEP_BASE_RUN")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> sweep_base_run;
		}
		else if (key == "REG_FRAC")
		{
			convert_ip(value, reg_frac);
		}
		/*else if (key == "USE_PARCOV_SCALING")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> use_parcov_scaling;
		}*/
		else if (key == "PARCOV_SCALE_FAC")
		{
			convert_ip(value, parcov_scale_fac);
		}
		else if (key == "JAC_SCALE")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> jac_scale;

		}

		else if (key == "UPGRADE_AUGMENT")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> upgrade_augment;

		}

		else if (key == "UPGRADE_BOUNDS")
		{
			if (value == "ROBUST")
				convert_ip(value, upgrade_bounds);
			else if (value == "CHEAP")
				convert_ip(value, upgrade_bounds);
			else
				throw runtime_error("unrecognozed 'upgrade_bounds' option: should 'robust' or 'cheap'");
			
		}

		else if (key == "GLOBAL_OPT")
		{
			if (value == "DE") global_opt = OPT_DE;
		}
		else if (key == "DE_F")
		{
			convert_ip(value, de_f);
		}
		else if (key == "DE_CR")
		{
			convert_ip(value, de_cr);
		}
		else if (key == "DE_POP_SIZE")
		{
			convert_ip(value, de_npopulation);
		}
		else if (key == "DE_MAX_GEN")
		{
			convert_ip(value, de_max_gen);
		}
		else if (key == "DE_DITHER_F")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> de_dither_f;
		}
		else if ((key == "OPT_OBJ_FUNC") || (key == "OPT_OBJECTIVE_FUNCTION"))
		{
			convert_ip(value,opt_obj_func);
		}
		else if (key == "OPT_COIN_LOG")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> opt_coin_log;
		}
		else if (key == "OPT_SKIP_FINAL")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> opt_skip_final;
		}

		else if ((key == "OPT_DEC_VAR_GROUPS") || (key == "OPT_DECISION_VARIABLE_GROUPS"))
		{
			opt_dec_var_groups.clear();
			vector<string> tok;
			tokenize(value, tok, ", ");
			for (const auto &name : tok)
			{
				opt_dec_var_groups.push_back(name);
			}
		}

		else if ((key == "OPT_EXT_VAR_GROUPS") || (key == "OPT_EXTERNAL_VARIABLE_GROUPS"))
		{
			opt_external_var_groups.clear();
			vector<string> tok;
			tokenize(value, tok, ", ");
			for (const auto &name : tok)
			{
				opt_external_var_groups.push_back(name);
			}
		}

		else if ((key == "OPT_CONSTRAINT_GROUPS"))
		{
			opt_constraint_groups.clear();
			vector<string> tok;
			tokenize(value, tok, ", ");
			for (const auto &name : tok)
			{
				opt_constraint_groups.push_back(name);
			}
		}
		
		else if (key == "OPT_RISK")
		{
			convert_ip(value, opt_risk);
		}

		else if (key == "OPT_ITER_DERINC_FAC")
		{
			convert_ip(value, opt_iter_derinc_fac);
		}

		else if (key == "OPT_DIRECTION")
		{
			string v;
			convert_ip(value,v);
			if (v == "MAX")
				opt_direction = -1;
			else if (v == "MIN")
				opt_direction = 1;
			else
				throw runtime_error("++opt_direction arg must be in {MAX,MIN}, not " + v);
		}

		else if (key == "OPT_ITER_TOL")
		{
			convert_ip(value, opt_iter_tol);
		}

		else if (key == "OPT_RECALC_FOSM_EVERY")
		{
			convert_ip(value, opt_recalc_fosm_every);
		}
		else if ((key == "IES_PAR_CSV") || (key == "IES_PARAMETER_CSV"))
		{
			convert_ip(value, ies_par_csv);
		}
		else if ((key == "IES_OBS_CSV") || (key == "IES_OBSERVATION_CSV"))
		{
			convert_ip(value, ies_obs_csv);
		}
		else if ((key == "IES_OBS_RESTART_CSV") || (key == "IES_OBSERVATION_RESTART_CSV"))
		{
			convert_ip(value, ies_obs_restart_csv);
		}

		else if ((key == "IES_USE_APPROXIMATE_SOLUTION") || (key == "IES_USE_APPROX"))
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_use_approx;
		}
		else if (key == "IES_LAMBDA_MULTS")
		{
			ies_lam_mults.clear();
			vector<string> tok;
			tokenize(value, tok, ",");
			for (const auto &iscale : tok)
			{
				ies_lam_mults.push_back(convert_cp<double>(iscale));
			}
		}
		else if ((key == "IES_INIT_LAM") || (key == "IES_INITIAL_LAMBDA"))
		{
			convert_ip(value, ies_init_lam);
		}
		else if (key == "IES_USE_APPROX") 
		{
			convert_ip(value, ies_obs_restart_csv);
		}
		else if (key == "IES_SUBSET_SIZE")
		{
			convert_ip(value, ies_subset_size);
		}

		else {

			throw PestParsingError(line, "Invalid key word \"" + key +"\"");
		}
	}
}

ostream& operator<< (ostream &os, const ParameterInfo& val)
{
	for(unordered_map<string, ParameterRec>::const_iterator s=val.parameter_info.begin(),
				e=val.parameter_info.end(); s!=e; ++s) {
			os << (*s).second;
		}
	return os;
}

ostream& operator<< (ostream &os, const ObservationGroupRec& val)
{
	os << "    gtarg = " << val.gtarg << endl;
	os << "    covfile = " << val.covfile << endl;
	return os;
}

ostream& operator<< (ostream &os, const ObservationRec& val)
{
	os << "    weight = " << val.weight << endl;
	os << "    group = " << val.group << endl;
	return os;
}

bool ObservationGroupRec::is_regularization(const string &grp_name)
{
	bool is_reg = false;
	int found = upper_cp(grp_name).find("REGUL");
	if (found == 0) is_reg = true;
	return is_reg;

}

const ObservationRec* ObservationInfo::get_observation_rec_ptr(const string &name) const
{
	const ObservationRec *ret_val = 0;
	unordered_map<string, ObservationRec>::const_iterator p_iter;

	p_iter = observations.find(name);
	if(p_iter != observations.end()) {
		ret_val = &((*p_iter).second);
	}
	return ret_val;
}

const ObservationGroupRec* ObservationInfo::get_group_rec_ptr(const string &name) const
{
	const ObservationGroupRec *ret_val = 0;
	unordered_map<string, ObservationGroupRec>::const_iterator p_iter;

	p_iter = groups.find(name);
	if(p_iter != groups.end()) {
		ret_val = &((*p_iter).second);
	}
	return ret_val;
}


bool ObservationRec::is_regularization() const
{
	bool is_reg = false;
	int found = upper_cp(group).find("REGUL");
	if (found == 0) is_reg=true;
	return is_reg;
}

int ObservationInfo::get_nnz_obs() const
{
	int nnz = 0;
	for (auto &obs : observations)
	{
		if ((!obs.second.is_regularization()) && (obs.second.weight > 0.0))
			nnz++;
	}
	return nnz;
}

int ObservationInfo::get_nnz_obs_and_reg() const
{
	int nnz = 0;
	for (auto &obs : observations)
	{
		if (obs.second.weight > 0.0)
			nnz++;
	}
	return nnz;
}

double ObservationInfo::get_weight(const string &obs_name) const
{
	return observations.find(obs_name)->second.weight;
}

void ObservationInfo::set_weight(const string &obs_name, double &value)
{
	if (observations.find(obs_name) == observations.end())
		throw PestError("ObservationInfo::set_weight() error: observation\
			    " + obs_name + " not found");
	observations[obs_name].weight = value;
}

string ObservationInfo::get_group(const string &obs_name) const
{
	return observations.find(obs_name)->second.group;
}

bool ObservationInfo::is_regularization(const string &obs_name) const
{
	return observations.find(obs_name)->second.is_regularization();
}

Observations ObservationInfo::get_regulatization_obs(const Observations &obs_in)
{
	Observations reg_obs;

	for (Observations::const_iterator b=obs_in.begin(), e=obs_in.end();
			b!=e; ++b)
	{
		if(is_regularization((*b).first)) {
			reg_obs.insert(*b);
		}
	}
	return reg_obs;
}


ostream& operator<< (ostream &os, const ObservationInfo& val)
{
	os << "PEST Observation Information" << endl;
	os << "PEST Observation Groups" << endl;
	for(unordered_map<string, ObservationGroupRec>::const_iterator b=val.groups.begin(),
		e=val.groups.end(); b!=e; ++b) {
			os << "    name = " << (*b).first << endl;
	}
	os << "PEST Observation Information" << endl;
	for(unordered_map<string, ObservationRec>::const_iterator b=val.observations.begin(),
		e=val.observations.end(); b!=e; ++b) {
			os << "  Observation Name = " << (*b).first << endl;
			os << (*b).second;
	}
	return os;
}

ostream& operator<< (ostream &os, const SVDInfo& val)
{
	os << "PEST SVD Information" << endl;
	os << "    maxsing = " << val.maxsing << endl;
	os << "    eigthresh = " << val.eigthresh << endl;
	return os;
}
