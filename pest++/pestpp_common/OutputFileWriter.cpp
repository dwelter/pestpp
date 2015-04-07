/*  
	� Copyright 2012, David Welter
	
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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include "OutputFileWriter.h"
#include "Transformable.h"
#include "Pest.h"
#include "ObjectiveFunc.h"
#include "PriorInformation.h"
#include "Jacobian.h"
#include "QSqrtMatrix.h"
#include "ModelRunPP.h"
#include "eigen_tools.h"
#include "utilities.h"

using namespace std;
using namespace Eigen;
using namespace::pest_utils;


OutputFileWriter::OutputFileWriter(FileManager &_file_manager, Pest &_pest_scenario, bool restart_flag, bool _save_rei, int _eigenwrite)
	: file_manager(_file_manager), pest_scenario(_pest_scenario),case_name(_file_manager.get_base_filename()), save_rei(_save_rei), eigenwrite(_eigenwrite)
{
	if (restart_flag)
	{
		ofstream &fout_sen = file_manager.open_ofile_ext("sen", ofstream::app);
		write_restart_header(fout_sen);
		ofstream &fout_svd = file_manager.open_ofile_ext("svd", ofstream::app);
		write_restart_header(fout_svd);
	}
	else
	{
		ofstream &fout_sen = file_manager.open_ofile_ext("sen");
		write_sen_header(fout_sen, case_name);
		file_manager.open_ofile_ext("svd");
	}
	if (pest_scenario.get_pestpp_options().get_iter_summary_flag())
	{
		prepare_iteration_summary_files(restart_flag);
	}
}

void OutputFileWriter::iteration_report(std::ostream &os, int iter, int nruns, string iteration_type, string svd_type, string mat_inv)
{
	os << "OPTIMISATION ITERATION NUMBER: " << iter << endl << endl;
	os << "  Iteration type: " << iteration_type << endl;
	os << "    SVD Package: " << svd_type << endl;
	os << "    Matrix Inversion: " << mat_inv << endl;
	os << "    Model calls so far : " << nruns << endl;
	os << endl;
}

void OutputFileWriter::prepare_iteration_summary_files(bool restart_flag)
{
	if (restart_flag)
	{
		file_manager.open_ofile_ext("ipar", ofstream::app);
		file_manager.open_ofile_ext("iobj", ofstream::app);
		file_manager.open_ofile_ext("isen", ofstream::app);
	}
	else
	{
		file_manager.open_ofile_ext("ipar");
		file_manager.open_ofile_ext("iobj");
		file_manager.open_ofile_ext("isen");
		ofstream &os_ipar = file_manager.get_ofstream("ipar");
		ofstream &os_isen = file_manager.get_ofstream("isen");
		os_ipar << "iteration";
		os_isen << "iteration";
		for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
		{
			os_ipar << ',' << lower_cp(par_name);
			os_isen << ',' << lower_cp(par_name);
		}
		os_ipar << endl;
		os_isen << endl;
		ofstream &os_iobj = file_manager.get_ofstream("iobj");
		os_iobj << "iteration,model_runs_completed,total_phi,measurement_phi,regularization_phi";
		for (auto &obs_grp : pest_scenario.get_ctl_ordered_obs_group_names())
		{
			os_iobj << ',' << lower_cp(obs_grp);
		}
		os_iobj << endl;
	}
}

void OutputFileWriter::write_sen_iter(int iter, map<string, double> &ctl_par_sens)
{
	ofstream &os = file_manager.get_ofstream("isen");
	os << iter;
	map<string, double>::iterator is;
	for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
	{
		is = ctl_par_sens.find(par_name);
		if (is == ctl_par_sens.end())
			os << ',' << -999;
		else
			os << ',' << ctl_par_sens.at(par_name);
	}
	os << endl;
}

void OutputFileWriter::write_par_iter(int iter, Parameters const &ctl_pars)
{
	ofstream &os = file_manager.get_ofstream("ipar");
	os << iter;
	for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
	{
		os << ',' << ctl_pars.get_rec(par_name);
	}
	os << endl;
}

void OutputFileWriter::write_obj_iter(int iter, int nruns, map<string, double> const &phi_report)
{
	ofstream &os = file_manager.get_ofstream("iobj");
	os << iter << ',' << nruns;
	os << ',' << phi_report.at("TOTAL");
	os << ',' << phi_report.at("MEAS");
	map<string, double>::const_iterator iregul = phi_report.find("REGUL");
	double val = 0.0;	
	if (iregul != phi_report.end())
	{
		val = phi_report.at("REGUL");
	}
	os << ',' << val;
	for (auto &obs_grp : pest_scenario.get_ctl_ordered_obs_group_names())
	{
		os << ',' << phi_report.at(obs_grp);
	}
	os << endl;
}

void OutputFileWriter::scenario_report(std::ostream &os)
{
	string mode = "estimation";
	if (pest_scenario.get_regul_scheme_ptr()->get_use_dynamic_reg())
	{
		mode = "regularization (with a \"z\")";
	}
	os << endl << "PEST++ run mode:- " << endl << "   " << mode << endl << endl;

	os << "Case dimensions:- " << endl;
	os << setw(0) << "    Number of parameters = " << pest_scenario.get_ctl_ordered_par_names().size() << endl;
	os << setw(0) << "    Number of adjustable parameters = " << pest_scenario.get_n_adj_par() << endl;
	os << setw(0) << "    Number of observations = " << pest_scenario.get_ctl_ordered_obs_names().size() << endl;
	os << setw(0) << "    Number of prior estimates = "  << pest_scenario.get_ctl_ordered_pi_names().size() << endl << endl;

	os << "Model command line(s):- " << endl;
	for (auto &cmd : pest_scenario.get_comline_vec())
	{
		os << "    " << cmd << endl;
	}
	os << endl;

	os << "Model interface files:-" << endl;
	os << "    template files:" << endl;
	for (auto &f : pest_scenario.get_tplfile_vec())
	{
		os << "      " << f << endl;
	}
	os << "    model input files:" << endl;
	for (auto &f : pest_scenario.get_inpfile_vec())
	{
		os << "      " << f << endl;
	}
	os << endl <<  "    instruction files:" << endl;
	for (auto &f : pest_scenario.get_insfile_vec())
	{
		os << "      " << f << endl;
	}
	os << "    model output files:" << endl;
	for (auto &f : pest_scenario.get_outfile_vec())
	{
		os << "      " << f << endl;
	}
	os << endl << endl;

	os << pest_scenario.get_control_info() << endl;
	os << pest_scenario.get_pestpp_options() << endl;

	const ParameterGroupRec *grp_rec;
	map<int, string> trans_type;
	trans_type[0] = "none";
	trans_type[1] = "fixed";
	trans_type[2] = "tied";
	trans_type[3] = "log";
	os << "Parameter group information" << endl;
	os << left << setw(15) << "NAME" << right << setw(15) << "INCREMENT TYPE" << setw(25) << "DERIVATIVE INCREMENT";
	os << setw(25) << "INCREMENT LOWER BOUND" << setw(15) << "FORCE CENTRAL" << setw(25) << "INCREMENT MULTIPLIER" << endl;
	for (auto &grp_name : pest_scenario.get_ctl_ordered_par_group_names())
	{
		grp_rec = pest_scenario.get_base_group_info().get_group_by_groupname(grp_name);
		os << left << setw(15) << lower_cp(grp_rec->name) << right << setw(15) << grp_rec->inctyp << setw(25) << grp_rec->derinc;
		os << setw(25) << grp_rec->derinclb << setw(15) << grp_rec->forcen << setw(25) << grp_rec->derincmul << endl;
	}
	os << endl << "Parameter information" << endl;
	os << left << setw(15) << "NAME" << setw(10) << "TRANSFORMATION" << right << setw(20) << "CHANGE LIMIT" << setw(15) << "INITIAL VALUE";
	os << setw(15) << "LOWER BOUND";	
	os << setw(15) << "UPPER BOUND" << setw(15) << "GROUP";
	
	os << setw(15) << "SCALE" << setw(15) << "OFFSET" << setw(20) << "DERIVATIVE COMMAND" << endl;
	const ParameterRec* par_rec;
	for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
	{
		par_rec = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(par_name);		
		os << left <<setw(15) << lower_cp(par_name);
		os << setw(10) << trans_type[static_cast<int>(par_rec->tranform_type)];
		os << right << setw(20) << par_rec->chglim;
		os << setw(15) << par_rec->init_value;
		os << setw(15) << par_rec->lbnd;
		os << setw(15) << par_rec->ubnd;
		os << setw(15) << lower_cp(par_rec->group);
		os << setw(15) << par_rec->scale;
		os << setw(15) << par_rec->offset;
		os << setw(20) << par_rec->dercom << endl;
	}
	os << endl << "Observation information" << endl;
	os << left << setw(25) << "NAME" << right << setw(20) << "VALUE" << setw(20) << "GROUP" << setw(20) << "WEIGHT" << endl;
	const ObservationRec* obs_rec;
	const Observations &obs = pest_scenario.get_ctl_observations();
	for (auto &obs_name : pest_scenario.get_ctl_ordered_obs_names())
	{
		obs_rec = pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(obs_name);		
		os << left << setw(25) << lower_cp(obs_name);
		os << right << setw(20) << obs.get_rec(obs_name);
		os << setw(20) << lower_cp(obs_rec->group);
		os << setw(20) << obs_rec->weight << endl;
	}
	const PriorInformation &pi = pest_scenario.get_prior_info();
	os << endl << "Prior information" << endl;
	if (pi.size() == 0)
	{
		os << endl << "   no prior information provided" << endl;
	}
	for (auto &pi_name : pi)
	{
		os << pi_name.first << "  " << pi_name.second;
	}
	os << endl << pest_scenario.get_svd_info() << endl;
	os << endl;
}


void OutputFileWriter::par_report(std::ostream &os, Parameters const &new_ctl_pars)
{
	double val;
	os << endl;
	os << "     Parameter            " << endl;
	os << "        Name         Value" << endl;
	os << "    ------------  ------------" << endl;
	for (auto &p_name : pest_scenario.get_ctl_ordered_par_names())
	{
		val = new_ctl_pars.get_rec(p_name);		
		os << left;
		os << "    " << setw(12) << lower_cp(p_name);
		os << right;
		os << "  " << setw(12) << val << endl;
	}
	os << endl;
}


void OutputFileWriter::par_report(std::ostream &os, int const iter, Parameters const &new_pars, Parameters const &old_pars,
	string par_type)
{	
	double p_old, p_new;
	double fac_change = -9999, rel_change = -9999;
	bool have_fac = false, have_rel = false;
	double max_fac_change = 0;
	double max_rel_change = 0;	
	string max_fac_par = "N/A";
	string max_rel_par = "N/A";
	
	os << "    Iteration "<<iter<<" Parameter Upgrades (" << par_type << " Parameters) " << endl;
	os << "      Parameter     Current       Previous       Factor       Relative" << endl;
	os << "        Name         Value         Value         Change        Change" << endl;
	os << "      ----------  ------------  ------------  ------------  ------------" << endl;
	vector<string> par_names;
	if (lower_cp(par_type) == "control file")
		par_names = pest_scenario.get_ctl_ordered_par_names();
	else
	{
		par_names = new_pars.get_keys();
		sort(par_names.begin(), par_names.end());
	}
	//for (const auto &ipar : new_ctl_pars)	
	for (auto &p_name : par_names)
	{
		Parameters::const_iterator pi = new_pars.find(p_name);
		if (pi == new_pars.end()) continue;
		p_new = new_pars.get_rec(p_name);
		p_old = old_pars.get_rec(p_name);		
		param_change_stats(p_old, p_new, have_fac, fac_change, have_rel, rel_change);
		if (have_fac && fac_change >= max_fac_change)
		{
			max_fac_change = fac_change;
			max_fac_par = p_name;
		}
		if (have_rel && abs(rel_change) >= abs(max_rel_change))
		{
			max_rel_change = rel_change;
			max_rel_par = p_name;
		}
		os << right;
		os << "    " << setw(12) << lower_cp(p_name);
		os << right;
		os << "  " << setw(12) << p_new;
		os << "  " << setw(12) << p_old;
		if (have_fac)
			os << "  " << setw(12) << fac_change;
		else
			os << "  " << setw(12) << "N/A";
		if (have_rel)
			os << "  " << setw(12) << rel_change;
		else
			os << "  " << setw(12) << "N/A";
		os << endl;
	}
	os << "       Maximum changes in \"" << par_type << "\" parameters:" << endl;
	os << "         Maximum relative change = " << max_rel_change << "   [" << lower_cp(max_rel_par) << "]" << endl;
	os << "         Maximum factor change = " << max_fac_change << "   [" << lower_cp(max_fac_par) << "]" << endl;	
	os << endl;	
	if ((lower_cp(par_type) == "control file") && (pest_scenario.get_pestpp_options().get_iter_summary_flag()))
	{
		write_par_iter(iter, new_pars);
	}
}


void OutputFileWriter::param_change_stats(double p_old, double p_new, bool &have_fac, 
	double &fac_change, bool &have_rel, double &rel_change)
{
	have_rel = have_fac = true;
	double a = max(abs(p_new), abs(p_old));
	double b = min(abs(p_new), abs(p_old));
	// compute relative change
	if (p_old == 0) {
		have_rel = false;
		rel_change = -9999;
	}
	else
	{
		rel_change = (p_old - p_new) / p_old;
	}
	//compute factor change
	if (p_old == 0.0 || p_new == 0.0) {
		have_fac = false;
		fac_change = -9999;
	}
	else {
		fac_change = a / b;
	}
}


void OutputFileWriter::phi_report(std::ostream &os, int const iter, int const nruns, map<string, double> const phi_comps, double const dynamic_reg_weight,bool final)
{
	map<string, double>::const_iterator it = phi_comps.find("REGUL");
	if ((!dynamic_reg_weight) || (it == phi_comps.end()))
	{
		if (final)
		{
			os << "    Final phi                                           Total : " << phi_comps.at("TOTAL") << endl;
		}
		else
		{
			os << "    Starting phi for this iteration                     Total : " << phi_comps.at("TOTAL") << endl;
		}		
	}
	else
	{
		if (final)
		{
			os << "    Final regularization weight factor                        : " << dynamic_reg_weight << endl;
			os << "    Final phi                                           Total : " << phi_comps.at("TOTAL") << endl;
			os << "    Final measurement phi for this iteration            Total : " << phi_comps.at("MEAS") << endl;
			os << "    Final regularization phi for this iteration         Total : " << phi_comps.at("REGUL") << endl;
		}
		else
		{
			os << "    Current regularization weight factor                      : " << dynamic_reg_weight << endl;
			os << "    Starting phi for this iteration                     Total : " << phi_comps.at("TOTAL") << endl;
			os << "    Starting measurement phi for this iteration         Total : " << phi_comps.at("MEAS") << endl;
			os << "    Starting regularization phi for this iteration      Total : " << phi_comps.at("REGUL") << endl;
		}
	}
	
	for (auto &gname : pest_scenario.get_ctl_ordered_obs_group_names())
	{
		os << "    Contribution to phi from observation group ";
		os << setw(17) << setiosflags(ios::right) << "\"" + lower_cp(gname) + "\" : ";
		os << phi_comps.at(gname) << endl;
	}
	/*if (pest_scenario.get_pestpp_options().get_iter_summary_flag())
	{
		write_obj_iter(iter, nruns,phi_comps);
	}*/
}


void OutputFileWriter::obs_report(ostream &os, const Observations &obs, const Observations &sim, const ObjectiveFunc &obj_func)
{
	
	os << setw(21) << " Name" << setw(13) << " Group" << setw(21) << " Measured" << setw(21) << " Modelled" << setw(21) << " Residual" << setw(21) << " Weight" << endl;
	//vector<string> obs_name_vec = obs.get_keys();
	vector<string> obs_name_vec = pest_scenario.get_ctl_ordered_obs_names();
	double obs_val, sim_val;
	//for(vector<string>::const_iterator b = obs_name_vec.begin(), 
	//	e = obs_name_vec.end(); b!=e; ++b)
	for (auto &b : obs_name_vec)
	{
		obs_val = obs.get_rec(b);
		sim_val = sim.get_rec(b);
		os << " " << setw(20) << lower_cp(b)
			<< " " << setw(12) << lower_cp(obj_func.get_obs_info_ptr()->get_observation_rec_ptr(b)->group)
			<< " " << showpoint << setw(20) << obs_val
			<< " " << showpoint << setw(20) << sim_val
			<< " " << showpoint << setw(20) << obs_val - sim_val
			<< " " << showpoint << setw(20) << obj_func.get_obs_info_ptr()->get_observation_rec_ptr(b)->weight << endl;
	}

}

void OutputFileWriter::write_rei(ofstream &fout, int iter_no, const Observations &obs, const Observations &sim, 
	const ObjectiveFunc &obj_func, const Parameters &pars)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << " MODEL OUTPUTS AT END OF OPTIMISATION ITERATION NO. " << iter_no << ":-" << endl;
	fout << endl << endl;
	obs_report(fout, obs, sim, obj_func);
	//process prior information
	const PriorInformation *prior_info_ptr = obj_func.get_prior_info_ptr();
	const PriorInformationRec *pi_rec_ptr;
	PriorInformation::const_iterator ipi;
	//for(PriorInformation::const_iterator b = prior_info_ptr->begin(), 
	//	e = prior_info_ptr->end(); b!=e; ++b)
	vector<string> obs_name_vec = pest_scenario.get_ctl_ordered_pi_names();
	double obs_val, residual, sim_val;
	for (auto &b : obs_name_vec)
	{
		ipi = prior_info_ptr->find(b);
		pi_rec_ptr = &(*ipi).second;
		obs_val = pi_rec_ptr->get_obs_value();
		residual = pi_rec_ptr->calc_residual(pars);
		sim_val = obs_val + residual;
		fout << " " << setw(20) << lower_cp(b)
		<< " " << setw(12) << lower_cp(pi_rec_ptr->get_group())
		<< " " << showpoint <<   setw(20) << obs_val
		<< " " << showpoint <<  setw(20) << sim_val
		<<  " " << showpoint <<  setw(20) << residual
		<< " " << showpoint << setw(20)  << sqrt(pi_rec_ptr->get_weight()) << endl;	
	}
}


void OutputFileWriter::write_par(ofstream &fout, const Parameters &pars, const TranOffset &offset_tran, const TranScale &scale_tran)
{
	Parameters::const_iterator it;
	pair<bool, double> val_pair;
	double scale, offset;

	fout.unsetf(ios::floatfield);
	fout.precision(numeric_limits<double>::digits10 + 1);
	fout << "single point" << endl;
	//for(Parameters::const_iterator b=pars.begin(), e=pars.end();
	//	b!=e; ++b)
	vector<string> par_name_vec = pest_scenario.get_ctl_ordered_par_names();
	for (auto &b : par_name_vec)
	{
		//name_ptr = &(*b).first;		
		val_pair = offset_tran.get_value(b);
		if (val_pair.first == true)
			offset = val_pair.second;
		else
			offset = 0.0;
		val_pair = scale_tran.get_value(b);
		if (val_pair.first == true)
			scale = val_pair.second;
		else
			scale = 1.0;

		fout << setw(14) << lower_cp(b) << setw(22) << " "
		<<  showpoint<< pars.get_rec(b) << " " << setw(20) << showpoint << scale << " " << setw(20) << showpoint << offset << endl;
	}
	fout.flush();
}

void OutputFileWriter::write_sen_header(std::ostream &fout, const string &case_name)
{
	fout << "                    PARAMETER SENSITIVITIES: CASE " << case_name << endl;
	fout << endl << endl;
}


void OutputFileWriter::write_restart_header(std::ostream &fout)
{
	fout << endl << endl;
	fout << "Restarting PEST++ .... " << endl;
	fout << endl << endl;
}


void OutputFileWriter::append_sen(std::ostream &fout, int iter_no, const Jacobian &jac,
	const ObjectiveFunc &obj_func, const ParameterGroupInfo &par_grp_info, const DynamicRegularization &regul,
	bool is_super)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << "NUMERIC PARAMETER SENSITIVITIES FOR OPTIMISATION ITERATION NO. " << setw(3) << iter_no << " ----->" << endl;
	fout << " Parameter name   Group        Current Value           CSS w/reg           CSS w/o reg" << endl;
	const vector<string> &par_list = jac.parameter_list();
	//const vector<string> &par_list = pest_scenario.get_ctl_ordered_par_names();
	const vector<string> &obs_list = jac.obs_and_reg_list();
	Parameters pars = jac.get_base_numeric_parameters();
	if (pars.size() == 0)
	{
		fout << "parameter values are not avaialble to compute CSS" << endl;
		fout << endl << endl;
	}
	else
	{
		VectorXd par_vec = pars.get_data_eigen_vec(par_list);
		MatrixXd par_mat_tmp = par_vec.asDiagonal();
		Eigen::SparseMatrix<double> par_mat = par_mat_tmp.sparseView();
		QSqrtMatrix Q_sqrt(obj_func.get_obs_info_ptr(), obj_func.get_prior_info_ptr());
		Eigen::SparseMatrix<double> q_sqrt_reg = Q_sqrt.get_sparse_matrix(obs_list, regul);
		Eigen::SparseMatrix<double> dss_mat_reg = q_sqrt_reg * jac.get_matrix(obs_list, par_list) * par_mat;
		Eigen::SparseMatrix<double> q_sqrt_no_reg = Q_sqrt.get_sparse_matrix(obs_list, DynamicRegularization::get_zero_reg_instance());
		Eigen::SparseMatrix<double> dss_mat_no_reg = q_sqrt_no_reg * jac.get_matrix(obs_list, par_list) * par_mat;

		//cout << q_sqrt_reg << endl << endl;
		//cout << q_sqrt_no_reg << endl << endl;

		int n_par = par_list.size();
		int n_nonzero_weights_reg = q_sqrt_reg.nonZeros();
		int n_nonzero_weights_no_reg = q_sqrt_no_reg.nonZeros();
		vector<string> par_names;
		if (is_super)
		{
			par_names = par_list;
			sort(par_names.begin(), par_names.end());
		}
		else
		{
			par_names = pest_scenario.get_ctl_ordered_par_names();
		}
		//for isen file
		map<string, double> par_sens;
		double val;
		//drop any names that aren't in par_list
		vector<string>::const_iterator is;
		int i;
		//for (int i = 0; i < n_par; ++i)
		for (auto &pname : par_names)
		{
			is = find(par_list.begin(), par_list.end(), pname);
			if (is == par_list.end()) continue;
			i = is - par_list.begin();
			fout << "   " << setw(15) << lower_cp(par_list[i])
				<< " " << setw(12) << lower_cp(par_grp_info.get_group_name(par_list[i]))
				<< " " << showpoint << setw(20) << pars.get_rec(par_list[i]);
			if (n_nonzero_weights_reg > 0)
			{
				fout << " " << showpoint << setw(20) << dss_mat_reg.col(i).norm() / pow(n_nonzero_weights_reg, 2.0);
			}
			else
			{
				fout << " " << showpoint << setw(20) << "NA";
			}
			if (n_nonzero_weights_no_reg > 0)
			{
				val = dss_mat_no_reg.col(i).norm() / pow(n_nonzero_weights_no_reg, 2.0);
				par_sens[pname] = val;
				fout << " " << showpoint << setw(20) << val;

			}
			else
			{
				fout << " " << showpoint << setw(20) << "NA";
				par_sens[pname] = -999;
			}
			fout << endl;
		}
		fout << endl << endl;
		if (!is_super)
		{
			write_sen_iter(iter_no, par_sens);
		}
	}
}


void OutputFileWriter::set_svd_output_opt(int _eigenwrite)
{
	eigenwrite = _eigenwrite;
}


void OutputFileWriter::write_svd(VectorXd &Sigma, Eigen::SparseMatrix<double> &Vt, double lambda, const Parameters &freeze_numeric_pars, VectorXd &Sigma_trunc)
{
		ofstream &fout_svd = file_manager.get_ofstream("svd");
		fout_svd<< "CURRENT VALUE OF MARQUARDT LAMBDA = " << lambda << " --------->" << endl << endl;
		fout_svd << "FROZEN PARAMETERS-" << endl;
		fout_svd << freeze_numeric_pars << endl << endl;
		fout_svd << "SINGULAR VALUES IN SOLUTION:-" << endl;
		print(Sigma, fout_svd, 7);
		fout_svd << "TRUNCATED SINGULAR VALUES:-" << endl;
		print(Sigma_trunc, fout_svd, 7);
		fout_svd << endl << endl;
		if (eigenwrite > 0)
		{
			fout_svd << "MATRIX OF EIGENVECTORS IN SOLUTION:-" << endl;
			print(Vt.toDense(), fout_svd, 7);
			fout_svd << endl;
		}
		fout_svd << "Number of singular values used in solution = " << Sigma.size() << endl << endl << endl;
}


void OutputFileWriter::write_svd_iteration(int iteration_no)
{
		ofstream &fout_svd = file_manager.get_ofstream("svd");
	// write head for SVD file
		fout_svd << "------------------------------------------------------------------------------" << endl;
		fout_svd << "OPTIMISATION ITERATION NO.        : " << iteration_no << endl << endl;
}
