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

void OutputFileWriter::par_report(std::ostream &os, Parameters const &new_ctl_pars, Parameters const &old_ctl_pars,string par_type)
{	
	double p_old, p_new;
	double fac_change = -9999, rel_change = -9999;
	bool have_fac = false, have_rel = false;
	double max_fac_change = 0;
	double max_rel_change = 0;	
	string max_fac_par = "N/A";
	string max_rel_par = "N/A";
	
	os << "    Parameter Upgrades (" << par_type << " Parameters)" << endl;
	os << "      Parameter     Current       Previous       Factor       Relative" << endl;
	os << "        Name         Value         Value         Change        Change" << endl;
	os << "      ----------  ------------  ------------  ------------  ------------" << endl;
	
	//for (const auto &ipar : new_ctl_pars)
	for (auto &p_name : pest_scenario.get_ctl_ordered_par_names())
	{
		Parameters::const_iterator pi = new_ctl_pars.find(p_name);
		if (pi == new_ctl_pars.end()) continue;
		p_new = new_ctl_pars.get_rec(p_name);
		p_old = old_ctl_pars.get_rec(p_name);		
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

void OutputFileWriter::phi_report(std::ostream &os, map<string, double> const phi_comps, double const dynamic_reg_weight,bool final)
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
	
	double val;
	for (auto &gname : pest_scenario.get_ctl_ordered_obs_group_names())
	{
		os << "    Contribution to phi from observation group ";
		os << setw(17) << setiosflags(ios::right) << "\"" + lower_cp(gname) + "\" : ";
		os << phi_comps.at(gname) << endl;
	}
}


void OutputFileWriter::write_rei(ofstream &fout, int iter_no, const Observations &obs, const Observations &sim, 
	const ObjectiveFunc &obj_func, const Parameters &pars)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << " MODEL OUTPUTS AT END OF OPTIMISATION ITERATION NO. " <<  iter_no << ":-" << endl;
	fout << endl << endl;
	fout <<  setw(21) << " Name" << setw(13) << " Group" <<  setw(21) << " Measured" << setw(21) << " Modelled" << setw(21) << " Residual" << setw(21) << " Weight" << endl;
	//vector<string> obs_name_vec = obs.get_keys();
	vector<string> obs_name_vec = pest_scenario.get_ctl_ordered_obs_names();
	double obs_val, sim_val, residual;
	//for(vector<string>::const_iterator b = obs_name_vec.begin(), 
	//	e = obs_name_vec.end(); b!=e; ++b)
	for (auto &b : obs_name_vec)
	{
		obs_val = obs.get_rec(b);
		sim_val = sim.get_rec(b);
		fout << " " << setw(20) << lower_cp(b)
			<< " " << setw(12) << lower_cp(obj_func.get_obs_info_ptr()->get_observation_rec_ptr(b)->group)
			<< " " << showpoint <<   setw(20) << obs_val 
			<< " " << showpoint <<  setw(20) << sim_val 
			<<  " " << showpoint <<  setw(20) << obs_val - sim_val
			<< " " << showpoint << setw(20)  << sqrt(obj_func.get_obs_info_ptr()->get_observation_rec_ptr(b)->weight) <<endl;
	}
	//process prior information
	const PriorInformation *prior_info_ptr = obj_func.get_prior_info_ptr();
	const PriorInformationRec *pi_rec_ptr;
	PriorInformation::const_iterator ipi;
	//for(PriorInformation::const_iterator b = prior_info_ptr->begin(), 
	//	e = prior_info_ptr->end(); b!=e; ++b)
	obs_name_vec = pest_scenario.get_ctl_ordered_pi_names();
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
	const string *name_ptr;
	Parameters::const_iterator it;
	pair<bool, double> val_pair;
	double scale, offset;

	fout.unsetf(ios::floatfield);
	fout.precision(15);
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

		fout << setw(14) << lower_cp(b) << setw(22) 
		<<  showpoint<< pars.get_rec(b) << " " << setw(20) << showpoint << scale << " " << setw(20) << showpoint << offset << endl;
	}
}


void OutputFileWriter::read_par(ifstream &fin, Parameters &pars)
{
	string line;
	string name;
	double value;
	vector<string> tokens;


	getline(fin, line);

	while (getline(fin, line))
	{
		strip_ip(line);
		tokens.clear();
		tokenize(line, tokens);
		name = tokens[0];
		convert_ip(tokens[1], value);
		pars[name] = value;
	}
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
	const ObjectiveFunc &obj_func, const ParameterGroupInfo &par_grp_info, const DynamicRegularization &regul)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << " NUMERIC PARAMETER SENSITIVITIES FOR OPTIMISATION ITERATION NO. " << setw(3) << iter_no << " ----->" << endl;
	fout << " Parameter name   Group        Current Value           CSS w/reg           CSS w/o reg" << endl;
	const vector<string> &par_list = jac.parameter_list();
	//const vector<string> &par_list = pest_scenario.get_ctl_ordered_par_names();
	const vector<string> &obs_list = jac.obs_and_reg_list();
	Parameters pars = jac.get_base_numeric_parameters();
	VectorXd par_vec = pars.get_data_eigen_vec(par_list);
	MatrixXd par_mat_tmp = par_vec.asDiagonal();
	Eigen::SparseMatrix<double> par_mat = par_mat_tmp.sparseView();
	QSqrtMatrix Q_sqrt(obj_func.get_obs_info_ptr(), obj_func.get_prior_info_ptr());
	Eigen::SparseMatrix<double> q_sqrt_reg = Q_sqrt.get_sparse_matrix(obs_list, regul);
	Eigen::SparseMatrix<double> dss_mat_reg = q_sqrt_reg * jac.get_matrix(obs_list, par_list) * par_mat;
	Eigen::SparseMatrix<double> q_sqrt_no_reg = Q_sqrt.get_sparse_matrix(obs_list, DynamicRegularization::get_zero_reg_instance());
	Eigen::SparseMatrix<double> dss_mat_no_reg = q_sqrt_no_reg * jac.get_matrix(obs_list, par_list) * par_mat;

	int n_par = par_list.size();
	int n_nonzero_weights_reg = q_sqrt_reg.nonZeros();
	int n_nonzero_weights_no_reg = q_sqrt_no_reg.nonZeros();
	const vector<string> ordered_par_names = pest_scenario.get_ctl_ordered_par_names();
	//drop any names that aren't in par_list

	vector<string>::const_iterator is;
	int i;
	//for (int i = 0; i < n_par; ++i)
	for (auto &pname : ordered_par_names)
	{
		is = find(par_list.begin(), par_list.end(), pname);
		if (is == par_list.end()) continue;		
		i = is - par_list.begin();
		fout << "   " << setw(15) << lower_cp(par_list[i])
			<< " " << setw(12) << lower_cp(par_grp_info.get_group_name(par_list[i]))
			<< " " << showpoint << setw(20) << pars.get_rec(par_list[i]);
		if (n_nonzero_weights_reg > 0)
		{
			fout << " " << showpoint << setw(20) << dss_mat_reg.col(i).norm() / double(n_nonzero_weights_reg);
		}
		else
		{
			fout << " " << showpoint << setw(20) << "NA";
		}
		if (n_nonzero_weights_no_reg > 0)
		{
			fout << " " << showpoint << setw(20) << dss_mat_no_reg.col(i).norm() / double(n_nonzero_weights_no_reg);
		}
		else
		{
			fout << " " << showpoint << setw(20) << "NA";
		}
		fout << endl;
	}
	fout << endl << endl;
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
