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

OutputFileWriter::OutputFileWriter(FileManager &_file_manager, const string &_case_name, bool restart_flag, bool _save_rei, int _eigenwrite)
	: file_manager(_file_manager), case_name(_case_name), save_rei(_save_rei), eigenwrite(_eigenwrite)
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

void OutputFileWriter::write_rei(ofstream &fout, int iter_no, const Observations &obs, const Observations &sim, 
	const ObjectiveFunc &obj_func, const Parameters &pars)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << " MODEL OUTPUTS AT END OF OPTIMISATION ITERATION NO. " <<  iter_no << ":-" << endl;
	fout << endl << endl;
	fout <<  setw(21) << " Name" << setw(13) << " Group" <<  setw(21) << " Measured" << setw(21) << " Modelled" << setw(21) << " Residual" << setw(21) << " Weight" << endl;
	vector<string> obs_name_vec = obs.get_keys();
	double obs_val, sim_val, residual;
	for(vector<string>::const_iterator b = obs_name_vec.begin(), 
		e = obs_name_vec.end(); b!=e; ++b)
	{
		obs_val = obs.get_rec(*b);
		sim_val = sim.get_rec(*b);
		fout << " " << setw(20) << *b
			<< " " << setw(12) << obj_func.get_obs_info_ptr()->get_observation_rec_ptr(*b)->group
			<< " " << showpoint <<   setw(20) << obs_val 
			<< " " << showpoint <<  setw(20) << sim_val 
			<<  " " << showpoint <<  setw(20) << obs_val - sim_val
			<< " " << showpoint << setw(20)  << sqrt(obj_func.get_obs_info_ptr()->get_observation_rec_ptr(*b)->weight) <<endl;
	}
	//process prior information
	const PriorInformation *prior_info_ptr = obj_func.get_prior_info_ptr();
	const PriorInformationRec *pi_rec_ptr;
	for(PriorInformation::const_iterator b = prior_info_ptr->begin(), 
		e = prior_info_ptr->end(); b!=e; ++b)
	{
		pi_rec_ptr = &(*b).second;
		obs_val = pi_rec_ptr->get_obs_value();
		residual = pi_rec_ptr->calc_residual(pars);
		sim_val = obs_val + residual;
		fout << " " << setw(20) << (*b).first
		<< " " << setw(12) << pi_rec_ptr->get_group()
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
	for(Parameters::const_iterator b=pars.begin(), e=pars.end();
		b!=e; ++b)
	{
		name_ptr = &(*b).first;
		val_pair = offset_tran.get_value(*name_ptr);
		if (val_pair.first == true)
			offset = val_pair.second;
		else
			offset = 0.0;
		val_pair = scale_tran.get_value(*name_ptr);
		if (val_pair.first == true)
			scale = val_pair.second;
		else
			scale = 1.0;

		fout << setw(14) << *name_ptr << setw(22) 
		<<  showpoint<< (*b).second << " " << setw(20) << showpoint << scale << " " << setw(20) << showpoint << offset << endl;
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

	size_t size_par_vec = pars.size();
	for (int i = 0; i < n_par; ++i)
	{
		fout << "   " << setw(15) << par_list[i]
			<< " " << setw(12) << par_grp_info.get_group_name(par_list[i])
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
