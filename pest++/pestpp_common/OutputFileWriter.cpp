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

OutputFileWriter::OutputFileWriter(FileManager &_file_manager, const string &_case_name, bool _save_rei, bool _save_svd)
	: file_manager(_file_manager), case_name(_case_name), save_rei(_save_rei), save_svd(_save_svd)
{

	ofstream &fout_sen = file_manager.open_ofile_ext("sen");
	write_sen_header(fout_sen, case_name);
	if (save_svd)
	{
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

void OutputFileWriter::append_sen(std::ostream &fout, int iter_no, const Jacobian &jac, const ObjectiveFunc &obj_func, const ParameterGroupInfo &par_grp_info)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << " OPTIMISATION ITERATION NO. " << setw(3) << iter_no << " ----->" << endl;
	fout << " Parameter name    Group        Current value        Sensitivity" << endl;
	const vector<string> &par_list = jac.parameter_list();
	const vector<string> &obs_list = jac.obs_and_reg_list();
	Parameters pars = jac.get_base_numeric_parameters();
	QSqrtMatrix Q_sqrt(obj_func.get_obs_info_ptr(), obj_func.get_prior_info_ptr(), 1.0);
	Eigen::SparseMatrix<double> q_sqrt = Q_sqrt.get_sparse_matrix(obs_list);
	Eigen::SparseMatrix<double> w_sen_mat = q_sqrt * jac.get_matrix(obs_list, par_list);
	int n_par = par_list.size();
	int n_nonzero_weights = q_sqrt.nonZeros();
	double nonzero_weights_fac = 0.0;
	if (n_nonzero_weights > 0) nonzero_weights_fac = 1.0 / double(n_nonzero_weights);
	size_t size_par_vec = pars.size();
	for (int i=0; i<n_par; ++i) 
	{
		 if (size_par_vec==0)
		 {
			fout << "   " << setw(15) << par_list[i]
			<< " " << setw(12) << par_grp_info.get_group_name(par_list[i])
			<< " " << showpoint <<   setw(20) << "NA"
			<< " " << showpoint <<   setw(20) << w_sen_mat.col(i).norm() * nonzero_weights_fac << endl;
		 }
		 else
		 {
			fout << "   " << setw(15) << par_list[i]
			<< " " << setw(12) << par_grp_info.get_group_name(par_list[i])
			<< " " << showpoint <<   setw(20) << pars.get_rec(par_list[i])
			<< " " << showpoint <<   setw(20) << w_sen_mat.col(i).norm() * nonzero_weights_fac << endl;
		 }
	}
	fout << endl << endl;
}

void OutputFileWriter::save_svd_output(bool _save_svd)
{
	save_svd = _save_svd;
	if (save_svd)
	{
		file_manager.open_ofile_ext("svd");
	}
}


void OutputFileWriter::write_svd(VectorXd &Sigma, Eigen::SparseMatrix<double> &Vt, double lambda, const Parameters &freeze_numeric_pars, VectorXd &Sigma_trunc)
{
	if (OutputFileWriter::save_svd)
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
		fout_svd << "MATRIX OF EIGENVECTORS IN SOLUTION:-" << endl;
		print(Vt.toDense(), fout_svd, 7);
		fout_svd << endl;
		fout_svd << "Number of singular values used in solution = " << Sigma.size() << endl << endl << endl;
	}
}

void OutputFileWriter::write_svd_iteration(int iteration_no)
{
	if (OutputFileWriter::save_svd)
	{
		ofstream &fout_svd = file_manager.get_ofstream("svd");
	// write head for SVD file
		fout_svd << "------------------------------------------------------------------------------" << endl;
		fout_svd << "OPTIMISATION ITERATION NO.        : " << iteration_no << endl << endl;
	}
}
