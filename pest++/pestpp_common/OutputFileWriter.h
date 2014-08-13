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
#ifndef OUTPUTFILEWRITER_H
#define OUTPUTFILEWRITER_H

#include <string>
#include <iostream>
#include <fstream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include "FileManager.h"
#include "Pest.h"

class Observations;
class ObjectiveFunc;
class Parameters;
class TranOffset;
class TranScale;
class Jacobian;
class QSqrtMatrix;
class ParameterGroupInfo;
class DynamicRegularization;

class OutputFileWriter
{
public:
	OutputFileWriter(FileManager &_file_manager, Pest &_pest_scenario, bool restart_flag = false, bool _save_rei = true, int _eigenwrite = 0);
	void write_rei(std::ofstream &fout, int iter_no, const Observations &obs,
		const Observations &sim, const ObjectiveFunc &obj_func, const Parameters &pars);
	void write_par(std::ofstream &fout, const Parameters &pars, const TranOffset &offset_tran, const TranScale &scale_tran);
	void read_par(std::ifstream &fin, Parameters &pars);
	void write_restart_header(std::ostream &fout);
	void write_sen_header(std::ostream &fout, const std::string &case_name);
	void set_svd_output_opt(int _eigenwrite);
	void append_sen(std::ostream &fout, int iter_no, const Jacobian &jac, const ObjectiveFunc &obj_func, const ParameterGroupInfo &par_grp_info, const DynamicRegularization &regul);
	void write_svd(Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double> &Vt, double lambda, const Parameters &freeze_numeric_pars, Eigen::VectorXd &Sigma_trunc);
	void write_svd_iteration(int iteration_no);
private:
	FileManager &file_manager;
	Pest &pest_scenario;
	std::string case_name;
	int eigenwrite;
	bool save_rei;
};
#endif /* OUTPUTFILEWRITER_H */
