#ifndef MORRISMETHOD_H_
#define MORRISMETHOD_H_

#include <vector>
#include <string>
#include <Eigen/Dense>
#include "Transformable.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ParamTransformSeq;
class RunManagerAbstract;
class ModelRun;

class MorrisMethod
{
public:
	MorrisMethod(const std::vector<std::string> &par_name_vec, const Parameters &fixed_pars, const Parameters &lower_bnd, 
		const Parameters &upper_bnd, int _p);
    void assemble_runs(RunManagerAbstract &rm, ParamTransformSeq &base_partran_seq);
	void calc_sen(RunManagerAbstract &rm, ModelRun model_run, std::ofstream &fout);
	~MorrisMethod(void);
private:
	int k; // number of parameters
	int p; // number of levels for each parameters
	int r;
	double delta;
    std::vector<std::string> par_name_vec;
	Parameters fixed_pars;
    Parameters lower_bnd;
    Parameters upper_bnd;
    MatrixXd b_star_mat;
	static MatrixXd create_B_mat(int k);
	static MatrixXd create_J_mat(int k);
	static MatrixXd create_D_mat(int k);
	static MatrixXd create_P_mat(int k);
    MatrixXd create_P_star_mat(int k);
	VectorXd create_x_vec(int k);
	Parameters get_ctl_parameters(int row);
	static int rand_plus_minus_1(void);
};


#endif /* MORRISMETHOD_H_ */