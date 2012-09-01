#ifndef MORRISMETHOD_H_
#define MORRISMETHOD_H_

#include <vector>
#include <string>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
class Parameters;

class MorrisMethod
{
public:
	MorrisMethod(int _p);
	MatrixXd create_P_star_mat(int k);
	~MorrisMethod(void);
private:
	int k; // number of parameters
	int p; // number of levels for each parameters
	double delta;
	static MatrixXd create_B_mat(int k);
	static MatrixXd create_J_mat(int k);
	static MatrixXd create_D_mat(int k);
	static MatrixXd create_P_mat(int k);
	VectorXd create_x_vec(int k);
	Parameters get_ctl_parameters(const MatrixXd &b_mat, int row,
		const std::vector<std::string> &par_name_vec, const Parameters &lower_bnd, 
		const Parameters &upper_bnd);
	static int rand_plus_minus_1(void);
};


#endif /* MORRISMETHOD_H_ */