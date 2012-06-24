#ifndef MORRISMETHOD_H_
#define MORRISMETHOD_H_

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

class MorrisMethod
{
public:
	MorrisMethod(int _p);
	MatrixXd create_P_star_mat(int k);
	~MorrisMethod(void);
private:
	int k;
	int p;
	double delta;
	static MatrixXd create_B_mat(int k);
	static MatrixXd create_J_mat(int k);
	static MatrixXd create_D_mat(int k);
	static MatrixXd create_P_mat(int k);
	VectorXd create_x_vec(int k);
	static int rand_plus_minus_1(void);
};


#endif /* MORRISMETHOD_H_ */