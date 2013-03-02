#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cassert>
#include <algorithm>
#include "sobol.h"
#include "Transformable.h"
#include "RunManagerAbstract.h"
#include "ParamTransformSeq.h"
#include "ModelRunPP.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


VectorXd Sobol::gen_uniform_rand_vec(long nsample, double min, double max)
{
  VectorXd v(nsample);
  double dy_2 = (max - min) / 2.0;
  double y_avg = (max + min) / 2.0;
  v.setRandom().array() =  v.setRandom().array() * dy_2 + y_avg ;
  return v;
}

MatrixXd Sobol::gen_N_matrix(const MatrixXd &m1, const MatrixXd &m2, const vector<int> &idx_vec)
{
  MatrixXd n = m2;
  for (int i : idx_vec)
  {
    n.col(i) = m1.col(i);
  }  
  return n;
}

void Sobol::assemble_runs(RunManagerAbstract &rm, ParamTransformSeq &base_partran_seq)
{
	MatrixXd m1;
	MatrixXd m2;
	long npar = par_name_vec.size();
	//generate random matrices
	double par_min;
	double par_max;
	VectorXd v1;
	for (int i=0; i<npar; ++i)
	{
		par_min = lower_bnd[par_name_vec[i]];
		par_max = upper_bnd[par_name_vec[i]];
		v1 = gen_uniform_rand_vec(n_sample, par_min, par_max);
		m1.col(i) = v1;
		v1 = gen_uniform_rand_vec(n_sample, par_min, par_max);
		m2.col(i) = v1;
	}
	// initialize run manager
	vector<string> all_par_names = par_name_vec;
	vector<string> fixed_par_names = fixed_ctl_pars.get_keys();
	all_par_names.insert(all_par_names.end(), fixed_par_names.begin(), fixed_par_names.end());


	//run_manager_ptr->initialize()
	//calculate first order runs a0,a1,....an
	for (int i=0; i<n_sample; ++i)
	{

	}


}

main()
{
  MatrixXd m1;
  MatrixXd m2;
  long npar = 3;
  long nsample = 10;
  m1 = MatrixXd::Zero(nsample,npar);
  for (int i=0; i<npar; ++i)
  {
    //m1.col(0).setConstant(2.0);
    //m1.col(1).setRandom();
    //cout << m1 << endl;
    VectorXd v1 = gen_uniform_rand_vec(nsample, 10, 100);
    m1.col(i) = v1;
  }
  cout << m1 << endl << endl;
  m2 = MatrixXd::Zero(nsample,npar);
  for (int i=0; i<npar; ++i)
  {
    //m1.col(0).setConstant(2.0);
    //m1.col(1).setRandom();
    //cout << m1 << endl;
    VectorXd v1 = gen_uniform_rand_vec(nsample, 10, 100);
    m2.col(i) = v1;
  }
  cout << m2 << endl;
  cout <<  gen_N_matrix(m1, m2, vector<int> idx_vec {1,2})

}
