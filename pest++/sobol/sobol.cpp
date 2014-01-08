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
using namespace Eigen;


Sobol::Sobol(RunManagerAbstract *_run_manager_ptr, const vector<string> &_adj_par_name_vec, const Parameters &_fixed__ctl_pars,
	const Parameters &_lower_bnd, const Parameters &_upper_bnd, ParamTransformSeq *_base_partran_seq, int _n_sample)
	: n_sample(_n_sample), adj_par_name_vec(_adj_par_name_vec), fixed_ctl_pars(_fixed__ctl_pars), lower_bnd(_lower_bnd),
	upper_bnd(_upper_bnd), run_manager_ptr(_run_manager_ptr), base_partran_seq(_base_partran_seq)
	{
	}
VectorXd Sobol::gen_uniform_rand_vec(long nsample, double min, double max)
{
  VectorXd v(nsample);
  double dy_2 = (max - min) / 2.0;
  double y_avg = (max + min) / 2.0;
  v.setRandom().array() =  v.setRandom().array() * dy_2 + y_avg ;
  return v;
}

void Sobol::gen_m1_m2()
{
	long npar = adj_par_name_vec.size();
	//generate random matrices
	double par_min;
	double par_max;
	VectorXd v1;
	m1 = MatrixXd::Zero(n_sample, npar);
	m2 = MatrixXd::Zero(n_sample, npar);
	for (int i=0; i<npar; ++i)
	{
		par_min = lower_bnd[adj_par_name_vec[i]];
		par_max = upper_bnd[adj_par_name_vec[i]];
		v1 = gen_uniform_rand_vec(n_sample, par_min, par_max);
		m1.col(i) = v1;
		v1 = gen_uniform_rand_vec(n_sample, par_min, par_max);
		m2.col(i) = v1;
	}
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

void Sobol::add_model_runs(const MatrixXd &n)
{
	for (int i=0; i<n_sample; ++i)
	{
		VectorXd tmp_vec =  n.row(i);
		Parameters tmp_pars(adj_par_name_vec, tmp_vec);
		tmp_pars.insert(fixed_ctl_pars.begin(), fixed_ctl_pars.end());
		base_partran_seq->ctl2model_ip(tmp_pars);
		run_manager_ptr->add_run(tmp_pars);
	}
}

void Sobol::assemble_runs(RunManagerAbstract &rm, ParamTransformSeq &base_partran_seq)
{
	MatrixXd n;
	run_manager_ptr->reinitialize();
	gen_m1_m2();
	
	//calculate a0
	int n_adj_par = adj_par_name_vec.size();

	add_model_runs(m1);

	//calculate first order runs a1,....an
	vector<int> idx_vec;
    for (int ai=0; ai<n_adj_par; ++ai)
	{
		idx_vec.clear();
		idx_vec.push_back(ai);
		n = gen_N_matrix(m1, m2, idx_vec);
		add_model_runs(n);
	}
	//calculate a012345..j
	add_model_runs(m2);
}

VectorXd Sobol::get_expected_value(int run_set)
{
	vector<string> par_name_vec = run_manager_ptr->get_par_name_vec();
	vector<string> obs_name_vec = run_manager_ptr->get_obs_name_vec();

	int npar = par_name_vec.size();
	int nobs = obs_name_vec.size();
	vector<double> par_value_vec;
	vector<double> obs_value_vec;
	par_value_vec.resize(npar, 0);
	obs_value_vec.resize(nobs, 0);

	int run_b = run_set * n_sample;
	int run_e = run_b + n_sample;
	Map<VectorXd> i_obs_vec(&obs_value_vec[0], nobs);
	VectorXd avg_vec = VectorXd::Zero(nobs);
	int nrun = 0;
	for(int run_id=run_b; run_id<run_e; ++run_id)
	{
		bool success = run_manager_ptr->get_run(run_id, &par_value_vec[0], npar, &obs_value_vec[0], nobs);
		if (success)
		{
			nrun++;
			avg_vec += i_obs_vec;
		}
	}
	avg_vec /= nrun;
	return avg_vec;
}

void Sobol::calc_sen()
{
	get_expected_value(3);
	//cout << get_expected_value(1) << endl << endl;
	//cout << get_expected_value(0) - get_expected_value(1) << endl;
}