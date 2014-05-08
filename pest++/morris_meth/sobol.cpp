#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <iterator>
#include "sobol.h"
#include "Transformable.h"
#include "RunManagerAbstract.h"
#include "ParamTransformSeq.h"
#include "ModelRunPP.h"
#include "Stats.h"

using namespace std;
using namespace Eigen;


Sobol::Sobol(const vector<string> &_adj_par_name_vec, const Parameters &_fixed__ctl_pars,
	const Parameters &_lower_bnd, const Parameters &_upper_bnd, int _n_sample,
	RunManagerAbstract *_rm_ptr, ParamTransformSeq *base_partran_seq_ptr,
	const std::vector<std::string> &_obs_name_vec, FileManager *_file_manager_ptr)
	: GsaAbstractBase(_rm_ptr, base_partran_seq_ptr, _adj_par_name_vec, 
	       _fixed__ctl_pars, _lower_bnd, _upper_bnd, _obs_name_vec, _file_manager_ptr), n_sample(_n_sample)
	{
	}
VectorXd Sobol::gen_rand_vec(long nsample, double min, double max)
{
	//std::uniform_real_distribution<double> distribution(min, max);
	std::normal_distribution<> distribution((max+min)/ 2.0 , (max-min)/4.0);
	VectorXd v(nsample);
	long v_len = v.size();
	for (long i = 0; i<v_len; ++i)
	{
		v[i] = distribution(rand_engine);
	}
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
		v1 = gen_rand_vec(n_sample, par_min, par_max);
		m1.col(i) = v1;
	}
	for (int i=0; i<npar; ++i)
	{
		par_min = lower_bnd[adj_par_name_vec[i]];
		par_max = upper_bnd[adj_par_name_vec[i]];
		v1 = gen_rand_vec(n_sample, par_min, par_max);
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
		base_partran_seq_ptr->ctl2model_ip(tmp_pars);
		run_manager_ptr->add_run(tmp_pars);
	}
}

void Sobol::assemble_runs()
{
	MatrixXd c;
	run_manager_ptr->reinitialize();
	gen_m1_m2();
	
	//calculate a0
	int n_adj_par = adj_par_name_vec.size();

	add_model_runs(m1);
	add_model_runs(m2);

	//cout << m1 << endl << endl;
	//cout << m2 << endl << endl;
	//calculate first order runs a1,....an
	vector<int> idx_vec;
	for (int ai=0; ai<n_adj_par; ++ai)
	{
		idx_vec.clear();

		idx_vec.push_back(ai);
		c = gen_N_matrix(m1, m2, idx_vec);
		//cout << c << endl << endl;
		add_model_runs(c);
	}
}

vector<double> Sobol::get_phi_vec(int run_set, ModelRun &model_run)
{
	ModelRun run0 = model_run;

	int run_b = run_set * n_sample;
	int run_e = run_b + n_sample;

	Parameters pars0;
	Observations obs0;
	int nrun = 0;
	vector<double> phi_vec = vector<double>(n_sample, MISSING_DATA);
	for(int run_id=run_b; run_id<run_e; ++run_id)
	{
		double phi = MISSING_DATA;
		bool success = run_manager_ptr->get_run(run_id, pars0, obs0);
		if (success)
		{
			run0.update_ctl(pars0, obs0);
			phi = run0.get_phi(0.0);
		}
		phi_vec[nrun] = phi;
		nrun++;
	}
	return phi_vec;
}

void Sobol::calc_sen(ModelRun model_run)
{
	vector<double> ya = get_phi_vec(0, model_run);
	vector<double> yb = get_phi_vec(1, model_run);
	
	//cout <<"######### A " << endl;
	//cout << "------YA--------" << endl;
	//cout << "pars" << endl;
	//cout << run_manager_ptr->get_model_parameters(0) << endl;
	//cout << "phi" << endl;
	//cout << ya << endl;
	//cout << "------YA_end--------" << endl;
	//cout <<"######### B " << endl;
	//cout << "------YB--------" << endl;
	//cout << yb << endl;
	//cout << "------YB_end--------" << endl;
	auto ya_stats = vec_calc_stats_missing_data(ya, MISSING_DATA);
	double mean_ya = ya_stats["mean"];
	double var_ya = sobol_u_missing_data(ya, ya, MISSING_DATA) - pow(mean_ya, 2.0);
	 
	cout << endl;
	cout << "mean_ya=" << mean_ya << endl;
	cout << "var_ya=" << var_ya << endl;
	int npar = adj_par_name_vec.size();
	for (int i=0; i<npar; ++i)
	{
		cout << "nruns = " << run_manager_ptr->get_nruns() << endl;
		cout <<"######### PAR = " << adj_par_name_vec[i] << endl;
		vector<double> yci = get_phi_vec(i+2, model_run);

		//{
		//	auto covar_stats  =  vec_covar_missing_data(ya, yci, MISSING_DATA);
		//	double mean_ya = covar_stats["mean_1"];
		//	double var_ya = covar_stats["var_1"];
		//	double mean_yci = covar_stats["mean_2"];
		//	double var_yci = covar_stats["var_2"];
		//	double covar = covar_stats["covar"];
		//	cout << " - mean_ya=" << mean_ya << endl;
	 //       cout << " - var_ya=" << var_ya << endl;
		//	cout << " - mean_yci=" << mean_yci << endl;
	 //       cout << " - var_yci=" << var_yci << endl;
		//	cout << " - covar=" << covar << endl;
		//	double si = covar / var_ya;
		//	double st =  1 - covar / var_ya;
		//	cout << "pararmeter " << adj_par_name_vec[i]  << "   si=" << si << ",   st=" << st << endl;
		//    cout <<"######### END" << endl;
		//}

		double sobol_uj = sobol_u_missing_data(ya, yci, MISSING_DATA);
		cout << " U=" << (sobol_uj - pow(mean_ya, 2.0)) << "  yar_ya= " << var_ya  <<  endl;
		double si = (sobol_uj -  pow(mean_ya, 2.0)) / var_ya;

		double sobol_umj = sobol_u_missing_data(yb, yci, MISSING_DATA);
		double st =  1 - (sobol_umj - pow(mean_ya, 2.0)) / var_ya;

		cout << "pararmeter " << adj_par_name_vec[i]  << "   si=" << si << ",   st=" << st << endl;
		cout <<"######### END" << endl;
	}
}
