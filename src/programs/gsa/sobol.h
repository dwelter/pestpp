#ifndef SOBOL_H_
#define SOBOL_H_

#include <vector>
#include <string>
#include <Eigen/Dense>
#include "GsaAbstractBase.h"
#include "Transformable.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ParamTransformSeq;
class RunManagerAbstract;
class ModelRun;

class Sobol : public GsaAbstractBase
{
public:
	Sobol(const std::vector<std::string> &_adj_par_name_vec, const Parameters &_fixed__ctl_pars,
		const Parameters &_lower_bnd, const Parameters &_upper_bnd, int n_sample,
		ParamTransformSeq *base_partran_seq,
		const std::vector<std::string> &_obs_name_vec, FileManager *_file_manager_ptr, PARAM_DIST dist);
	void assemble_runs(RunManagerAbstract &run_manager);
	void calc_sen(RunManagerAbstract &run_manager, ModelRun model_run);
	void calc_sen_single(RunManagerAbstract &run_manager, ModelRun model_run, std::ofstream &fout_sbl, const std::string &obs_name);
private:
	VectorXd gen_rand_vec(long nsample, double min, double max, bool log_transform=false);
	void gen_m1_m2();
	MatrixXd gen_N_matrix(const MatrixXd &m1, const MatrixXd &m2, const vector<int> &idx_vec);
	void add_model_runs(RunManagerAbstract &run_manager, const MatrixXd &n);
	vector<double> get_obs_vec(RunManagerAbstract &run_manager, int run_set, ModelRun &model_run, const string &obs_name );
	vector<double> get_phi_vec(RunManagerAbstract &run_manager, int run_set, ModelRun &model_run);
	int n_sample;
	Eigen::MatrixXd m1;
	Eigen::MatrixXd m2;
};
#endif /* SOBOL_H_ */
