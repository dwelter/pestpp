#ifndef GSA_ABSTRACT_BASE_H_
#define GSA_ABSTRACT_BASE_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <random>
#include "Transformable.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ParamTransformSeq;
class RunManagerAbstract;
class ModelRun;
class FileManager;


class GsaAbstractBase
{
public:
	static const double MISSING_DATA;
	static mt19937_64 rand_engine;
	GsaAbstractBase(RunManagerAbstract *rm_ptr, ParamTransformSeq *base_partran_seq_ptr,
		const std::vector<std::string> &par_name_vec, const Parameters &_fixed_ctl_pars,
		const Parameters &_lower_bnd, const Parameters &_upper_bnd,
		const std::vector<std::string> &_obs_name_vec, FileManager *_file_manager_ptr);
	virtual void assemble_runs() = 0;
	virtual void calc_sen(ModelRun model_run, std::ofstream &fout_raw, std::ofstream &fout_sen, std::ofstream &fout_orw) = 0;
	static std::map<std::string, std::string>  process_gsa_file(std::ifstream &fin, FileManager &file_manager);
	std::string log_name(const string &name) const;
	//std::vector<double> calc_reduced_normal_intervals(int n_interval, double min, double max);
	std::vector<double> calc_interval_midpoints(int n_interval, double min, double max);
	double ltqnorm(double p);
	double mean_sq(const std::vector<double> &ao);
	double mean_sq_1order(const std::vector<double> &a0, const std::vector<double> &a1);
	double comp_var(const std::vector<double> &a0, const std::vector<double> &a1, double sample_mean);
	std::pair<double, double> max_min(const std::vector<double> &ao);
	virtual ~GsaAbstractBase(void);

protected:
	std::vector<std::string> adj_par_name_vec;
	std::vector<std::string> obs_name_vec;
	Parameters fixed_ctl_pars;
	Parameters lower_bnd;
	Parameters upper_bnd;
	std::set<std::string> log_trans_pars;
	RunManagerAbstract *run_manager_ptr;
	ParamTransformSeq *base_partran_seq_ptr;
	FileManager *file_manager_ptr;
	static void parce_line(const std::string &line, std::map<std::string, std::string> &arg_map);
};

std::ostream& operator<< (std::ostream& out, const std::vector<double> &rhs);
std::ostream& operator<< (std::ostream& out, const std::pair<double, double> &pr);

#endif /* GSA_ABSTRACT_BASE_H_ */
