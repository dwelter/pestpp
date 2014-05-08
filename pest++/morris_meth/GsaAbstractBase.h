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
	enum class PARAM_DIST{ normal, uniform };
	static const double MISSING_DATA;
	static mt19937_64 rand_engine;
	GsaAbstractBase(RunManagerAbstract *rm_ptr, ParamTransformSeq *base_partran_seq_ptr,
		const std::vector<std::string> &par_name_vec, const Parameters &_fixed_ctl_pars,
		const Parameters &_lower_bnd, const Parameters &_upper_bnd,
		const std::vector<std::string> &_obs_name_vec, FileManager *_file_manager_ptr,
		PARAM_DIST _par_dist = PARAM_DIST::normal);
	virtual void assemble_runs() = 0;
	virtual void calc_sen(ModelRun model_run) = 0;
	static std::map<std::string, std::string>  process_gsa_file(std::ifstream &fin, FileManager &file_manager);
	std::string log_name(const string &name) const;
	bool is_log_trans_par(const string &name) const;
	std::vector<double> calc_interval_midpoints(int n_interval, double min, double max);
	double ltqnorm(double p);
	virtual ~GsaAbstractBase(void);

protected:
	PARAM_DIST par_dist;
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
	map<string, double> calc_parameter_norm_std_dev();
	map<string, double> calc_parameter_unif_std_dev();
};

std::ostream& operator<< (std::ostream& out, const std::vector<double> &rhs);
std::ostream& operator<< (std::ostream& out, const std::pair<double, double> &pr);

#endif /* GSA_ABSTRACT_BASE_H_ */
