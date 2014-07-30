#ifndef MORRISMETHOD_H_
#define MORRISMETHOD_H_

#include <vector>
#include <string>
#include <map>
#include <set>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include "GsaAbstractBase.h"
#include "Transformable.h"
#include "GsaAbstractBase.h"
#include "pest_data_structs.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class ParamTransformSeq;
class RunManagerAbstract;
class ModelRun;
class FileManager;
class RunningStats;

class MorrisObsSenFile
{
public:
	void initialize(const std::vector<std::string> &par_names_vec, const std::vector<std::string> &obs_names_vec, double _no_data, const GsaAbstractBase *_gsa_abstract_base);
	void add_sen_run_pair(const std::string &par_name, double p1, Observations &obs1, double p2, Observations &obs2);
	void calc_pooled_obs_sen(std::ofstream &fout_obs_sen, map<string, double> &obs_2_sen_weight, map<string, double> &par_2_sen_weight);
private:
	double no_data;
	vector<string> par_names_vec;
	vector<string> obs_names_vec;
	const GsaAbstractBase *gsa_abstract_base;
	map< pair<string,string>, RunningStats> map_obs_stats;
	map<string, int> parname_to_indexmap;
};

class MorrisMethod : public GsaAbstractBase
{
public:
	MorrisMethod(const std::vector<std::string> &_adj_par_name_vec, const Parameters &_fixed_ctl_pars, const Parameters &lower_bnd,
		const Parameters &upper_bnd, const set<string> &_log_trans_pars, int _p, int _r, 
		ParamTransformSeq *base_partran_seq,
		const std::vector<std::string> &_obs_name_vec, FileManager *_file_manager_ptr,
		const ObservationInfo *_obs_info_ptr, double norm, bool _calc_pooled_obs, double _delta);
	void process_pooled_var_file();
	void initialize(const set<string> &_log_trans_pars, int _p, int _r, double _delta);
	void assemble_runs(RunManagerAbstract &run_manager);
	void calc_sen(RunManagerAbstract &run_manager, ModelRun model_run);
	~MorrisMethod(void);
private:
	int k; // number of parameters
	int p; // number of levels for each parameters
	int r;
	bool calc_obs_sen;
	double delta;
	double norm;
	MatrixXd b_star_mat;
	static MatrixXd create_B_mat(int k);
	static MatrixXd create_J_mat(int k);
	static MatrixXd create_D_mat(int k);
	static MatrixXd create_P_mat(int k);
		MatrixXd create_P_star_mat(int k);
	VectorXd create_x_vec(int k);
	Parameters get_ctl_parameters(int row);
	static int rand_plus_minus_1(void);
	MorrisObsSenFile obs_sen_file;
	const ObservationInfo *obs_info_ptr;
	map<std::string, std::string> group_2_pool_group_map;
};




#endif /* MORRISMETHOD_H_ */
