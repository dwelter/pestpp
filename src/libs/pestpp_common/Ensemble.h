#ifndef ENSEMBLE_H_
#define ENSEMBLE_H_

#include <map>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunStorage.h"
#include "covariance.h"
#include "RunManagerAbstract.h"



class Ensemble
{
public:
	static mt19937_64 rand_engine;
	Ensemble(Pest &_pest_scenario, FileManager &_file_manager,
		OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, unsigned int seed = 1);
	
	Mat to_matrix();
	Mat to_matrix(vector<string> &row_names, vector<string> &col_names);

	void to_csv(string &file_name);
	void from_eigen_mat(Eigen::MatrixXd mat, vector<string> real_names, vector<string> var_names);
	pair<int, int> shape() { return pair<int, int>(reals.rows(), reals.cols()); }
	void throw_ensemble_error(string message);
	void throw_ensemble_error(string message,vector<string> vec);
	const vector<string> get_var_names() const { return var_names; }
	const vector<string> get_real_names() const { return real_names; }

	Eigen::VectorXd get_real_vector(int ireal);
	Eigen::VectorXd get_real_vector(const string &real_name);
	Eigen::MatrixXd get_eigen(vector<string> row_names, vector<string> col_names);

	const Eigen::MatrixXd get_eigen() const { return reals; }
	const Eigen::MatrixXd* get_eigen_ptr() const { return &reals; }
	void set_reals(Eigen::MatrixXd _reals);

	Eigen::MatrixXd get_eigen_mean_diff();
	Eigen::MatrixXd get_eigen_mean_diff(vector<string> &_real_names);
	
	void reorder(vector<string> &_real_names, vector<string> *_var_names);

	~Ensemble();
protected:
	Pest &pest_scenario;
	FileManager &file_manager;
	ObjectiveFunc *obj_func_ptr;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	Eigen::MatrixXd reals;
	vector<string> var_names;
	vector<string> real_names;	
	void read_csv(int num_reals,ifstream &csv);
	vector<string> prepare_csv(const vector<string> &names, ifstream &csv, bool forgive);
};

class ParameterEnsemble : public Ensemble
{
	
public:
	enum transStatus { CTL, NUM, MODEL };
	ParameterEnsemble(const ParamTransformSeq &_par_transform, Pest &_pest_scenario,
		FileManager &_file_manager,OutputFileWriter &_output_file_writer,
		PerformanceLog *_performance_log, unsigned int seed = 1);
	void from_csv(string &file_name);

	void enforce_bounds();
	void to_csv(string &file_name);
	Pest* get_pest_scenario_ptr() { return &pest_scenario; }
	const transStatus get_trans_status() const { return tstat; }
	void set_trans_status(transStatus _tstat) { tstat = _tstat; }
	const ParamTransformSeq get_par_transform() const { return par_transform; }
	void transform_ip(transStatus to_tstat);
	//ParameterEnsemble get_mean_diff();
private:
	ParamTransformSeq par_transform;
	transStatus tstat;
};

class ObservationEnsemble : public Ensemble
{
public: 
	ObservationEnsemble(ObjectiveFunc *_obj_func, Pest &_pest_scenario, FileManager &_file_manager,
    OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, unsigned int seed = 1);
	void update_from_obs(int row_idx, Observations &obs);
	void from_csv(string &file_name);
	vector<double> get_phi_vec();
	//ObservationEnsemble get_mean_diff();
private: 
	ObjectiveFunc *obj_func_ptr;
};

class EnsemblePair
{
public:
	EnsemblePair(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	void queue_runs(RunManagerAbstract *run_mgr_ptr);
	void run(RunManagerAbstract *run_mgr_ptr);
	void process_runs(RunManagerAbstract *run_mgr_ptr);
	/*Eigen::MatrixXd get_active_oe_eigen();
	Eigen::MatrixXd get_active_pe_eigen();
	Eigen::MatrixXd get_active_oe_mean_diff();
	Eigen::MatrixXd get_active_pe_mean_diff();*/

	ParameterEnsemble* get_pe_ptr() { return &pe; }
	ObservationEnsemble* get_oe_ptr() { return &oe; }
	const vector<int> get_act_real_indices() const { return active_real_indices; }
	void set_act_real_indices(vector<int> _act_real_indices) { active_real_indices = _act_real_indices; }
	vector<string> get_pe_active_names();
	vector<string> get_oe_active_names();
	//EnsemblePair get_mean_diff();

private:
	ParameterEnsemble &pe;
	ObservationEnsemble &oe;
	vector<int> active_real_indices;
	map<int, int> real_run_ids;
};
#endif 