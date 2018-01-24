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
	//Ensemble(Pest &_pest_scenario, FileManager &_file_manager,
	//	OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, unsigned int seed = 1);
	Ensemble(Pest* _pest_scenario);
	Ensemble() { ; }
	
	//Ensemble get(vector<string> &_real_names, vector<string> &_var_names);

	void to_csv(string &file_name);
	void from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names);
	pair<int, int> shape() { return pair<int, int>(reals.rows(), reals.cols()); }
	void throw_ensemble_error(string message);
	void throw_ensemble_error(string message,vector<string> vec);
	const vector<string> get_var_names() const { return var_names; }
	const vector<string> get_real_names() const { return real_names; }

	const vector<string> get_real_names(vector<int> &indices);

	void add_to_cols(Eigen::MatrixXd &_reals, const vector<string> &_var_names);

	

	Eigen::VectorXd get_real_vector(int ireal);
	Eigen::VectorXd get_real_vector(const string &real_name);

	Eigen::MatrixXd get_eigen(vector<string> row_names, vector<string> col_names);
	const Eigen::MatrixXd get_eigen() const { return reals; }
	const Eigen::MatrixXd* get_eigen_ptr() const { return &reals; }
	void set_eigen(Eigen::MatrixXd _reals);

	Eigen::MatrixXd get_eigen_mean_diff();
	Eigen::MatrixXd get_eigen_mean_diff(const vector<string> &_real_names, const vector<string> &_var_names);
	
	void append_other_rows(Ensemble &other);
	void append(string real_name, Transformable &trans);
	Covariance get_diagonal_cov_matrix();

	void reorder(vector<string> &_real_names, vector<string> &_var_names);
	void drop_rows(vector<int> &row_idxs);
	void keep_rows(vector<int> &row_idxs);
	Pest* get_pest_scenario_ptr() { return pest_scenario_ptr; }
	Pest get_pest_scenario() { return *pest_scenario_ptr; }
	void set_pest_scenario(Pest *_pest_scenario) { pest_scenario_ptr = _pest_scenario; }
	void set_real_names(vector<string> &_real_names);

	void draw(int num_reals, Covariance &cov, Transformable &tran, const vector<string> &draw_names);
	~Ensemble();
protected:
	Pest* pest_scenario_ptr;
	//FileManager &file_manager;
	//ObjectiveFunc *obj_func_ptr;
	//OutputFileWriter &output_file_writer;
	//PerformanceLog *performance_log;

	Eigen::MatrixXd reals;
	vector<string> var_names;
	vector<string> real_names;	
	void read_csv(int num_reals,ifstream &csv, map<string,int> header_info);
	void from_binary(string &file_name, vector<string> &names,  bool transposed);
	map<string,int> prepare_csv(const vector<string> &names, ifstream &csv, bool forgive);
};

class ParameterEnsemble : public Ensemble
{
	
public:
	enum transStatus { CTL, NUM, MODEL };
	/*ParameterEnsemble(const ParamTransformSeq &_par_transform, Pest &_pest_scenario,
		FileManager &_file_manager,OutputFileWriter &_output_file_writer,
		PerformanceLog *_performance_log, unsigned int seed = 1);
	*/
	ParameterEnsemble(Pest *_pest_scenario_ptr);
	ParameterEnsemble() { ; }
	
	//ParameterEnsemble get_new(const vector<string> &_real_names, const vector<string> &_var_names);

	
	//void from_csv(string &file_name,const vector<string> &ordered_names);
	void from_csv(string &file_name);
	void from_binary(string &file_name);// { Ensemble::from_binary(file_name, false); }
	void from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names,
		transStatus _tstat = transStatus::NUM);
	void enforce_bounds();
	void to_csv(string &file_name);
	//Pest* get_pest_scenario_ptr() { return &pest_scenario; }
	const transStatus get_trans_status() const { return tstat; }
	void set_trans_status(transStatus _tstat) { tstat = _tstat; }
	const ParamTransformSeq get_par_transform() const { return par_transform; }
	void transform_ip(transStatus to_tstat);
	void set_pest_scenario(Pest *_pest_scenario);
	map<int,int> add_runs(RunManagerAbstract *run_mgr_ptr,vector<int> &real_idxs=vector<int>());

	void draw(int num_reals, Covariance &cov);
	Covariance get_diagonal_cov_matrix();
	//ParameterEnsemble get_mean_diff();
private:
	ParamTransformSeq par_transform;
	transStatus tstat;
};

class ObservationEnsemble : public Ensemble
{
public: 
	/*ObservationEnsemble(ObjectiveFunc *_obj_func, Pest &_pest_scenario, FileManager &_file_manager,
    OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, unsigned int seed = 1);
	*/
	ObservationEnsemble(Pest *_pest_scenario_ptr);
	ObservationEnsemble() { ; }
	void update_from_obs(int row_idx, Observations &obs);
	void update_from_obs(string real_name, Observations &obs);
	//void from_csv(string &file_name, const vector<string> &ordered_names);
	void from_csv(string &file_name);
	void from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names);
	void from_binary(string &file_name);// { Ensemble::from_binary(file_name, true); }
	vector<int> update_from_runs(map<int,int> &real_run_ids, RunManagerAbstract *run_mgr_ptr);
	void draw(int num_reals, Covariance &cov);

	//ObservationEnsemble get_mean_diff();
};


#endif 