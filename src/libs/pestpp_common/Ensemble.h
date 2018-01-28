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

	void to_csv(string file_name);
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

	void reorder(const vector<string> &_real_names, const vector<string> &_var_names);
	void drop_rows(const vector<int> &row_idxs);
	void keep_rows(const vector<int> &row_idxs);
	Pest* get_pest_scenario_ptr() { return pest_scenario_ptr; }
	Pest get_pest_scenario() { return *pest_scenario_ptr; }
	void set_pest_scenario(Pest *_pest_scenario) { pest_scenario_ptr = _pest_scenario; }
	void set_real_names(vector<string> &_real_names);
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
	void read_csv(int num_reals,ifstream &csv);
	vector<string> prepare_csv(const vector<string> &names, ifstream &csv, bool forgive);
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

	
	void from_csv(string file_name,const vector<string> &ordered_names);
	void from_csv(string file_name);
	void from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names,
		transStatus _tstat = transStatus::NUM);
	void enforce_bounds();
	void to_csv(string file_name);
	//Pest* get_pest_scenario_ptr() { return &pest_scenario; }
	transStatus get_trans_status() const { return tstat; }
	void set_trans_status(transStatus _tstat) { tstat = _tstat; }
	ParamTransformSeq get_par_transform() const { return par_transform; }
	void transform_ip(transStatus to_tstat);

	map<int,int> add_runs(RunManagerAbstract *run_mgr_ptr,const vector<int> &real_idxs=vector<int>());

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
	void from_csv(string file_name, const vector<string> &ordered_names);
	void from_csv(string file_name);
	void from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names);
	vector<int> update_from_runs(map<int,int> &real_run_ids, RunManagerAbstract *run_mgr_ptr);
	//ObservationEnsemble get_mean_diff();
};

class EnsemblePair
{
public:
	EnsemblePair(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	void queue_runs(RunManagerAbstract *run_mgr_ptr);
	void run(RunManagerAbstract *run_mgr_ptr);
	vector<int> process_runs(RunManagerAbstract *run_mgr_ptr);
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
	int num_active() { return active_real_indices.size(); }
	//void set_pe(ParameterEnsemble &_pe);
	//void set_oe(ObservationEnsemble &_oe);
	map<int, int> get_real_run_ids() { return real_run_ids; }
	void set_real_run_ids(map<int, int> &_real_run_ids) { real_run_ids = _real_run_ids; }

private:
	ParameterEnsemble &pe;
	ObservationEnsemble &oe;
	vector<int> active_real_indices;
	map<int, int> real_run_ids;
};
#endif 