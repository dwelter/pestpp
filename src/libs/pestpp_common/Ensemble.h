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
	pair<int, int> shape() { return pair<int, int>(reals.rows(), reals.cols()); }
	void throw_ensemble_error(string message);
	void throw_ensemble_error(string message,vector<string> vec);
	const vector<string> get_var_names() const { return var_names; }
	const vector<string> get_real_names() const { return real_names; }

	Eigen::VectorXd get_real_vector(int ireal);
	Eigen::VectorXd get_real_vector(const string &real_name);
	Eigen::VectorXd get_var_values(int icol);
	Eigen::VectorXd get_var_values(const string &var_name);
	const Eigen::MatrixXd* eptr() const { return &reals; }
	Eigen::MatrixXd get_eigen(vector<string> row_names, vector<string> col_names);
	const Eigen::MatrixXd get_reals() const { return reals; }


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
	void from_csv(int num_reals,ifstream &csv);
	vector<string> prepare_csv(vector<string> names, ifstream &csv, bool forgive);
};

class ParameterEnsemble : public Ensemble
{
	
public:
	enum transStatus { CTL, NUM, MODEL };
	ParameterEnsemble(const ParamTransformSeq &_par_transform, Pest &_pest_scenario,
		FileManager &_file_manager,OutputFileWriter &_output_file_writer,
		PerformanceLog *_performance_log, unsigned int seed = 1);
	void initialize_with_csv(string &file_name);

	void enforce_bounds();
	void to_csv(string &file_name);
	Pest* get_pest_scenario_ptr() { return &pest_scenario; }
	const transStatus get_trans_status() const { return tstat; }
	void set_trans_status(transStatus _tstat) { tstat = _tstat; }
	const ParamTransformSeq get_par_transform() const { return par_transform; }
	
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
	void initialize_with_csv(string &file_name);

private: 
	ObjectiveFunc *obj_func_ptr;
};

class EnsemblePair
{
public:
	EnsemblePair(ParameterEnsemble &_pe, ObservationEnsemble &_oe);

	void run(RunManagerAbstract *run_mgr_ptr);

private:
	ParameterEnsemble &pe;
	ObservationEnsemble &oe;
	vector<int> active_real_indices;

};
#endif 