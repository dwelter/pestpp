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
	Mat to_matrix(vector<string> row_names, vector<string> col_names);

	void write_csv(string file_name);
	void initialize_with_csv(string &file_name);

	pair<int, int> shape();
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
	
	void from_csv(ifstream &csv);
	vector<string> prepare_csv(vector<string> names, ifstream &csv, bool forgive);

};

class ParameterEnsemble : private Ensemble
{
public:
	ParameterEnsemble(Pest &_pest_scenario, FileManager &_file_manager,
		const ParamTransformSeq &_par_transform, OutputFileWriter &_output_file_writer,
		PerformanceLog *_performance_log, unsigned int seed = 1);
	void enforce_bounds();
	void initialize_with_draws(int num_reals);

private:
	ParamTransformSeq par_transform;
};

class ObservationEnsemble : private Ensemble
{
public: 
	ObservationEnsemble(Pest &_pest_scenario, FileManager &_file_manager, ObjectiveFunc *_obj_func,
    OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, unsigned int seed = 1);
	bool update_from_runs(RunManagerAbstract *run_mgmt_ptr);
	void initialize_with_draws(int num_reals);

private: 
	ObjectiveFunc *obj_func_ptr;
};

#endif 