#ifndef DIFFERENTIALEVOLUTION_H_
#define DIFFERENTIALEVOLUTION_H_

#include <unordered_map>
#include <random>
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunStorage.h"

class Pest;
class RunManagerAbstract;
class RestartController;


class ParameterInfoDE
{
public:
	ParameterInfoDE(double _lower_bnd = 0, double _upper_bnd = 0, bool _log_transform = false);
	ParameterInfoDE(const ParameterInfoDE &rhs);
	ParameterInfoDE& operator=(const ParameterInfoDE &rhs);
	~ParameterInfoDE();
public:
	double lower_bnd;
	double upper_bnd;
	bool log_transform;
};

class DifferentialEvolution
{
public:
	static mt19937_64 rand_engine;
	DifferentialEvolution(Pest &_pest_scenario, FileManager &_file_manager, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, OutputFileWriter &_output_file_writer,
		PerformanceLog *_performance_log, unsigned int seed = 1);
	void initialize_population(RunManagerAbstract &run_manager, int d);
	void solve(RunManagerAbstract &run_manager, RestartController &restart_controller,
		int max_gen, double f, double cr, ModelRun &cur_run);
	~DifferentialEvolution();
private:
	FileManager &file_manager;
	ObjectiveFunc *obj_func_ptr;
	std::unordered_map<std::string, ParameterInfoDE> parameter_info;
	const ParameterInfo *ctl_par_info_ptr;
	const ParameterGroupInfo *par_group_info_ptr;
	ParamTransformSeq par_transform;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	const ObservationInfo *obs_info_ptr;
	const PriorInformation *prior_info_ptr;
	std::vector<std::string> par_list;
	RunStorage gen_1;
	double f;
	int d;
	int best_run_idx;

	void initialize_vector(Parameters &ctl_pars);
	void mutation(RunManagerAbstract &run_manager, double f, double cr);
	int recombination(RunManagerAbstract &run_manager);
};

#endif //DIFFERENTIALEVOLUTION_H_