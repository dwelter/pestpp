#ifndef TORNADO_PLOT_H_
#define TORNADO_PLOT_H_

#include "GsaAbstractBase.h"

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
class TornadoPlot : public GsaAbstractBase
{
public:
	TornadoPlot(Pest &_pest_scenario, FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
		const ParamTransformSeq &_par_transform, bool _calc_obs_sen);
	~TornadoPlot();

	void assemble_runs(RunManagerAbstract &run_manager);
	void tornado_calc(RunManagerAbstract &run_manager, ModelRun model_run, std::ofstream &fout, const string obs_name = "");
	void calc_sen(RunManagerAbstract &run_manager, ModelRun model_run);
private:
	bool calc_obs_sen;
	Parameters init_pars;
	const ObservationInfo *obs_info_ptr;
	Parameters max_ctl_pars;
	Parameters min_ctl_pars;
};
#endif //TORNADO_PLOT_H_

