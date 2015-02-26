#ifndef TORNADO_H_
#define TORNADO_H_

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
#include "GsaAbstractBase.h"
class Tornado : public GsaAbstractBase
{
public:
	Tornado(const std::vector<std::string> &_adj_par_name_vec, const Parameters &_fixed_ctl_pars, const Parameters &lower_bnd,
		const Parameters &upper_bnd, const Parameters &inital_pars,
		ParamTransformSeq *base_partran_seq, const std::vector<std::string> &_obs_name_vec, FileManager *_file_manager_ptr);
	void assemble_runs(RunManagerAbstract &run_manager);
	//void calc_sen(RunManagerAbstract &run_manager, ModelRun model_run);
	virtual ~Tornado();
private:
	Parameters inital_pars;
};

#endif /* TORNADO_H_ */