#include <vector>
#include <Eigen/Dense>
#include <cassert>
#include <numeric>
#include <math.h>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <regex>
#include "Tornado.h"
#include "Transformable.h"
#include "RunManagerAbstract.h"
#include "ParamTransformSeq.h"
#include "ModelRunPP.h"
#include "utilities.h"
#include "FileManager.h"
#include "Stats.h"

using namespace std;
using namespace pest_utils;
using Eigen::MatrixXd;
using Eigen::VectorXd;


Tornado::Tornado(const std::vector<std::string> &_adj_par_name_vec, const Parameters &_fixed_ctl_pars, const Parameters &_lower_bnd,
	const Parameters &_upper_bnd, const Parameters &_inital_pars,
	ParamTransformSeq *_base_partran_seq_ptr, const std::vector<std::string> &_obs_name_vec, FileManager *_file_manager_ptr)
	: GsaAbstractBase(_base_partran_seq_ptr, _adj_par_name_vec, _fixed_ctl_pars, _lower_bnd, _upper_bnd,
	_obs_name_vec, file_manager_ptr, PARAM_DIST::uniform), inital_pars(_inital_pars)
{
}

void Tornado::assemble_runs(RunManagerAbstract &run_manager)
{
	int run_id = -999;
	Parameters base_pars = base_partran_seq_ptr->ctl2model_cp(inital_pars);
	run_id = run_manager.add_run(base_pars, "base_run", Parameters::no_data);
	for (const auto &ipar : adj_par_name_vec)
	{
	    string run_name = ipar + "_low";
		Parameters pars = inital_pars;
		auto iter = lower_bnd.find(ipar);
		if (iter != lower_bnd.end())
		{
			pars[ipar] = iter->second;
		}
		else
		{
			ostringstream str;
			str << "Error Tornado::assemble runs can not find lower bound for  parameter: " << ipar;
			throw(str.str());
		}
		// convert control parameters to model parameters
		base_partran_seq_ptr->ctl2model_ip(pars);
		run_id = run_manager.add_run(pars, run_name, Parameters::no_data);
		 
		run_name = ipar + "_high";
		pars = inital_pars;
		iter = upper_bnd.find(ipar);
		if (iter != upper_bnd.end())
		{
			pars[ipar] = iter->second;
		}
		else
		{
			ostringstream str;
			str << "Error Tornado::assemble runs can not find upper bound for parameter: " << ipar;
			throw(str.str());
		}
		// convert control parameters to model parameters
		base_partran_seq_ptr->ctl2model_ip(pars);
		run_id = run_manager.add_run(pars, run_name, Parameters::no_data);
	}
}


Tornado::~Tornado()
{
}
