/*  
    © Copyright 2012, David Welter
    
    This file is part of PEST++.
   
    PEST++ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PEST++ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/

#ifndef PEST_H_
#define PEST_H_

#include <string>
#include <map>
#include <set>
#include <vector>
#include "ParamTransformSeq.h"
#include "pest_data_structs.h"
#include "PriorInformation.h"
#include "Regularization.h"

using namespace std;
class FileManager;
class PriorInformation;

class Pest {
public:
	friend ostream& operator<< (ostream &os, const Pest& val);
	Pest();
	void set_defaults();
	void check_inputs();
	int process_ctl_file(ifstream &fin, FileManager &file_manager);
	const Parameters& get_ctl_parameters() const {return ctl_parameters;}
	const Observations& get_ctl_observations() const {return observation_values;}
	const ParameterInfo& get_ctl_parameter_info()const {return ctl_parameter_info;}
	const ParameterGroupInfo& get_base_group_info() const {return  base_group_info;}
	const ObservationInfo &get_ctl_observation_info() const {return observation_info;}
	const PriorInformation &get_prior_info() {return prior_info;}
	const PriorInformation *get_prior_info_ptr() {return &prior_info;}
	const SVDInfo& get_svd_info() const {return svd_info;}
	const ControlInfo&  get_control_info() const {return control_info;}
	const ParamTransformSeq& get_base_par_tran_seq() const {return base_par_transform;}
	const vector<string> &get_ctl_ordered_par_names() {return ctl_ordered_par_names;}
	const vector<string> &get_ctl_ordered_obs_names() {return ctl_ordered_obs_names;}
	const ModelExecInfo &get_model_exec_info() {return model_exec_info;}
	const  vector<string> &get_comline_vec();
	const  vector<string> &get_tplfile_vec();
	const  vector<string> &get_inpfile_vec();
	const  vector<string> &get_insfile_vec();
	const  vector<string> &get_outfile_vec();
	const PestppOptions &get_pestpp_options() const {return pestpp_options;}
	const Regularization* get_regul_scheme_ptr() {return regul_scheme_ptr;}
	vector<string> get_nonregul_obs() const;
	virtual ~Pest();
private:
	ControlInfo control_info;
	SVDInfo svd_info;
	Parameters ctl_parameters;
	ParameterInfo ctl_parameter_info;
	ParameterGroupInfo base_group_info;
	Observations observation_values;
	ObservationInfo observation_info;
	PriorInformation prior_info;
	ModelExecInfo model_exec_info;
	PestppOptions pestpp_options;
	ParamTransformSeq base_par_transform;
	vector<string> ctl_ordered_par_names;
	vector<string> ctl_ordered_obs_names;
	Regularization *regul_scheme_ptr;
};
ostream& operator<< (ostream &os, const Pest& val);
#endif /* PEST_H_ */
