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

#ifndef JACOBIAN_H_
#define JACOBIAN_H_
#include<map>
#include<unordered_map>
#include<vector>
#include<set>
#include<list>
#include<Eigen\Dense>
#include "Transformable.h"

class ParamTransformSeq;
class ParameterInfo;
class ParameterGroupInfo;
class RunManagerAbstract;
class ObjectiveFunc;
class ModelRun;
class FileManager;
class PriorInformation;

class Jacobian {
protected:
	class JacobianRun;

public:
	Jacobian(FileManager &_file_manager);
	virtual const vector<string>& parameter_list() const{return base_numeric_par_names;}
	virtual const vector<string>& observation_list() const {return  base_sim_obs_names;}
	virtual vector<string> obs_and_reg_list() const;
	virtual const Parameters &get_base_numeric_parameters() const{return base_numeric_parameters;};
	virtual const Observations &get_base_sim_observations() const {return base_sim_observations;}
	virtual Eigen::MatrixXd get_matrix(const vector<string> & par_name_vec, const vector<string> &obs_names) const;
	virtual void calculate(ModelRun &model_run, vector<string> numeric_par_names, vector<string> obs_names, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, bool phiredswh_flag=false, bool calc_init_obs=true);
	virtual void calculate(ModelRun &model_run, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, bool phiredswh_flag=false, bool calc_init_obs=true);
	virtual void save(const std::string &ext="jco") const;
	virtual const set<string>& get_failed_parameter_names() const;
	virtual ~Jacobian();
protected:
	vector<string> base_numeric_par_names;  //ordered names of base parameters used to calculate the jacobian
	Parameters base_numeric_parameters;  //values of base parameters used to calculate the jacobian
	set<string> failed_parameter_names;
	vector< string>  base_sim_obs_names;  //names of base observations used to calculate the jacobian
	Observations  base_sim_observations;  //values of base observations used to calculate the jacobian
	Eigen::MatrixXd matrix;
	map<string, map<string, double>> prior_info_sen;
	FileManager &file_manager;  // filemanger used to get name of jaobian file

	virtual void calc_derivative(const string &numeric_par_name, int jcol, list<ModelRun> &run_list,  const ParameterGroupInfo &group_info,
		const ParameterInfo &ctl_par_info, const PriorInformation &prior_info);
	virtual bool forward_diff(const string &par_name, const Parameters &pest_parameters, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, double &new_par);
	virtual bool central_diff(const string &par_name, const Parameters &pest_parameters, 
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, vector<double> &new_par, 
		vector<Parameters> &model_par_vec);
	virtual bool out_of_bounds(const Parameters &model_parameters, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, vector<string> &out_of_bound_par_vec) const;
	virtual double derivative_inc(const string &name, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,   double cur_par_value,  bool central = false);
	virtual bool get_derivative_parameters(const string &par_name, ModelRun &init_model_run, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, 
		vector<JacobianRun> &del_numeric_par_vec, bool phiredswh_flag);
	virtual void calc_prior_info_sen(const string &par_name, ModelRun &run1, ModelRun &run2, const PriorInformation &prior_info);
	virtual int size_prior_info_sen() const;
	virtual unordered_map<string, int> get_par2col_map() const;
	virtual unordered_map<string, int> get_obs2row_map() const;
	
	class JacobianRun {
	public:
		JacobianRun(const string &_par_name, double _numeric_value) : par_name(_par_name), numeric_value(_numeric_value) {}
		JacobianRun(const JacobianRun &rhs) : par_name(rhs.par_name), numeric_value(rhs.numeric_value) {}
		~JacobianRun(){}
		JacobianRun & operator= (const JacobianRun &rhs) {par_name=rhs.par_name; numeric_value=rhs.numeric_value; return *this;}
		string par_name;
		double numeric_value;
	};
};

#endif /* JACOBIAN_H_ */
