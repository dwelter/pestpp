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

#ifndef MODELRUNPP_H_
#define MODELRUNPP_H_

#include <list>
#include <utility>
#include <string>
#include "Transformable.h"
#include "ObjectiveFunc.h"
#include "ParamTransformSeq.h"

using namespace std;
class ModelRunShallowCopy;

class ModelRunAbstractBase
{
public:
	class Compare
	{
	public:
		enum FIELD{NUMERIC_PAR, MODEL_PAR, CTL_PAR, SIM_OBS};
		Compare::Compare(const string &_name, FIELD _field)
			: name(_name), field(_field) {}
		bool operator()(ModelRunAbstractBase &run1, ModelRunAbstractBase &run2);
	private:
		string name;
		FIELD field;
	};
public:
	ModelRunAbstractBase() {}
	virtual void set_numeric_parameters(const Parameters &pars) = 0;
	virtual void set_ctl_parameters(const Parameters &pars) = 0;
	virtual void set_model_parameters(const Parameters &pars) = 0;
	virtual void set_observations(const Observations &observations) = 0;
	virtual const Parameters &get_numeric_pars() = 0;
	virtual Parameters get_ctl_pars() = 0;
	virtual const Parameters &get_model_pars() = 0;
	virtual const Observations &get_obs() = 0;
	virtual Observations get_obs_template() const = 0;
	virtual const ParamTransformSeq &get_par_tran() = 0;
	virtual double get_phi(double regul_weight=1.0) = 0;
	virtual PhiComponets get_phi_comp() = 0;
	virtual vector<double> get_residuals_vec(const vector<string> &obs_names) = 0;
	virtual void phi_report(ostream &os) = 0;
	virtual bool obs_valid() const = 0;
	virtual bool numeric_pars_valid() const = 0;
	virtual bool model_pars_valid() const = 0;
	virtual bool phi_valid() const = 0;
	virtual ~ModelRunAbstractBase(){};
};
class ModelRun : public ModelRunAbstractBase
{
public:
	friend class ModelRunShallowCopy;
	enum CopyType {DeepCopy, ShallowCopy};
	//ModelRun(const ModelRun &rhs, CopyType=DeepCopy);
	ModelRun(const ObjectiveFunc *_objectiveFunc, const ParamTransformSeq &_par_tran, Observations _sim_obs);
	ModelRun(const ObjectiveFunc *_objectiveFunc, const ParamTransformSeq &_par_tran, const Parameters &_numeric_pars, const Observations &_sim_obs);
	ModelRun& operator=(const ModelRun &rhs);
	//const ModelRun& operator=(const ModelRunShallowCopy &rhs);
	virtual void set_numeric_parameters(const Parameters &pars);
	virtual void set_ctl_parameters(const Parameters &pars);
	virtual void set_model_parameters(const Parameters &pars);
	virtual void set_observations(const Observations &observations);
	virtual const Parameters &get_numeric_pars();
	virtual Parameters get_ctl_pars();
	virtual const Parameters &get_model_pars();
	virtual const Observations &get_obs();
	virtual Observations get_obs_template() const;
	virtual const ParamTransformSeq &get_par_tran()  {return par_tran;}
	const ObjectiveFunc *get_obj_func_ptr()  {return obj_func_ptr;}
	void freeze_parameters(const map<string,double> par_map);
	void thaw_parameters();
	virtual double get_phi(double regul_weight=1.0);
	virtual PhiComponets get_phi_comp();
	virtual vector<double> get_residuals_vec(const vector<string> &obs_names);
	virtual void phi_report(ostream &os);
	virtual bool obs_valid() const;
	virtual bool numeric_pars_valid() const;
	virtual bool model_pars_valid() const;
	virtual bool phi_valid() const;
	virtual ~ModelRun();
protected:
	PhiComponets phi_comp;
	bool numeric_pars_is_valid;
	bool model_pars_is_valid;
	bool obs_is_valid;
	bool phi_is_valid;
	const ObjectiveFunc *obj_func_ptr;
	Parameters numeric_pars;
	Parameters model_pars;
	Observations sim_obs;
private:
	ParamTransformSeq par_tran;
};

//class ModelRunShallowCopy  : public ModelRunAbstractBase
//{
//public:
//	friend class ModelRun;
//	ModelRunShallowCopy(const ModelRun &model_run);
//	//ModelRunShallowCopy(const ModelRunShallowCopy &model_run);
//	//const ModelRunShallowCopy& operator=(const ModelRunShallowCopy &rhs);
//	const ModelRunShallowCopy& ModelRunShallowCopy::operator=(const ModelRun &rhs);
//	virtual void set_numeric_parameters(const Parameters &pars);
//	virtual void set_ctl_parameters(const Parameters &pars);
//	virtual void set_model_parameters(const Parameters &pars);
//	virtual void set_observations(const Observations &observations);
//	virtual const Parameters &get_numeric_pars();
//	virtual Parameters get_ctl_pars();
//	virtual const Parameters &get_model_pars();
//	virtual const Observations &get_obs();
//	virtual Observations get_obs_template() const;
//	virtual const ParamTransformSeq &get_par_tran()  {return *(base_model_run_ptr->par_tran);}
//	virtual double get_phi(double regul_weight=1.0);
//	virtual PhiComponets get_phi_comp();
//	virtual vector<double> get_residuals_vec(const vector<string> &obs_names);
//	virtual void phi_report(ostream &os);
//	virtual bool obs_valid() const;
//	virtual bool numeric_pars_valid() const;
//	virtual bool model_pars_valid() const;
//	virtual bool phi_valid() const;
//	const ModelRun *get_base_model_run_ptr()const;
//	virtual ~ModelRunShallowCopy(){};
//private:
//	const ModelRun *base_model_run_ptr;
//	Parameters numeric_pars;
//	Parameters model_pars;
//	Observations sim_obs;
//	PhiComponets phi_comp;
//	bool numeric_pars_is_valid;
//	bool model_pars_is_valid;
//	bool obs_is_valid;
//	bool phi_is_valid;
//};


#endif /* MODELRUNPP_H_ */
