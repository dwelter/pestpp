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
#ifndef PEST_DATAS_STRUCTS_H_
#define PEST_DATAS_STRUCTS_H_

#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include "Transformable.h"


//class TransformSeq<;
class LaTridiagMatDouble;

class ControlInfo {
public:
	double relparmax;
	double facparmax;
	double facorig;
	double phiredswh;
	int noptmax;
	int jacfile;
	int numcom;
	double phiredstp;
	int nphistp;
	int nphinored;
	double relparstp;
	int nrelpar;
	int noptswitch;
	double splitswh;
	ControlInfo() : relparmax(0.0), facparmax(0.0), facorig(0.0), phiredswh(0.0), noptmax(0),
		phiredstp(0.0), nphistp(0), nphinored(0), relparstp(0.0), nrelpar(0), noptswitch(0),
		splitswh(0.0) {}
};
ostream& operator<< (ostream &os, const ControlInfo& val);


class SVDInfo {
public:
	int maxsing;
	int eigwrite;
	double eigthresh;
	SVDInfo() : maxsing(0), eigwrite(0), eigthresh(1.0e-7) {}
};
ostream& operator<< (ostream &os, const SVDInfo& val);

class ParameterGroupRec {
public:
	string name;
	string inctyp;
	double derinc;
	double derinclb;
	string forcen;
	double derincmul;
	string dermthd;
	double splitthresh;
	double splitreldiff;
	ParameterGroupRec(const string &_name="", const string &_inctyp="", double _derinc=0.0, double _derinclb=0.0, 
		const string &_forcen="", double _derincmul=0.0, 
		const string &_dermthd = "", double _splitthresh = 0.0, double _splitreldiff=.50)
		: name(_name), inctyp(_inctyp), derinc(_derinc), derinclb(_derinclb), forcen(_forcen),
		derincmul(_derincmul), dermthd(_dermthd), splitthresh(_splitthresh), 
		splitreldiff(_splitreldiff){}
	ParameterGroupRec(const ParameterGroupRec &rhs) {*this=rhs;}
	ParameterGroupRec& operator=(const ParameterGroupRec &rhs);
};

ostream& operator<< (ostream &os, const ParameterGroupRec& val);
ostream& operator<< (ostream &os, const map<string, ParameterGroupRec> &val);

class ParameterGroupInfo {
	friend ostream& operator<< (ostream &os, const ParameterGroupInfo &val);
public:
	ParameterGroupInfo() {}
	ParameterGroupInfo(const ParameterGroupInfo&rhs) {*this=rhs;}
	void insert_group(const string &group_name, ParameterGroupRec &rec);

	/** @brief Creates a link from a parameter to its parameters group.
	
	This method add a record to a hash table to link the specified parameter
	to the specified group
	*/
	void insert_parameter_link(const string &parameter_name, const string & group_name);
	const ParameterGroupRec* get_group_rec_ptr(const string &par_name) const;
	const ParameterGroupRec* get_group_by_groupname(const string &group_name) const { return groups.at(group_name); }
	string get_group_name(const string &par_name) const;
	const ParameterGroupInfo& operator=(const ParameterGroupInfo &rhs);
	~ParameterGroupInfo();
private:
	unordered_map<string, ParameterGroupRec*> groups;
	unordered_map<string, ParameterGroupRec*> parameter2group;
	
};


class ParameterRec {
public:
	enum class TRAN_TYPE {NONE, FIXED, TIED, LOG};
	string chglim;
	double lbnd;
	double ubnd;
	double init_value;
	double scale;
	double offset;
	string group;
	bool dercom;
	TRAN_TYPE tranform_type;
	ParameterRec() : chglim(""), lbnd(0.0), ubnd(0.0), init_value(0.0), group(""),
		dercom(false), tranform_type(TRAN_TYPE::NONE), scale(1.0), offset(0.0){}
	bool is_active() const { return !(tranform_type == TRAN_TYPE::FIXED || tranform_type == TRAN_TYPE::TIED); }
};
ostream& operator<< (ostream &os, const ParameterRec& val);

class ParameterInfo {
	friend ostream& operator<< (ostream &os, const ParameterInfo& val);
public:
	Parameters get_low_bnd(const vector<string> &keys) const;
	Parameters get_up_bnd(const vector<string> &keys) const;
	Parameters get_init_value(const vector<string> &keys) const;
	const ParameterRec* get_parameter_rec_ptr(const string &name) const;
	void insert(const string &name, const ParameterRec &rec) {parameter_info[name] = rec;}
	ParameterInfo() {}
	~ParameterInfo() {}	
private:
	unordered_map<string, ParameterRec> parameter_info;
};

ostream& operator<< (ostream &os, const ParameterInfo& val);

class ObservationGroupRec {
public:
	double gtarg;  // optional
	string covfile; // optional
	ObservationGroupRec() : gtarg(0.0), covfile(""){};
	static bool is_regularization(const std::string &grp_name);

};
ostream& operator<< (ostream &os, const ObservationGroupRec& val);

class ObservationRec {
public:
	double weight;
	string group;
	ObservationRec(double _weight=0.0,const string &_group="") : weight(_weight), group(_group) {}
	bool is_regularization() const;
};
ostream& operator<< (ostream &os, const ObservationRec& val);

class ObservationInfo {
public:
	ObservationInfo() {}
	ObservationInfo(const ObservationInfo & rhs) : groups(rhs.groups), observations(rhs.observations){}
	bool is_regularization(const string &obs_name) const;
	unordered_map<string, ObservationGroupRec> groups;
	unordered_map<string, ObservationRec> observations;
	double get_weight(const string &obs_name) const;
	void set_weight(const string &obs_name, double &value);
	string get_group(const string &obs_name) const;
	const ObservationRec* get_observation_rec_ptr(const string &name) const;
	const ObservationGroupRec* get_group_rec_ptr(const string &name) const;
	Observations get_regulatization_obs(const Observations &obs_in);	
};

class ModelExecInfo {
public:
	std::vector<std::string> comline_vec;
	std::vector<std::string> tplfile_vec;
	std::vector<std::string> inpfile_vec;
	std::vector<std::string> insfile_vec;
	std::vector<std::string> outfile_vec;
};



class PestppOptions {
public:
	enum SVD_PACK{EIGEN, PROPACK};
	enum MAT_INV{Q12J, JTQJ};
	PestppOptions(int _n_iter_base = 50, int _n_iter_super=0, int _max_n_super = 50, 
		double _super_eigthres = 1.0E-6, SVD_PACK _svd_pack = PestppOptions::EIGEN,
		MAT_INV _mat_inv = PestppOptions::JTQJ, double _auto_norm = -999,
		double _super_relparmax = 0.1, int max_run_fail=3,
		bool iter_summary_flag = true, bool der_forgive = true);
	void parce_line(const string &line);
	int get_max_n_super() const{return max_n_super;}
	double get_super_eigthres() const{return super_eigthres;}
	int get_n_iter_base() const{return n_iter_base;}
	int get_n_iter_super() const{return n_iter_super;}
	SVD_PACK get_svd_pack() const {return svd_pack;}
	MAT_INV get_mat_inv() const { return mat_inv;}
	double get_auto_norm() const{return auto_norm;}
	double get_super_relparmax() const{ return super_relparmax; }
	int get_max_run_fail() const{ return max_run_fail; }
	int get_max_super_frz_iter()const { return max_super_frz_iter; }
	int get_max_reg_iter()const { return max_reg_iter; }
	const vector<double>& get_base_lambda_vec() const {return base_lambda_vec;}	
	bool get_iter_summary_flag() const { return iter_summary_flag;  }
	bool get_der_forgive() const { return der_forgive; }
	void set_max_n_super(int _max_n_super) {max_n_super = _max_n_super;}
	void set_super_eigthres(double _super_eigthres) {super_eigthres = _super_eigthres;}
	void set_n_iter_base(int _n_iter_base) {n_iter_base = _n_iter_base;}
	void set_n_iter_super(int _n_iter_super){n_iter_super = _n_iter_super;}
	void set_svd_pack(const SVD_PACK _svd_pack=EIGEN) {svd_pack = _svd_pack;}
	void set_mat_inv(const MAT_INV _mat_inv = JTQJ) { mat_inv = _mat_inv; }
	void set_auto_norm(double _auto_norm) {auto_norm = _auto_norm;};
	void set_super_relparmax(double _super_relparmax) { super_relparmax = _super_relparmax; };
	void set_max_run_fail(int _max_run_fail){ max_run_fail = _max_run_fail; }
	void set_max_super_frz_iter(int n) { max_super_frz_iter = n; }
	void set_max_reg_iter(int n) { max_reg_iter = n; }	
	void set_iter_summary_flag(bool _iter_summary_flag){iter_summary_flag = _iter_summary_flag;}
	void set_uncert_flag(bool _flag){uncert = _flag; }
	bool get_uncert_flag()const { return uncert; }
	void set_prediction_names(vector<string> _names){ prediction_names = _names; }
	vector<string> get_prediction_names()const { return prediction_names; }
	void set_parcov_filename(string _filename){ parcov_filename = _filename; }
	string get_parcov_filename()const { return parcov_filename; }
private:
	int n_iter_base;
	int n_iter_super;
	int max_n_super;
	double super_eigthres;
	SVD_PACK svd_pack;
	MAT_INV mat_inv;
	double auto_norm;
	double super_relparmax;
	int max_run_fail;
	int max_super_frz_iter;
	int max_reg_iter;
	vector<double> base_lambda_vec;	
	bool iter_summary_flag;
	bool der_forgive;
	bool uncert;
	vector<string> prediction_names;
	string parcov_filename;
};
ostream& operator<< (ostream &os, const PestppOptions& val);
ostream& operator<< (ostream &os, const ObservationInfo& val);
#endif  /* PEST_DATAS_STRUCTS_H_ */