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
	ControlInfo() : relparmax(0.0), facparmax(0.0), facorig(0.0), phiredswh(0.0), noptmax(0),
		phiredstp(0.0), nphistp(0), nphinored(0), relparstp(0.0), nrelpar(0){}
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
	ParameterGroupRec(const string &_name="", const string &_inctyp="", double _derinc=0.0, double _derinclb=0.0, 
		const string &_forcen="", double _derincmul=0.0, 
		const string &_dermthd="") : name(_name), inctyp(_inctyp), derinc(_derinc), derinclb(_derinclb), forcen(_forcen),
			derincmul(_derincmul), dermthd(_dermthd) {}
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
	string group;
	bool dercom;
	TRAN_TYPE tranform_type;
	ParameterRec() : chglim(""), lbnd(0.0), ubnd(0.0), init_value(0.0), group(""),
		dercom(false), tranform_type(TRAN_TYPE::NONE){}
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
	PestppOptions(int _n_iter_base = 50, int n_iter_super=0, int _max_n_super = 50, double _super_eigthres = 1.0E-6, SVD_PACK _svd_pack = PestppOptions::EIGEN,
		MAT_INV _mat_inv = PestppOptions::JTQJ, double _auto_norm = -999, double _super_relparmax = 0.1);
	void parce_line(const string &line);
	int get_max_n_super() const{return max_n_super;}
	double get_super_eigthres() const{return super_eigthres;}
	int get_n_iter_base() const{return n_iter_base;}
	int get_n_iter_super() const{return n_iter_super;}
	SVD_PACK get_svd_pack() const {return svd_pack;}
	MAT_INV get_mat_inv() const { return mat_inv;}
	double get_auto_norm() const{return auto_norm;}
	double get_super_relparmax() const{ return super_relparmax; }

	void set_max_n_super(int _max_n_super) {max_n_super = _max_n_super;}
	void set_super_eigthres(double _super_eigthres) {super_eigthres = _super_eigthres;}
	void set_n_iter_base(int _n_iter_base) {n_iter_base = _n_iter_base;}
	void set_n_iter_super(int _n_iter_super){n_iter_super = _n_iter_super;}
	void set_svd_pack(const SVD_PACK _svd_pack=EIGEN) {svd_pack = _svd_pack;}
	void set_mat_inv(const MAT_INV _mat_inv = JTQJ) { mat_inv = _mat_inv; }
	void set_auto_norm(double _auto_norm) {auto_norm = _auto_norm;};
	void set_super_relparmax(double _super_relparmax) { super_relparmax = _super_relparmax; };

private:
	int max_n_super;;
	double super_eigthres;
	int n_iter_base;
	int n_iter_super;
	SVD_PACK svd_pack;
	MAT_INV mat_inv;
	double auto_norm;
	double super_relparmax;
};

ostream& operator<< (ostream &os, const ObservationInfo& val);
#endif  /* PEST_DATAS_STRUCTS_H_ */