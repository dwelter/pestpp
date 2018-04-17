#ifndef ENSEMBLESMOOTHER_H_
#define ENSEMBLESMOOTHER_H_

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
#include "ObjectiveFunc.h"
 
class PhiHandler
{
public:

	enum phiType { MEAS, COMPOSITE, REGUL, ACTUAL };
	PhiHandler() { ; }
	PhiHandler(Pest *_pest_scenario, FileManager *_file_manager, 
		       ObservationEnsemble *_oe_base, ParameterEnsemble *_pe_base,
		       Covariance *_parcov, double *_reg_factor);
	void update(ObservationEnsemble &oe, ParameterEnsemble &pe);
	double get_mean(phiType pt);
	double get_std(phiType pt);
	double get_max(phiType pt);
	double get_min(phiType pt);
	map<string, double>* get_phi_map(PhiHandler::phiType &pt);
	void report();
	void write(int iter_num, int total_runs);
	vector<int> get_idxs_greater_than(double bad_phi, ObservationEnsemble &oe);

	Eigen::MatrixXd get_obs_resid(ObservationEnsemble &oe);
	Eigen::MatrixXd get_par_resid(ParameterEnsemble &pe);
	Eigen::MatrixXd get_actual_obs_resid(ObservationEnsemble &oe);
	Eigen::VectorXd get_q_vector();
	vector<string> get_lt_obs_names() { return lt_obs_names; }
	vector<string> get_gt_obs_names() { return gt_obs_names; }

	void apply_ineq_constraints(Eigen::MatrixXd &resid);

private:
	map<string, double> get_summary_stats(phiType pt);
	string get_summary_string(phiType pt);
	string get_summary_header();
	void prepare_csv(ofstream &csv,vector<string> &names);
	void prepare_group_csv(ofstream &csv, vector<string> extra = vector<string>());

	map<string, Eigen::VectorXd> calc_meas(ObservationEnsemble &oe);
	map<string, Eigen::VectorXd> calc_regul(ParameterEnsemble &pe);
	map<string, Eigen::VectorXd> calc_actual(ObservationEnsemble &oe);
	map<string, double> calc_composite(map<string,double> &_meas, map<string,double> &_regul);
	//map<string, double>* get_phi_map(PhiHandler::phiType &pt);
	void write_csv(int iter_num, int total_runs,ofstream &csv, phiType pt,
		           vector<string> &names);
	void write_group_csv(int iter_num, int total_runs, ofstream &csv, 
		map<string, double> extra = map<string, double>());

	double *reg_factor;
	vector<string> oreal_names,preal_names;
	Pest* pest_scenario;
	FileManager* file_manager;
	ObservationEnsemble* oe_base;
	ParameterEnsemble* pe_base;
	//Covariance parcov_inv;
	Eigen::VectorXd parcov_inv_diag;
	map<string, double> meas;
	map<string, double> regul;
	map<string, double> composite;
	map<string, double> actual;

	vector<string> lt_obs_names;
	vector<string> gt_obs_names;

	map<string, vector<int>> obs_group_idx_map;
	map<string, vector<int>> par_group_idx_map;
	
	map<string, double> get_obs_group_contrib(Eigen::VectorXd &phi_vec);
	map<string, double> get_par_group_contrib(Eigen::VectorXd &phi_vec);	

};


class IterEnsembleSmoother
{
public:
	IterEnsembleSmoother(Pest &_pest_scenario, FileManager &_file_manager,
		OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log,
		RunManagerAbstract* _run_mgr_ptr);
	void initialize();
	void iterate_2_solution();
	void pareto_iterate_2_solution();
	void finalize();
	void throw_ies_error(string message);


private:
	int  verbose_level;
	Pest &pest_scenario;
	FileManager &file_manager;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;
	PhiHandler ph;
	Covariance parcov;
	double reg_factor;

	int iter,subset_size;
	bool use_subset;
	double last_best_lam, last_best_mean,last_best_std;
	double lambda_max, lambda_min;
	int warn_min_reals, error_min_reals;
	vector<double> lam_mults;

	//string fphi_name;
	//ofstream fphi;
	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;

	ParameterEnsemble pe, pe_base;
	ObservationEnsemble oe, oe_base;
	//Eigen::MatrixXd prior_pe_diff;
	Eigen::MatrixXd Am;
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> obscov_inv_sqrt, parcov_inv_sqrt;

	void solve();
	void adjust_pareto_weight(string &obsgroup, double wfac);

	//EnsemblePair run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	vector<int> run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs=vector<int>());
	vector<ObservationEnsemble> run_lambda_ensembles(vector<ParameterEnsemble> &pe_lams, vector<double> &lam_vals);
	//map<string, double> get_phi_vec_stats(map<string,PhiComponets> &phi_info);
	//map<string,PhiComponets> get_phi_info(ObservationEnsemble &_oe);
	void report_and_save();
	void save_mat(string prefix, Eigen::MatrixXd &mat);
	bool initialize_pe(Covariance &cov);
	bool initialize_oe(Covariance &cov);
	void drop_bad_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	//void check_ensembles(ObservationEnsemble &oe, ParameterEnsemble &pe);
	template<typename T, typename A>
	void message(int level, string &_message, vector<T, A> _extras);
	void message(int level, string &_message);
	
	template<typename T, typename A>
	void message(int level, const char* _message, vector<T, A> _extras);// { message(level, string(_message), _extras); }
	void message(int level, const char* _message);// { message(level, string(_message)); }

	template<typename T>
	void message(int level, string &_message, T extra);

	template<typename T>
	void message(int level, const char* _message, T extra);

	void sanity_checks();

	void add_bases();
};

#endif 