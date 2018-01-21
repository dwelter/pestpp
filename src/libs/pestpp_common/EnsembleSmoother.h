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
		       Covariance *_parcov_inv, double *_reg_factor);
	void update(ObservationEnsemble &oe, ParameterEnsemble &pe);
	double get_mean(phiType pt);
	double get_std(phiType pt);
	double get_max(phiType pt);
	double get_min(phiType pt);
	map<string, double>* get_phi_map(PhiHandler::phiType &pt);
	void report();
	void write(int iter_num, int total_runs);
	vector<int> get_idxs_greater_than(double bad_phi, ObservationEnsemble &oe);

private:
	map<string, double> get_summary_stats(phiType pt);
	string get_summary_string(phiType pt);
	string PhiHandler::get_summary_header();
	void prepare_csv(ofstream &csv,vector<string> &names);
	map<string, double> calc_meas(ObservationEnsemble &oe);
	map<string, double> calc_regul(ParameterEnsemble &pe);
	map<string, double> calc_actual(ObservationEnsemble &oe);
	map<string, double> calc_composite(map<string,double> &_meas, map<string,double> &_regul);
	//map<string, double>* get_phi_map(PhiHandler::phiType &pt);
	void write_csv(int iter_num, int total_runs,ofstream &csv, phiType pt,
		           vector<string> &names);

	double *reg_factor;
	vector<string> oreal_names,preal_names;
	Pest* pest_scenario;
	FileManager* file_manager;
	ObservationEnsemble* oe_base;
	ParameterEnsemble* pe_base;
	Covariance* parcov_inv;
	map<string, double> meas;
	map<string, double> regul;
	map<string, double> composite;
	map<string, double> actual;


};


class IterEnsembleSmoother
{
public:
	IterEnsembleSmoother(Pest &_pest_scenario, FileManager &_file_manager,
		OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log,
		RunManagerAbstract* _run_mgr_ptr);
	void initialize();
	void solve();
	void finalize();
	void throw_ies_error(string &message);

private:
	int  verbose_level;
	Pest &pest_scenario;
	FileManager &file_manager;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;
	PhiHandler ph;
	Covariance parcov_inv;
	double reg_factor;

	int iter,subset_size;
	bool use_subset;
	double last_best_lam, last_best_mean,last_best_std;
	double lambda_max, lambda_min;
	vector<double> lam_mults;

	//string fphi_name;
	//ofstream fphi;
	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;

	ParameterEnsemble pe, pe_base;
	ObservationEnsemble oe, oe_base;
	//Eigen::MatrixXd prior_pe_diff;
	Eigen::MatrixXd Am;
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> obscov_inv_sqrt;


	//EnsemblePair run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	vector<int> run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, vector<int> &real_idxs=vector<int>());
	vector<ObservationEnsemble> run_lambda_ensembles(vector<ParameterEnsemble> &pe_lams, vector<double> &lam_vals);
	//map<string, double> get_phi_vec_stats(map<string,PhiComponets> &phi_info);
	//map<string,PhiComponets> get_phi_info(ObservationEnsemble &_oe);
	void report_and_save();
	void save_mat(string prefix, Eigen::MatrixXd &mat);
	void initialize_pe(Covariance &cov);
	void initialize_oe(Covariance &cov);
	void drop_bad_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	//void check_ensembles(ObservationEnsemble &oe, ParameterEnsemble &pe);
};

#endif 