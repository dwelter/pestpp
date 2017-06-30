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

class PhiStats
{
public:
	PhiStats(const map<string, double> &_phi_info);
	PhiStats() { ; }
	void rec_report(ofstream &f_rec);
	void csv_report(ofstream &csv, const int iter, const vector<string> &names);
	void update(const map<string, double> &_phi_info);
	const double get_mean() const { return mean; }
	const double get_std() const { return std;  }
	static void initialize_csv(ofstream &csv, const vector<string> &names);

private:
	map<string, double> phi_map;
	double mean, std, min, max;
	int size;
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
	Pest &pest_scenario;
	FileManager &file_manager;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;

	int iter,subset_size;
	double last_best_lam, last_best_mean,last_best_std;

	vector<double> lam_mults;

	string fphi_name;
	ofstream fphi;
	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;

	ParameterEnsemble pe, pe_base;
	ObservationEnsemble oe, oe_base;
	//Eigen::MatrixXd prior_pe_diff;
	Eigen::MatrixXd Am;
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> obscov_inv_sqrt;


	//EnsemblePair run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	vector<int> run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	vector<ObservationEnsemble> run_lambda_ensembles(vector<ParameterEnsemble> &pe_lams);
	//map<string, double> get_phi_vec_stats(map<string,PhiComponets> &phi_info);
	//map<string,PhiComponets> get_phi_info(ObservationEnsemble &_oe);
	PhiStats report_and_save();
	void lam_test_report(double lambda, PhiStats &phistats);
	map<string, double> get_phi_map(ObservationEnsemble &_oe);
};


#endif 