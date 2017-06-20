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

	int iter;
	double last_best_lam;
	string fphi_name;
	ofstream fphi;
	vector<string> phi_stat_names;
	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;

	ParameterEnsemble pe;
	ObservationEnsemble oe, oe_base;
	Eigen::MatrixXd prior_pe_diff;
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> obscov_inv_sqrt;


	EnsemblePair run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	map<string, double> get_phi_vec_stats();
	map<string,PhiComponets> get_phi_info();
	void report_and_save();

};


#endif 