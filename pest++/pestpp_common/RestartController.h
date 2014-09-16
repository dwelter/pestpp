#ifndef RESTART_CTL_H_
#define RESTART_CTL_H_

#include <string>
#include <fstream>
#include <iostream>

class TerminationController;
class SVDSolver;

class RestartController
{
public:
	enum class RestartOption {NONE, REUSE_JACOBIAN, RESUME_NEW_ITERATION, RESUME_JACOBIAN_RUNS};
	enum class IterationType{BASE, SUPER};
	RestartController(void);
	static void write_start_iteration(std::ostream &fout, const SVDSolver &svd_solver, int _iter_num, int _global_iter_num);
	static void write_start_parameters_updated(std::ostream &fout, const std::string &parameter_filename);
	static void write_finish_parameters_updated(std::ostream &fout, const std::string &parameter_filename);
	static void write_jac_runs_built(std::ostream &fout); 
	static void RestartController::write_iteration_complete(std::ostream &fout);
	RestartOption get_restart_option() const { return restart_option; }
	RestartOption& get_restart_option() { return restart_option; }
	void process_rst_file(std::ifstream &fin);
	void update_termination_ctl(TerminationController &term_ctl);
	//void update_best_par_file(const std::string &_best_par_file) { best_par_file = _best_par_file; }
	~RestartController(void);
private:
	int global_iter_no;
	int local_iter_no;
	RestartOption restart_option; 
	IterationType iteration_type;
	int nopt_count;
	int nphinored_count;
	int nrelpar_count;
	std::vector<double> lowest_phi;
	std::string best_par_file;
};

#endif /* RESTART_CTL_H_ */