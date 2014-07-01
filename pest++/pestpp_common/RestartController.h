#ifndef RESTART_CTL_H_
#define RESTART_CTL_H_

#include <string>
#include <fstream>
#include <iostream>

class TerminationController;

class RestartController
{
public:
	enum class RestartOption {NONE, REUSE_JACOBIAN, RESUME_NEW_ITERATION, RESUME_JACOBIAN_RUNS, RESUME_UPGRADE_RUNS};
	enum class IterationType{BASE, SUPER};
	RestartController(void);
	RestartOption get_restart_option() const { return restart_option; }
	RestartOption& get_restart_option() { return restart_option; }
	void process_rst_file(std::ifstream &fin);
	void update_termination_ctl(TerminationController &term_ctl);
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
};

#endif /* RESTART_CTL_H_ */