#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "RestartController.h"
#include "TerminationController.h"
#include "utilities.h"

using namespace std;
using namespace::pest_utils;

RestartController::RestartController(void)
{
	restart_option = RestartOption::NONE;
}

void RestartController::process_rst_file(std::ifstream &fin, TerminationController &term_ctl)
{
	string line;
	vector<string> tokens;
	while (getline(fin, line))
	{
		tokens.clear();
		tokenize(line, tokens);

		if (tokens.empty())
		{
		}
		else if (tokens[0] == "start_iteration")
		{
			convert_ip(tokens[1], local_iter_no);
			convert_ip(tokens[1], global_iter_no);
		}
		else if (tokens[0] == "termination_info_1")
		{
			convert_ip(tokens[1], term_ctl.noptmax);
			convert_ip(tokens[2], term_ctl.nopt_count);
			convert_ip(tokens[3], term_ctl.nphinored);
			convert_ip(tokens[4], term_ctl.nphinored_count);
			convert_ip(tokens[5], term_ctl.nrelpar);

		}
		else if (tokens[0] == "termination_info_2")
		{
			convert_ip(tokens[1], term_ctl.nrelpar);
			convert_ip(tokens[2], term_ctl.nphistp);
			convert_ip(tokens[3], term_ctl.phiredstp);
			convert_ip(tokens[4], term_ctl.relparstp);

		}
		else if (tokens[0] == "termination_info_3")
		{
			for (int i = 1; i<tokens.size(); ++i)
			{
				double val = convert_cp<double>(tokens[1]);
				term_ctl.lowest_phi.push_back(val);
			}
		}
		else if (tokens[0] == "native_par_iteration")
		{
			iteration_type = IterationType::BASE;
		}
		else if (tokens[0] == "super_par_iteration")
		{
			iteration_type = IterationType::SUPER;
		}
		else if (tokens[0] == "jacobian_model_runs_built")
		{
			restart_option = RestartOption::RESUME_JACOBIAN_RUNS;
		}
		else if (tokens[0] == "upgrade_model_runs_built")
		{
			restart_option = RestartOption::RESUME_UPGRADE_RUNS;
		}
	}
}


RestartController::~RestartController(void)
{
}
