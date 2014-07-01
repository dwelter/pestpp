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

void RestartController::process_rst_file(std::ifstream &fin)
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
			restart_option = RestartOption::RESUME_NEW_ITERATION;
		}
		else if (tokens[0] == "termination_info_1")
		{
			//tokens[1] is noptmax;
			convert_ip(tokens[2], nopt_count);
			//tokens[3] is nphinored);
			convert_ip(tokens[4], nphinored_count);
			//tokens[5] is term_ctl.nrelpar);

		}
		else if (tokens[0] == "termination_info_2")
		{
			convert_ip(tokens[1], nrelpar_count);
			//tokens[2] is nphistp;
			//tokens[3] is phiredstp;
			//tokens[4] is relparstp;

		}
		else if (tokens[0] == "termination_info_3")
		{
			for (int i = 1; i<tokens.size(); ++i)
			{
				double val = convert_cp<double>(tokens[1]);
				lowest_phi.push_back(val);
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
		//else if (tokens[0] == "upgrade_model_runs_built")
		//{
		//	restart_option = RestartOption::RESUME_UPGRADE_RUNS;
		//}
	}
}

void RestartController::update_termination_ctl(TerminationController &term_ctl)
{
	term_ctl.nopt_count = nopt_count;
	term_ctl.nphinored_count = nphinored_count;
	term_ctl.nrelpar_count = nrelpar_count;
	term_ctl.lowest_phi = lowest_phi;
}

RestartController::~RestartController(void)
{
}
