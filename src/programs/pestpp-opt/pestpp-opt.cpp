// pestpp-opt.cpp : Defines the entry point for the console application.
//
#include <iomanip>
#include "ClpSimplex.hpp"
#include "ClpPresolve.hpp"

void echo_results(ClpSimplex model)
{
	double phi = model.getObjValue();
	const double* decision_var_solution = model.getColSolution();
	const double* constraint_solution = model.getRowPrice();
	const double* obj_coefs_solution = model.getObjCoefficients();
	const std::vector<std::string> decision_var_names = *model.columnNames();
	const std::vector<std::string> constraint_names = *model.rowNames();

	std::cout << "phi: " << phi << std::endl;
	std::cout << std::endl << std::endl << "optimal decision variable values: " << std::endl;
	for (int i = 0; i < decision_var_names.size(); i++)
	{
		std::cout << std::setw(20) << std::left << decision_var_names[i] << ": " << std::setw(20) << std::left << decision_var_solution[i] << std::endl;
	}
	std::cout << "shadow prices on binding constraints:" << std::endl;
	for (int i = 0; i < constraint_names.size(); i++)
	{
		if (constraint_solution[i] == 0.0)
			continue;
		std::cout << std::setw(20) << constraint_names[i] << ": " << std::setw(20) << constraint_solution[i] << std::endl;
	}


	return;
}

int main()
{
	ClpSimplex model;
	int status = 0;
	status = model.readMps("model.mps",true,false);
	if (status)
	{
		std::cout << "error reading mps file" << std::endl;
		return 1;
	}
	ClpPresolve presolve_info;
	ClpSimplex* presolved_model = presolve_info.presolvedModel(model);
	if (presolved_model)
	{
		presolved_model->dual();
		presolve_info.postsolve(true);
	}
	model.checkSolution();
	echo_results(model);
    return 0;
}

