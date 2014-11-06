#include <vector>
#include <string>
#include "Pest.h"
#include "covariance.h"
#include "linear_analysis.h"


vector<string> get_common(vector<string> v1, vector<string> v2)
{
	vector<string> common;
	for (auto e1 : v1)
		if (find(v2.begin(), v2.end(), e1) != v2.end()) common.push_back(e1);
	return common;
}



linear_analysis::linear_analysis(string &jco_filename)
{

}

linear_analysis::linear_analysis(string &jco_filename, string &parcov_filename, string &obscov_filename)
{

}

linear_analysis::linear_analysis(Mat _jacobian, Pest pest_scenario)
{

}

linear_analysis linear_analysis::get(vector<string> &new_par_names, vector<string> &new_obs_names)
{

}

vector<string> linear_analysis::check_and_align()
{
	vector<string> warnings, errors,common_pars,common_obs;

	//TODO: align everything with the jco...use it as the base
	if (jacobian.get_col_names() != parcov.get_col_names())
		common_pars = get_common(jacobian.get_col_names(), parcov.get_col_names());
	if (jacobian.get_row_names() != obscov.get_col_names())
		common_obs = get_common(jacobian.get_row_names(), obscov.get_col_names());





	

	errors.insert(errors.end(),warnings.begin(), warnings.end());
	return errors;
}