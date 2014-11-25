#include <vector>
#include <string>
#include "Pest.h"
#include "utilities.h"
#include "eigen_tools.h"
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
	if (!pest_utils::check_exist_in(jco_filename))
		throw runtime_error("linear_analysis::linear_analysis() error: jco_filename does not exists");
	
	string pst_filename = jco_filename;
	pst_filename.replace(jco_filename.size() - 4, jco_filename.size() - 1, ".pst");
	if (pest_utils::check_exist_in(pst_filename))
	{
		ifstream ipst(pst_filename);
		Pest pest_scenario;
		pest_scenario.process_ctl_file(ipst,pst_filename);
		ipst.close();

		Mat jco;
		pest_utils::upper_ip(jco_filename);
		if ((jco_filename.find(".JCO") != string::npos) || (jco_filename.find("JCB") != string::npos))
			jco.from_binary(jco_filename);
		else
			jco.from_ascii(jco_filename);

		linear_analysis(jco, pest_scenario);
	}
	else
	{
		pest_utils::upper_ip(jco_filename);
		if ((jco_filename.find(".JCO") != string::npos) || (jco_filename.find("JCB") != string::npos))
			jacobian.from_binary(jco_filename);
		else
			jacobian.from_ascii(jco_filename);
	}
}

linear_analysis::linear_analysis(string &jco_filename, string &parcov_filename, string &obscov_filename)
{
	if (!pest_utils::check_exist_in(jco_filename))
		throw runtime_error("linear_analysis::linear_analysis() error: jco_filename does not exists");
	if (!pest_utils::check_exist_in(parcov_filename))
		throw runtime_error("linear_analysis::linear_analysis() error: parcov_filename does not exists");
	if (!pest_utils::check_exist_in(obscov_filename))
		throw runtime_error("linear_analysis::linear_analysis() error: obscov_filename does not exists");

	pest_utils::upper_ip(jco_filename);
	if ((jco_filename.find(".JCO") != string::npos) || (jco_filename.find("JCB") != string::npos))
		jacobian.from_binary(jco_filename);
	else
		jacobian.from_ascii(jco_filename);

	pest_utils::upper_ip(parcov_filename);
	if (parcov_filename.find(".UNC") != string::npos)
		parcov.from_uncertainty_file(parcov_filename);
	else if (parcov_filename.find(".PST") != string::npos)
		parcov.from_parameter_bounds(parcov_filename);
	else
		parcov.from_ascii(parcov_filename);

	pest_utils::upper_ip(obscov_filename);
	if (obscov_filename.find(".UNC") != string::npos)
		obscov.from_uncertainty_file(obscov_filename);
	else if (obscov_filename.find(".PST") != string::npos)
		obscov.from_observation_weights(obscov_filename);
	else
		obscov.from_ascii(obscov_filename);
}

linear_analysis::linear_analysis(Mat _jacobian, Pest pest_scenario)
{
	jacobian = _jacobian;
	parcov.from_parameter_bounds(pest_scenario);
	obscov.from_observation_weights(pest_scenario);

}

//linear_analysis linear_analysis::get(vector<string> &new_par_names, vector<string> &new_obs_names)
//{
//
//}

void linear_analysis::align()
{
	vector<string> common;
	if (jacobian.get_col_names() != parcov.get_col_names())
	{
		//common = get_common(jacobian.get_col_names(), parcov.get_row_names());
		//parcov = parcov.get(common);
		//jacobian = jacobian.get(jacobian.get_row_names(), common);
		parcov = parcov.get(jacobian.get_col_names());
	}
	if (jacobian.get_row_names() != obscov.get_col_names())
	{
		//common = get_common(jacobian.get_row_names(), obscov.get_row_names());
		//obscov = obscov.get(common);
		//jacobian = jacobian.get(jacobian.get_row_names(), common);
		obscov = obscov.get(jacobian.get_row_names());
	}
	for (int i = 0; i < predictions.size(); i++)
	{
		if (jacobian.get_col_names() != predictions[i].get_row_names())
		{
			predictions[i] = predictions[i].get(jacobian.get_col_names(), predictions[i].get_row_names());
		}
	}
}

double linear_analysis::prior_parameter_variance(string &par_name)
{
	int ipar = find(parcov.rn_ptr()->begin(), parcov.rn_ptr()->end(), par_name) - parcov.rn_ptr()->begin();
	if (ipar == parcov.nrow())
		throw runtime_error("lienar_analysis::prior_parameter_variance() error: parameter: " + par_name + " not found");
	return parcov.get_matrix().diagonal()[ipar];	
}

double linear_analysis::posterior_parameter_variance(string &par_name)
{
	if (posterior.nrow() == 0) calc_posterior();
	int ipar = find(posterior.rn_ptr()->begin(), posterior.rn_ptr()->end(), par_name) - posterior.rn_ptr()->begin();
	if (ipar == posterior.nrow())
		throw runtime_error("lienar_analysis::posterior_parameter_variance() error: parameter: " + par_name + " not found");
	return posterior.get_matrix().diagonal()[ipar];
}

Mat linear_analysis::posterior_parameter_matrix()
{
	if (posterior.nrow() == 0) calc_posterior();
	return posterior;
}

void linear_analysis::calc_posterior()
{
	align();
	posterior = Covariance(*parcov.rn_ptr(), ((*jacobian.transpose().eptr() * *obscov.inv().eptr() * *jacobian.eptr()) + *parcov.inv().eptr())).inv();
}



void linear_analysis::set_predictions(vector<string> preds)
{
	const vector<string>* obs_names = jacobian.rn_ptr();
	for (auto pred : preds)
	{
		pest_utils::upper_ip(pred);
		if (find(obs_names->begin(), obs_names->end(), pred) != obs_names->end())
		{
			Mat mpred = jacobian.extract(pred, vector<string>());
			mpred.transpose_ip();
			predictions.push_back(mpred);
		}
		else
		{
			if (!pest_utils::check_exist_in(pred))
			{
				throw runtime_error("linear_analysis::set_predictions() error: pred: " + pred + " not found in jco rows and is not an accessible file");
			}
			Mat mpred;
			mpred.from_ascii(pred);
			if (mpred.ncol() != 1)
			{
				if (mpred.nrow() == 1)
				{
					mpred.transpose_ip();
				}
				else
				{
					throw runtime_error("linear_analysis::set_predictions() error: pred: " + pred + "must be shape (1,npar)");
				}
			}
			//if the pred vector is not completely aligned with the JCO
			if (*mpred.rn_ptr() != *jacobian.cn_ptr())
			{
				vector<string> missing;
				const vector<string> *mpred_par_names = mpred.rn_ptr();
				for (auto jco_par : *jacobian.cn_ptr())
				{
					if (find(mpred_par_names->begin(), mpred_par_names->end(), jco_par) == mpred_par_names->end())
					{
						missing.push_back(jco_par);
					}
				}
				if (missing.size() > 0)
				{
					stringstream ss;
					for (auto m : missing) ss << m << ",";
					throw runtime_error("linear_analysis::set_predictions() error: parameters missing from pred: " + pred + " : " + ss.str());
				}
				mpred = mpred.get(*jacobian.cn_ptr(), *mpred.cn_ptr());
			}
			predictions.push_back(mpred);

		}
	}

}

void linear_analysis::set_predictions(vector<Mat> preds)
{

}

