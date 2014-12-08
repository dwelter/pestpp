#include <vector>
#include <string>
#include <sstream>
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
		parcov = parcov.get(jacobian.get_col_names());
	}
	if (jacobian.get_row_names() != obscov.get_col_names())
	{
		obscov = obscov.get(jacobian.get_row_names());
	}
	for (auto &p : predictions)
	{
		if (*jacobian.cn_ptr() != *p.second.rn_ptr())
		{
			Mat new_pred = p.second.get(*jacobian.cn_ptr(),*p.second.rn_ptr());
			predictions[p.first] = new_pred;
		}
	}
}

map<string, double> linear_analysis::worth(vector<string> &obs_names)
{
	if (predictions.size() == 0)
		throw runtime_error("linear_analysis::posterior_predictive_worth() error: no predictions are set");

	align();
	for (int i = 0; i != obs_names.size(); i++)
		pest_utils::upper_ip(obs_names[i]);

	//check the inputs
	vector<string> errors;
	const vector<string>* jobs_names = jacobian.rn_ptr();
	for (auto &oname : obs_names)
	{
		if (find(jobs_names->begin(), jobs_names->end(), oname) == jobs_names->end())
			errors.push_back("obs not found in jacobian: "+oname);
	}
	if (errors.size() > 0)
	{
		stringstream ss;
		for (auto &e : errors)
			ss << e << ',';
		throw runtime_error("linear_analysis::posterior_predictive_worth() errors: " + ss.str());
	}


	

	//get a list of remaining obs
	vector<string> keep_obs_names;
	for (auto &oname : *jobs_names)
	{
		if (find(obs_names.begin(), obs_names.end(), oname) == obs_names.end())
		{
			keep_obs_names.push_back(oname);
		}
	}
	//get a new linear analysis object with only the keep names
	Covariance* parcov_ptr = &parcov;
	map<string,Mat>* pred_ptr = &predictions;
	linear_analysis keep(jacobian.get(keep_obs_names, *jacobian.cn_ptr()), *parcov_ptr, obscov.get(keep_obs_names),*pred_ptr);
	
	//calculate posterior with all obs and with only keep obs
	map<string, double> org_post = posterior_prediction_variance();
	map<string, double> keep_post = keep.posterior_prediction_variance();
	
	map<string, double> results;
	for (auto &pred : predictions)
		results[pred.first] = keep_post[pred.first] - org_post[pred.first];

	return results;
}


map<string,pair<double, double>> linear_analysis::contribution(vector<string> &cond_par_names)
{
	if (predictions.size() == 0)
		throw runtime_error("linear_analysis::predictive_contribution() error: no predictions are set");

	align();
	for (int i = 0; i != cond_par_names.size(); i++)
		pest_utils::upper_ip(cond_par_names[i]);
	//check the inputs
	vector<string> errors;
	const vector<string>* jpar_names = jacobian.cn_ptr();
	for (auto par_name : cond_par_names)
	{
		if (find(jpar_names->begin(), jpar_names->end(), par_name) == jpar_names->end())
			errors.push_back("par not found in jacobian: " + par_name);
	}
	if (errors.size() > 0)
	{
		stringstream ss;
		for (auto e : errors)
			ss << e << ',';
		throw runtime_error("linear_analysis::predictive_contribution() errors: " + ss.str());
	}

	//the parameters that will remain uncertain
	vector<string> keep_par_names;
	for (auto pname : *parcov.rn_ptr())
	{
		if (find(cond_par_names.begin(), cond_par_names.end(), pname) == cond_par_names.end())
			keep_par_names.push_back(pname);
	}

	map<string,Mat> cond_preds;
	for (auto &pred : predictions)
		cond_preds[pred.first] = pred.second.get(keep_par_names, *pred.second.cn_ptr());

	//get a new linear analysis object - we can use the same obscov b/c it doesn't change
	Covariance* obscov_ptr = &obscov;
	linear_analysis cond_la(jacobian.get(*jacobian.rn_ptr(), keep_par_names), condition_on(keep_par_names, cond_par_names), *obscov_ptr,cond_preds);

	//calculate the prior and posterior with all parameters
	map<string, double> org_prior = prior_prediction_variance();
	map<string, double> org_post = posterior_prediction_variance();

	//calculate the prior and posterior with some conditioned parameters
	map<string, double> cond_prior = cond_la.prior_prediction_variance();
	map<string, double> cond_post = cond_la.posterior_prediction_variance();

	//calculate the reduction in prior and posterior uncertainty for each prediction
	map<string, pair<double, double>> results;
	for (auto pred : predictions)
	{
		pair<double, double> reduction(org_prior[pred.first] - cond_prior[pred.first], org_post[pred.first] - cond_post[pred.first]);
		results[pred.first] = reduction;
	}
	return results;
}


Covariance linear_analysis::condition_on(vector<string> &keep_par_names, vector<string> &cond_par_names)
{
	//C11
	Covariance new_parcov = parcov.get(keep_par_names);
	//if parcov is diagonal, then there is no conditioning
	if (parcov.get_mattype() == Mat::MatType::DIAGONAL)
	{
		return new_parcov;
	}		
	//the portion of parcov that is becoming "known": C22
	Covariance cond_parcov = parcov.get(cond_par_names);
	//C22^-1
	cond_parcov.inv_ip();
	//C12
	Mat upper_off_diag = parcov.get(keep_par_names, cond_par_names);
	//C11 - (C12*C22^-1*C21)
	Eigen::SparseMatrix<double> mcond = *new_parcov.eptr() - (*upper_off_diag.eptr() * *cond_parcov.eptr() * *upper_off_diag.transpose().eptr());
	Covariance new_cond_parcov(keep_par_names, mcond);
	return new_cond_parcov;
}

double linear_analysis::prior_parameter_variance(string &par_name)
{
	int ipar = find(parcov.rn_ptr()->begin(), parcov.rn_ptr()->end(), par_name) - parcov.rn_ptr()->begin();
	if (ipar == parcov.nrow())
		throw runtime_error("lienar_analysis::prior_parameter_variance() error: parameter: " + par_name + " not found");
	const Eigen::SparseMatrix<double>* ptr = parcov.eptr();
	return ptr->diagonal()[ipar];	
}

double linear_analysis::posterior_parameter_variance(string &par_name)
{
	if (posterior.nrow() == 0) calc_posterior();
	int ipar = find(posterior.rn_ptr()->begin(), posterior.rn_ptr()->end(), par_name) - posterior.rn_ptr()->begin();
	if (ipar == posterior.nrow())
		throw runtime_error("lienar_analysis::posterior_parameter_variance() error: parameter: " + par_name + " not found");
	const Eigen::SparseMatrix<double>* ptr = posterior.eptr();
	return ptr->diagonal()[ipar];
}


Mat linear_analysis::posterior_parameter_matrix()
{
	if (posterior.nrow() == 0) calc_posterior();
	return posterior;
}

Mat* linear_analysis::posterior_parameter_ptr()
{
	if (posterior.nrow() == 0) calc_posterior();
	Mat* ptr = &posterior;
	return ptr;
}

double linear_analysis::prior_prediction_variance(string &pred_name)
{
	pest_utils::upper_ip(pred_name);
	map<string, Mat>::iterator p_iter = predictions.find(pred_name);
	if (p_iter == predictions.end())
		throw runtime_error("linear_analysis::prior_pred_variance() error: pred:" + pred_name + " not found in predicitons");
	Eigen::SparseMatrix<double> result = (*p_iter->second.transpose().eptr() * *parcov.eptr() * *p_iter->second.eptr());
	return result.valuePtr()[0];

}

map<string, double> linear_analysis::prior_prediction_variance()
{
	map<string, double> result;
	for (auto &pred : predictions)
	{
		string pname(pred.first);
		result[pname] = prior_prediction_variance(pname);
	}
	return result;
}

double linear_analysis::posterior_prediction_variance(string &pred_name)
{
	pest_utils::upper_ip(pred_name);
	map<string, Mat>::iterator p_iter = predictions.find(pred_name);
	if (p_iter == predictions.end())
		throw runtime_error("linear_analysis::prior_pred_variance() error: pred:" + pred_name + " not found in predicitons");
	if (posterior.nrow() == 0) calc_posterior();
	Eigen::SparseMatrix<double> result = (*p_iter->second.transpose().eptr() * *posterior.eptr() * *p_iter->second.eptr());
	return result.valuePtr()[0];

}

map<string, double> linear_analysis::posterior_prediction_variance()
{
	map<string, double> result;
	for (auto &pred : predictions)
	{
		string pname(pred.first);
		result[pname] = posterior_prediction_variance(pname);
	}
	return result;
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
			if (predictions.find(pred) != predictions.end())
				throw runtime_error("linear_analysis::set_predictions() error: pred:" + pred + " already in predictions");
			Mat mpred = jacobian.extract(pred, vector<string>());
			mpred.transpose_ip();
			predictions[pred] = mpred;
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
			string pname = mpred.get_col_names()[0];
			if (predictions.find(pname) != predictions.end())
				throw runtime_error("linear_analysis::set_predictions() error: pred:" + pred + " already in predictions");
			predictions[pname] = mpred;

		}
	}

}

void linear_analysis::set_predictions(vector<Mat> preds)
{
	throw runtime_error("not implemented");
}

