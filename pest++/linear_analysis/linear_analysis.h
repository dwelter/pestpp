#ifndef LINEAR_ANALYSIS_H_
#define LINEAR_ANALYSIS_H_

//#include <vector>
//#include <string>

#include "Pest.h"
#include "covariance.h"

using namespace std;
class linear_analysis
{
public:
	//the easiest constructor, builds parcov and obscov from the pst associated with jco_filename
	linear_analysis(string &jco_filename);

	//loads parcov and obscov from files, can .pst, .mat or .unc files
	linear_analysis(string &jco_filename, string &parcov_filename, string &obscov_filename);

	//load parcov and obscov from parameter bounds and observation weights
	linear_analysis(Mat _jacobian, Pest pest_scenario);

	//directly from Mat objects
	linear_analysis(Mat _jacobian, Mat _parcov, Mat _obscov, map<string,Mat> _predictions)
	{
		jacobian = _jacobian; parcov = _parcov; obscov = _obscov,predictions=_predictions;
	}

	void set_predictions(vector<string> preds);
	void set_predictions(vector<Mat> preds);

	//get a new linear analysis object consisting of a subset of par and obs names
	linear_analysis get(vector<string> &new_par_names, vector<string> &new_obs_names);

	//exposed schur functionality
	//from the diagonal of parcov
	double prior_parameter_variance(string &par_name);
	//a map of <par_name,variance> from diagonal of parcov
	map<string, double> prior_parameter_variance();
	Mat prior_parameter_matrix() { return parcov; }

	//from the diagonal of schur's complement
	double posterior_parameter_variance(string &par_name);
	//a map of <par_name,variance> from the diagonal of schur's complement
	map<string, double> posterior_parameter_variance();
	//the full matrix
	Mat posterior_parameter_matrix();
	Mat* posterior_parameter_ptr();
	

	//prior predictive variance from parcov
	double prior_prediction_variance(string &pred_name);
	//map <pred_name,variance> from parcov
	map<string, double> prior_prediction_variance();

	//posterior predictive variance from schur's complement
	double posterior_prediction_variance(string &pred_name);
	//map <pred_name,variance> from schur's complement
	map<string, double> posterior_prediction_variance();

	//the reduction in predictive variance from some obs
	double posterior_predictive_worth(string &pred_name, vector<string> &obs_names);
	//<pred_name,variance_reduction> from some obs
	map<string, double> posterior_predictive_worth(vector<string> &obs_names);

	//reduction in prior and posterior predictive variance from perfect knowledge of some pars
	//pair<double, double> predictive_contribution(string &pred_name, vector<string> &par_names);
	//<pred_name,prior and posterior variance_reduction> from perfect knowledge of some pars
	map<string, pair<double, double>> predictive_contribution(vector<string> &par_names);




	//exposed error variance functionality

	//<err_var_component("null","solution","omitted"),error_variance> for a set of singular values and a parameter
	map<string, vector<double>> parameter_error_variance_components(vector<int> sing_vals, string &par_name);
	//<err_var_component("null","solution","omitted"),error_variance> for a set of singular values and a prediction
	map<string, vector<double>> prediction_error_variance_components(vector<int> sing_vals, string &pred_name);

	//extract elements from the jacobian, parcov, and predictions and set them as omitted
	void extract_omitted(vector<string> &omitted_par_names);
	void extract_omitted(string &omitted_par_name){ extract_omitted(vector<string>{omitted_par_name}); }


	//other stuff

	//perform the svd and set the components
	void SVD();

	//<par_names,ident>
	map<string, double> parameter_ident(int sing_val);

	//<singular_value,vector<sup_obs>>
	map<int, vector<double>> super_obs(vector<int> sing_vals);
	map<int, vector<double>> super_obs(vector<int> sing_vals, vector<string> &sup_obs_names);

	//<singular_value,vector<sup_par>>
	map<int, vector<double>> super_par(vector<int> sing_vals);
	map<int, vector<double>> super_par(vector<int> sing_vals, vector<string> &sup_par_names);

	Mat get_jacobian(){ return jacobian; }
	Mat get_omitted_jacobian(){ return omitted_jacobian; }
	Covariance get_parcov(){ return parcov; }
	Covariance get_omitted_parcov(){ return omitted_parcov; }
	Covariance get_obscov(){ return obscov; }

	Mat get_U(){ return U; }
	Mat get_S(){ return S; }
	Mat get_V(){ return V; }

	Mat* get_jacobian_ptr(){ return &jacobian; }
	Covariance* get_parcov_ptr(){ return &parcov; }
	Covariance* get_obscov_ptr(){ return &obscov; }
	Mat* get_omitted_jacobian_ptr(){ return &omitted_jacobian; }
	Covariance* get_omitted_parcov_ptr(){ return &omitted_parcov; }

	Mat* get_U_ptr(){ return &U; }
	Mat* get_S_ptr(){ return &S; }
	Mat* get_V_ptr(){ return &V; }

	map<string,Mat> get_predictions(){ return predictions; }
	vector<Mat> get_omitted_predictions(){ return omitted_predictions; }
	
	//returns a list of warnings and errors, aligns the different linear components, sets isaligned to true
	void align();
	bool get_isaligned(){ return isaligned; }

private:
	Mat jacobian;
	Covariance parcov;
	Covariance obscov;
	//vector<Mat> predictions;
	map<string, Mat> predictions;
	bool isaligned;

	Covariance posterior;
	void calc_posterior();
	Covariance condition_on(vector<string> &keep_par_names,vector<string> &cond_par_names);

	//svd stuff
	Mat U,S,V;
	
	//error variance stuff
	Mat R, G, I_minus_R;
	vector<Mat> omitted_predictions;
	Mat omitted_jacobian;
	Covariance omitted_parcov;
	//scale the jacobian by parcov
	void kl_scale();








};
#endif