#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Stats.h"


using namespace std;





double vec_mean(const vector<double> &data_vec)
{
	auto n = data_vec.size();

	double tmp_mean = 0;
	for (const auto &x : data_vec)
	{
		tmp_mean += x / n;
	}
	return tmp_mean;
}

double vec_mean_missing_data(const vector<double> &data_vec, double missing_val)
{
	auto n_est = data_vec.size();
	size_t  n_act = 0;

	double tmp_mean = 0;
	for (const auto &x : data_vec)
	{
		if (x != missing_val)
		{
			tmp_mean += x / n_est;
			++n_act;
		}
	}
	tmp_mean = tmp_mean * n_est / n_act;
	return tmp_mean;
}


double vec_var(const vector<double> &data_vec)
{
	double mk = 0;
	double qk = 0;
	double xk;
	int n = data_vec.size();

	mk = data_vec[0];
	qk = 0;
	for (int k=2; k<=n; ++k)
	{
		xk =  data_vec[k-1];
		qk += (k - 1) * pow((xk - mk), 2.0) / k;// must happend before mk is updated
		mk += (xk - mk) / k;
	}
	return qk / (n-1);
}

map<string, double>  vec_covar(const vector<double> &x_vec, const vector<double> &y_vec)
{
	map<string, double> ret_val;

	double m_x_k = 0;
	double m_y_k = 0;
	double m_xy_k = 0;
	double x_k, y_k, xy_k; 

	assert (x_vec.size() == y_vec.size());
	int n = x_vec.size();

	m_x_k = x_vec[0];
	m_y_k = y_vec[0];
	m_xy_k = m_x_k * m_y_k;

	double q_x_k = 0;
	double q_y_k = 0;
	double co_q_xy_k = 0;

	for (int k=2; k<=n; ++k)
	{
		x_k =  x_vec[k-1];
		y_k =  y_vec[k-1];
		xy_k =  x_k * y_k;

		// must happend before m_x_k and m_y_k are updated
		q_x_k += (k - 1) * pow((x_k - m_x_k), 2.0) / k;
		q_y_k += (k - 1) * pow((y_k - m_y_k), 2.0) / k;

		m_x_k += (x_k - m_x_k) / k;
		m_y_k += (y_k - m_y_k) / k;
		m_xy_k += (xy_k - m_xy_k) / k;
	}

	co_q_xy_k = m_xy_k - m_x_k * m_y_k;
	ret_val["mean_1"] = m_x_k;
	ret_val["mean_2"] = m_y_k;
	ret_val["var_1"] = q_x_k / (n-1);
	ret_val["var_2"] = q_y_k / (n-1);
	ret_val["covar"] = co_q_xy_k / (n-1);

	return ret_val;
}


map<string, double>  vec_covar_missing_data(const vector<double> &x_vec, const vector<double> &y_vec, double missing_val)
{
	map<string, double> ret_val;

	double m_x_k = 0;
	double m_y_k = 0;
	double m_xy_k = 0;
	double x_k, y_k, xy_k; 
	size_t n_actual = 0;

	assert (x_vec.size() == y_vec.size());
	int n = x_vec.size();
	
	int k;
	for (k=1; k<=n; ++k)
	{
		if (x_vec[k-1] != missing_val && x_vec[k-1] != missing_val)
		{
			m_x_k = x_vec[k-1];
			m_y_k = y_vec[k-1];
			m_xy_k = x_vec[0] * y_vec[0];
			++n_actual;
			break;
		}
	}

	double q_x_k = 0;
	double q_y_k = 0;
	double co_q_xy_k = 0;
	
	for (k+=1; k<=n; ++k)
	{
		x_k =  x_vec[k-1];
		y_k =  y_vec[k-1];

		if (x_k != missing_val && x_k != missing_val)
		{
			xy_k =  x_k * y_k;

			// must happend before m_x_k and m_y_k are updated
			q_x_k += (k - 1) * pow((x_k - m_x_k), 2.0) / k;
			q_y_k += (k - 1) * pow((y_k - m_y_k), 2.0) / k;

			m_x_k += (x_k - m_x_k) / k;
			m_y_k += (y_k - m_y_k) / k;
			m_xy_k += (xy_k - m_xy_k) / k;

			++n_actual;
		}
	}
	co_q_xy_k = m_xy_k - m_x_k * m_y_k;
	ret_val["mean_1"] = m_x_k;
	ret_val["mean_2"] = m_y_k;
	ret_val["var_1"] = q_x_k / (n_actual - 1);
	ret_val["var_2"] = q_y_k / (n_actual - 1);
	ret_val["covar"] = co_q_xy_k / (n_actual - 1);

	return ret_val;
}

map<string, double>  vec_calc_stats(const vector<double> &data_vec)
{
	map<string, double> ret_val;

	double mk = 0;
	double mk_abs = 0;
	double qk = 0;
	double xk;
	int n = data_vec.size();

	mk = data_vec[0];
	mk_abs = abs(data_vec[0]);
	qk = 0;
	for (int k=2; k<=n; ++k)
	{
		xk =  data_vec[k-1];
		qk += (k - 1) * pow((xk - mk), 2.0) / k; // must happend before mk is updated
		mk += (xk - mk) / k;
		mk_abs += (abs(xk) - mk_abs) / k;
	}
	ret_val["mean"] = mk;
	ret_val["mean_abs"] = mk_abs;
	ret_val["var"] =  qk/ (n-1);
	ret_val["std_dev"] = sqrt(qk/(n-1));
	return ret_val;
}

map<string, double>  vec_calc_stats_missing_data(const vector<double> &data_vec, double missing_val)
{
	map<string, double> ret_val;

	double mk = 0;
	double mk_prev = 0;
	double mk_abs = 0;
	double qk = 0;
	double xk;
	int n = data_vec.size();
	size_t n_actual = 0;
	
	int k;
	for (k=1; k<=n; ++k)
	{
		if (data_vec[k-1] != missing_val)
		{
			mk = data_vec[k-1];
			mk_abs = abs(mk);
			qk = 0;
			++n_actual;
			break;
		}
	}

	for (k+=1; k<=n; ++k)
	{
		xk =  data_vec[k-1];
		if (xk != missing_val)
		{
			mk_prev = mk;
			mk += (xk - mk_prev) / k;
			mk_abs += (abs(xk) - mk_abs) / k;
			qk += (k - 1) * pow((xk - mk_prev), 2.0) / k;
			++n_actual;
		}
	}
	ret_val["mean"] = mk;
	ret_val["mean_abs"] = mk_abs;
	ret_val["var"] =  qk / (n_actual - 1);
	ret_val["std_dev"] = sqrt(qk / (n_actual - 1));
	return ret_val;
}



double sobol_u_missing_data(const std::vector<double> &x_vec, const std::vector<double> &y_vec, double missing_val)
{
	double m_xy_k = 0;
	double x_k, y_k, xy_k; 

	assert (x_vec.size() == y_vec.size());
	size_t n = x_vec.size();
	size_t n_actual = 0;
	
	int k;
	for (k=1; k<=n; ++k)
	{
		x_k =  x_vec[k-1];
		y_k =  y_vec[k-1];
		if (x_k != missing_val && y_k != missing_val)
		{
			xy_k =  x_k * y_k;
			m_xy_k = xy_k;
			++n_actual;
			break;
		}
	}
	
	for (k+=1; k<=n; ++k)
	{
		x_k =  x_vec[k-1];
		y_k =  y_vec[k-1];

		if (x_k != missing_val && x_k != missing_val)
		{
			xy_k =  x_k * y_k;
			m_xy_k += (xy_k - m_xy_k) / k;
			++n_actual;
		}
	}
	return m_xy_k * n / (n_actual-1);
}

void RunningStats::reset()
{
    n = 0;
    mk = 0.0;
	mk_abs = 0.0;
    qk = 0.0;
}

void RunningStats::add(double sample)
{
    ++n;
    if (n==1)
	{
        mk = sample;
		mk_abs = abs(sample);
        qk = 0.0;
	}
	else
	{
		qk += (n - 1) * pow((sample - mk), 2.0) / n; // must happend before mk is updated
		mk += (sample - mk) / n;
		mk_abs += (abs(sample) - mk_abs) / n;
	}
}

void RunningStats::add(const std::vector<double> &sample)
{
    for (auto &s : sample)
	{
       add(s); 
	}
}

double RunningStats::comp_var()
{
    return qk / (n-1);
}

double RunningStats::comp_sigma()
{
    return sqrt(comp_var());
}

double RunningStats::comp_mean()
{
return mk;
}

double RunningStats::comp_abs_mean()
{
	return mk_abs;
}

int RunningStats::comp_nsamples()
{
	return n;
}