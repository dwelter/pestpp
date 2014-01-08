#include <vector>
#include <map>
#include <string>

double vec_mean(const std::vector<double> &data_vec);
double vec_mean_missing_data(const std::vector<double> &data_vec, double missing_val);

double vec_var(const std::vector<double> &data_vec);

std::map<std::string, double>  vec_covar_missing_data(const std::vector<double> &x_vec, const std::vector<double> &y_vec, double missing_val);

std::map<std::string, double>  vec_calc_stats(const std::vector<double> &data_vec);
std::map<std::string, double>  vec_calc_stats_missing_data(const std::vector<double> &data_vec, double missing_val);

double sobol_u_missing_data(const std::vector<double> &x_vec, const std::vector<double> &y_vec, double missing_val);
