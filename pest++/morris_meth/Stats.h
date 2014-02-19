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

class RunningStats
{
    //Based on paper "Computing the standard deviation efficiently" by Mark Hoemmen, Berkeley
public:
	RunningStats() : n(0), mk(0.0), mk_abs(0.0), qk(0.0) {}
    void reset();
    void add(double sample);
    void add(const std::vector<double> &sample);
    double comp_var();
    double comp_sigma();
    double comp_mean();
	double comp_abs_mean();
	int comp_nsamples();
	~RunningStats() {};
private:
    long n;
    double mk;
	double mk_abs;
    double qk;
};