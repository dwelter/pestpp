
#include <string>
#include <vector>
#include<Eigen/Sparse>

#include "Pest.h"

using namespace std;

class Mat
{
	
public:
	Mat(){};
	Mat(vector<string> _row_names, vector<string> _col_names);
	Mat(vector<string> _row_names, vector<string> _col_names,Eigen::SparseMatrix<double> _matrix);
	enum class MatType{ DENSE, DIAGONAL, BANDED };
	
	void to_ascii(const string &filename);
	void from_ascii(const string &filename);
	void to_binary(const string &filename);
	void from_binary(const string &filename);
	
	void align(vector<string> &other_row_names, vector<string> &other_col_names);
	
	Mat get(vector<string> &other_row_names, vector<string> &other_col_names);
	Mat extract(vector<string> &other_row_names, vector<string> &other_col_names);
	void drop(vector<string> &other_row_names, vector<string> &other_col_names);

	int nrow(){ return matrix.rows(); }
	int ncol(){ return matrix.cols(); }	

	bool isAligned(Mat other_mat);

	Mat operator *(Mat &other_mat);
	Mat operator +(Mat &other_mat);
	Mat operator -(Mat &other_mat);

protected:
	Eigen::SparseMatrix<double> matrix;
	vector<string> row_names;
	vector<string> col_names;
	int icode = 999;
	MatType mattype;
	vector<string> read_namelist(ifstream &in, int &nitems);
};


class Covariance : public Mat
{
public:
	Covariance(vector<string> &names);
	Covariance();

	void from_uncertainty_file(const string &filename);
	void from_parameter_bounds(Pest &pest_scenario);
	void from_observation_weights(Pest &pest_scenario);

	void from_parameter_bounds(const string &pst_filename);
	void from_observation_weights(const string &pst_filename);

	void to_uncertainty_file(const string &filename);
};
