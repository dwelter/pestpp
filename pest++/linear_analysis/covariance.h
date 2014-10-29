
#include <string>
#include <vector>
#include<Eigen/Sparse>

#include "Pest.h"

using namespace std;

class Mat
{
	
public:
	enum class MatType{ DENSE, DIAGONAL, BLOCK };
	Mat(){};
	Mat(vector<string> _row_names, vector<string> _col_names)
		{ row_names = _row_names,col_names = _col_names; }
	Mat(vector<string> _row_names, vector<string> _col_names, Eigen::SparseMatrix<double> _matrix);
	Mat(vector<string> _row_names, vector<string> _col_names, Eigen::SparseMatrix<double> _matrix, MatType _mattype);
	
	vector<string> get_row_names(){ return row_names; }
	vector<string> get_col_names(){ return col_names; }
	Eigen::SparseMatrix<double> get_matrix(){ return matrix; }
	MatType get_mattype(){ return mattype; }

	void to_ascii(const string &filename);
	void from_ascii(const string &filename);
	void to_binary(const string &filename);
	void from_binary(const string &filename);
	

	void align(vector<string> &other_row_names, vector<string> &other_col_names);
	void transpose();
	
	
	Mat get(vector<string> &other_row_names, vector<string> &other_col_names);
	Mat get(){ return Mat(row_names, col_names, matrix); }
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
	int icode = 2;
	MatType mattype;
	vector<string> read_namelist(ifstream &in, int &nitems);
};


class Covariance : public Mat
{
public:
	Covariance(vector<string> &names);
	Covariance();
	Covariance(Mat _mat);
	Covariance(vector<string> _row_names, Eigen::SparseMatrix<double> _matrix);
	Covariance get(vector<string> &other_names);
	void from_uncertainty_file(const string &filename);
	void from_parameter_bounds(Pest &pest_scenario);
	void from_observation_weights(Pest &pest_scenario);

	void from_parameter_bounds(const string &pst_filename);
	void from_observation_weights(const string &pst_filename);

	void to_uncertainty_file(const string &filename);
};
