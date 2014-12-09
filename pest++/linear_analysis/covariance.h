#ifndef QSQRT_MATRIX_H_
#define QSQRT_MATRIX_H_
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include<Eigen/Sparse>

#include "Pest.h"
#include "logger.h"

using namespace std;



class Mat
{
public:
	enum class MatType{ DIAGONAL, SPARSE, DENSE };
	Mat(){ autoalign = true; };
	Mat(string filename);
	
	Mat(vector<string> _row_names, vector<string> _col_names,
		Eigen::DiagonalMatrix<double,Eigen::Dynamic> _matrix, bool _autoalign = true);

	Mat(vector<string> _row_names, vector<string> _col_names, 
		Eigen::SparseMatrix<double> _matrix,bool _autoalign=true);
	
	vector<string> get_row_names(){ return row_names; }
	vector<string> get_col_names(){ return col_names; }
	const vector<string>* rn_ptr();
	const vector<string>* cn_ptr();
	Eigen::SparseMatrix<double> get_matrix(){ return matrix; }
	
	const Eigen::SparseMatrix<double>* eptr();
	const Eigen::SparseMatrix<double>* get_U_ptr();
	const Eigen::SparseMatrix<double>* get_V_ptr();
	const Eigen::VectorXd* get_s_ptr();

	Mat get_U();
	Mat get_V();
	Mat get_s();

	MatType get_mattype(){ return mattype; }

	void to_ascii(const string &filename);
	void from_ascii(const string &filename);
	void to_binary(const string &filename);
	void from_binary(const string &filename);

	void transpose_ip();
	Mat transpose();
	Mat T();
	Mat inv();
	void inv_ip();
	void SVD();

	Mat identity();

	Mat get(const vector<string> &other_row_names, const vector<string> &other_col_names);
	Mat leftCols(const int idx);
	Mat rightCols(const int idx);
	Mat extract(const vector<string> &extract_row_names, const vector<string> &extract_col_names);
	Mat extract(const string &extract_row_name, const vector<string> &extract_col_names);
	Mat extract(const vector<string> &extract_row_names, const string &extract_col_name);
	void drop_rows(const vector<string> &drop_row_names);
	void drop_cols(const vector<string> &drop_col_names);

	int nrow(){ return row_names.size(); }
	int ncol(){ return col_names.size(); }	


protected:
	bool autoalign;
	Eigen::SparseMatrix<double> matrix;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> diagonal;
	Eigen::SparseMatrix<double> U;
	Eigen::SparseMatrix<double> V;
	Eigen::VectorXd s;
	Eigen::SparseMatrix<double> lower_chol;
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
	Covariance(string filename);
	Covariance(Mat _mat);
	Covariance(vector<string> _row_names, Eigen::SparseMatrix<double> _matrix);
	
	Covariance get(vector<string> &other_names);
	Mat get(vector<string> &other_row_names, vector<string> &other_col_names){ return Mat::get(other_row_names, other_col_names); }
	void drop(vector<string> &drop_names);
	Covariance extract(vector<string> &extract_names);

	void from_uncertainty_file(const string &filename);
	void from_parameter_bounds(Pest &pest_scenario);
	void from_observation_weights(Pest &pest_scenario);

	void from_parameter_bounds(const string &pst_filename);
	void from_observation_weights(const string &pst_filename);

	void to_uncertainty_file(const string &filename);

	vector<Eigen::VectorXd> draw(int ndraws);
	vector<double> standard_normal(default_random_engine gen);
	void cholesky();

private:
	Eigen::SparseMatrix<double> lower_cholesky;
};

ostream& operator<< (std::ostream &os, Mat mat);
#endif