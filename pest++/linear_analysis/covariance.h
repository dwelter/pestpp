
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
	enum class MatType{ DENSE, DIAGONAL, BLOCK };
	Mat(){ autoalign = true; };
	Mat(string filename);
	Mat(vector<string> _row_names, vector<string> _col_names)
		{ row_names = _row_names,col_names = _col_names; }
	Mat(vector<string> _row_names, vector<string> _col_names, 
		Eigen::SparseMatrix<double> _matrix,bool _autoalign=true);
	Mat(vector<string> _row_names, vector<string> _col_names, 
		Eigen::SparseMatrix<double> _matrix, MatType _mattype,bool _autoalign=true);
	
	vector<string> get_row_names(){ return row_names; }
	vector<string> get_col_names(){ return col_names; }
	Eigen::SparseMatrix<double> get_matrix(){ return matrix; }
	
	const Eigen::SparseMatrix<double>* get_matrix_ptr();
	const Eigen::SparseMatrix<double>* get_U_ptr();
	const Eigen::SparseMatrix<double>* get_V_ptr();
	const Eigen::VectorXd* get_s_ptr();
	bool get_autoalign(){ return autoalign; }
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
	void SVD();

	Mat get(vector<string> &other_row_names, vector<string> &other_col_names);
	Mat get(){ return Mat(row_names, col_names, matrix); }
	//Mat extract(vector<string> &extract_row_names, vector<string> &extract_col_names);
	Mat extract(vector<string> &extract_row_names, vector<string> &extract_col_names);
	void drop_rows(vector<string> &drop_row_names);
	void drop_cols(vector<string> &drop_col_names);

	int nrow(){ return matrix.rows(); }
	int ncol(){ return matrix.cols(); }	

	bool isAligned(Mat other_mat);

	Mat operator *(Mat &other_mat);
	Mat operator *(double val);
	Mat operator +(Mat &other_mat);
	Mat operator -(Mat &other_mat);
	


protected:
	bool autoalign;
	Eigen::SparseMatrix<double> matrix;
	Eigen::SparseMatrix<double> U;
	Eigen::SparseMatrix<double> V;
	Eigen::VectorXd s;
	Eigen::SparseMatrix<double> lower_chol;
	vector<string> row_names;
	vector<string> col_names;
	int icode = 2;
	MatType mattype = MatType::DENSE;
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