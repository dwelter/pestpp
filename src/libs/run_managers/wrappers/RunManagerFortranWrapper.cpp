#include "RunManagerFortranWrapper.h"
#include "utilities.h"
#include "RunManagerPanther.h"
#include "RunManagerSerial.h"
#include "pest_error.h"

typedef class RunManagerAbstract RunManagerAbstract;

static RunManagerAbstract *_run_manager_ptr_ = nullptr;
static string _run_manager_error;
static ofstream fout_run_manager_log_file;

using namespace pest_utils;

extern "C"
{

int rmif_create_serial_(char *f_comline, int  *comline_str_len, int *comline_array_len,
	char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
	char *f_inp, int  *inp_str_len, int *inp_array_len,
	char *f_ins, int  *ins_str_len, int *ins_array_len,
	char *f_out, int  *out_str_len, int *out_array_len,
	char *f_storfile, int *storfile_len,
	char *f_rundir, int *rundir_len, int *n_max_fail)

{
	_run_manager_error.clear();
	vector<string> comline_vec =  fortran_str_array_2_vec(f_comline, *comline_str_len, *comline_array_len);
	vector<string> tpl_vec =  fortran_str_array_2_vec(f_tpl, *tpl_str_len, *tpl_array_len);
	vector<string> inp_vec =  fortran_str_array_2_vec(f_inp, *inp_str_len, *inp_array_len);
	vector<string> ins_vec =  fortran_str_array_2_vec(f_ins, *ins_str_len, *ins_array_len);
	vector<string> out_vec =  fortran_str_array_2_vec(f_out, *out_str_len, *out_array_len);
	string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
	string rundir =  fortran_str_2_string(f_rundir, *rundir_len);
	int err = 0;
	try {
		_run_manager_ptr_ = new RunManagerSerial(comline_vec, tpl_vec, inp_vec, ins_vec,
			out_vec, storfile, rundir, *n_max_fail);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch(char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}


int rmif_create_panther_(
	char *f_storfile, int *storfile_len,
	char *f_port, int *f_port_len,
	char *f_info_filename, int *info_filename_len, int *n_max_fail,
	double *overdue_reched_fac, double *overdue_giveup_fac)
{
	_run_manager_error.clear();
	int err = 0;
	try {
		string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
		string port =  fortran_str_2_string(f_port, *f_port_len);
		string info_filename =  fortran_str_2_string(f_info_filename, *info_filename_len);
		fout_run_manager_log_file.open(info_filename);
		_run_manager_ptr_ = new RunManagerPanther(storfile, 
			port, fout_run_manager_log_file, *n_max_fail, *overdue_reched_fac, *overdue_giveup_fac);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}


int rmif_err_msg_(char *fortran_str, int *str_len)
{
	int err = 0;
	string_to_fortran_char(_run_manager_error, fortran_str, *str_len);
	return err;
}

int rmif_add_run_(double *parameter_data, int *npar, int *model_exe_index, int *id)
{
	_run_manager_error.clear();
	int err = 0;
	try {
		vector<double> data(parameter_data, parameter_data+*npar);
		*id = _run_manager_ptr_->add_run(data, *model_exe_index);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_add_run_with_info_(double *parameter_data, int *npar, int *model_exe_index, int *id,
	char *f_info_txt, int  *info_txt_len, double *info_value)
{
	_run_manager_error.clear();
	int err = 0;
	try {
		string info_txt = fortran_str_2_string(f_info_txt, *info_txt_len);
		vector<double> data(parameter_data, parameter_data + *npar);
		*id = _run_manager_ptr_->add_run(data, *model_exe_index, info_txt, *info_value);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (...)
	{
		err = 1;
	}
	return err;
}

int rmif_initialize_(char *f_pname, int  *pname_str_len, int *pname_array_len,
				 char *f_oname, int  *oname_str_len, int *oname_array_len)

{
	_run_manager_error.clear();
	int err = 0;
	try {
		vector<string> pname_vec =  fortran_str_array_2_vec(f_pname, *pname_str_len, *pname_array_len);
		vector<string> oname_vec =  fortran_str_array_2_vec(f_oname, *oname_str_len, *oname_array_len);
		_run_manager_ptr_->initialize(pname_vec, oname_vec);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_initialize_restart_(char *f_storfile, int *storfile_len)
{
	_run_manager_error.clear();
	int err = 0;
	try {
		string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
		_run_manager_ptr_->initialize_restart(storfile);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_reinitialize_()
{
	_run_manager_error.clear();
    int err = 0;
	try {
	    _run_manager_ptr_->reinitialize();
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_run_()
{
	_run_manager_error.clear();
	int err = 0;
	try {
		_run_manager_ptr_->run();
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_run_until_(int *condition, int *no_ops, double *time_sec, int *return_cond)
{
	_run_manager_error.clear();
	int err = 0;
	RunManagerAbstract::RUN_UNTIL_COND enum_input_cond;
	RunManagerAbstract::RUN_UNTIL_COND enum_return_cond;
	enum_input_cond = static_cast<RunManagerAbstract::RUN_UNTIL_COND>(*condition);
	try {
		enum_return_cond = _run_manager_ptr_->run_until(enum_input_cond, *no_ops, *time_sec);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (...)
	{
		err = 1;
	}
	*return_cond = static_cast<int>(enum_return_cond);
	return err;
}

int rmif_cancel_run_(int *run_id)
{
	_run_manager_error.clear();
	int err = 0;
	try {
		_run_manager_ptr_->cancel_run(*run_id);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (...)
	{
		err = 1;
	}
	return err;
}

int rmif_get_run_(int *run_id, double *parameter_data, int *npar, double *obs_data, int *nobs)
{
	_run_manager_error.clear();
	int err = 1;
	bool success = false;
	size_t n_par = *npar;
	size_t n_obs = *nobs;
	try {
		success = _run_manager_ptr_->get_run(*run_id, parameter_data, n_par, obs_data, n_obs);
		if (success) err = 0;
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch(...) {
		err = 1;
	}
	return err;
}


int rmif_get_run_with_info_(int *run_id, double *parameter_data, int *npar, double *obs_data, int *nobs,
	char *f_info_txt, int  *info_txt_len, double *info_value)
{
	_run_manager_error.clear();
	int err = 1;
	bool success = false;
	size_t n_par = *npar;
	size_t n_obs = *nobs;
	try {
		string info_txt;
		success = _run_manager_ptr_->get_run(*run_id, parameter_data, n_par, obs_data, n_obs, info_txt, *info_value);
		string_to_fortran_char(info_txt, f_info_txt, *info_txt_len);
		if (success) err = 0;
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (...) {
		err = 1;
	}
	return err;
}


int rmif_get_run_info_(int *run_id, int *run_status, int *model_exe_index, char *f_info_txt, int  *info_txt_len, double *info_value)
{
	_run_manager_error.clear();
	int err = 1;
	try {
		string info_txt;
		_run_manager_ptr_->get_info(*run_id, *run_status, *model_exe_index, info_txt, *info_value);
		err = 0;
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (...) {
		err = 1;
	}
	return err;
}

int rmif_get_run_status_info_(int *run_id, int *run_status, double *max_runtime, int* n_concurrent_runs)
{
	_run_manager_error.clear();
	int err = 1;
	try {
		_run_manager_ptr_->get_run_status_info(*run_id, *run_status, *max_runtime, *n_concurrent_runs);
		err = 0;
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (...) {
		err = 1;
	}
	return err;
}

int rmif_delete_()
{
	_run_manager_error.clear();
	int err = 0;
	try {
		delete _run_manager_ptr_;
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
		catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_get_num_failed_runs_(int *nfail)
{
	_run_manager_error.clear();
    int err = 0;
	*nfail = -999;
	try
	{
        const std::set<int> &fail_set = _run_manager_ptr_->get_failed_run_ids();
        *nfail = fail_set.size();
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
    catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_get_failed_run_ids_(int *run_id_array, int *len_run_id_array)
{
	_run_manager_error.clear();
    int err = 0;
    try
	{
		std::fill_n(run_id_array, *len_run_id_array, -999);
        const std::set<int> &fail_set = _run_manager_ptr_->get_failed_run_ids();
        size_t n_failed_tmp = min(fail_set.size(), size_t(*len_run_id_array));
        std::copy_n(fail_set.begin(), n_failed_tmp, run_id_array);
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
    catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_get_num_total_runs_(int *nruns)
{
	_run_manager_error.clear();
    int err = 0;
    try {
        *nruns = _run_manager_ptr_->get_total_runs();
	}
	catch (const exception &e)
	{
		err = 1;
		_run_manager_error = e.what();
	}
	catch (char const *e)
	{
		err = 1;
		_run_manager_error = e;
	}
	catch (const string e)
	{
		err = 1;
		_run_manager_error = e;
	}
    catch(...)
	{
		err = 1;
	}
	return err;
}

}
