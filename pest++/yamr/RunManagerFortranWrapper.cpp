#include "RunManagerFortranWrapper.h"
#include "utilities.h"
#include "RunManagerYAMR.h"
#include "RunManagerSerial.h"
#include "RunManagerGenie.h"
#include "pest_error.h"

typedef struct RunManagerAbstract RunManagerAbstract;

static RunManagerAbstract *_run_manager_ptr_ = nullptr;
static ofstream fout_run_manager_log_file;

using namespace pest_utils;

extern "C"
{

int RMIF_CREATE_SERIAL(char *f_comline, int  *comline_str_len, int *comline_array_len,
	char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
	char *f_inp, int  *inp_str_len, int *inp_array_len,
	char *f_ins, int  *ins_str_len, int *ins_array_len,
	char *f_out, int  *out_str_len, int *out_array_len,
	char *f_storfile, int *storfile_len,
	char *f_rundir, int *rundir_len, int *n_max_fail)

{
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
	catch(...)
	{
		err = 1;
	}
	return err;
}


int RMIF_CREATE_YAMR(char *f_comline, int  *comline_str_len, int *comline_array_len,
	char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
	char *f_inp, int  *inp_str_len, int *inp_array_len,
	char *f_ins, int  *ins_str_len, int *ins_array_len,
	char *f_out, int  *out_str_len, int *out_array_len,
	char *f_storfile, int *storfile_len,
	char *f_port, int *f_port_len, 
	char *f_info_filename, int *info_filename_len, int *n_max_fail)

{
	int err = 0;
	try {
		vector<string> comline_vec =  fortran_str_array_2_vec(f_comline, *comline_str_len, *comline_array_len);
		vector<string> tpl_vec =  fortran_str_array_2_vec(f_tpl, *tpl_str_len, *tpl_array_len);
		vector<string> inp_vec =  fortran_str_array_2_vec(f_inp, *inp_str_len, *inp_array_len);
		vector<string> ins_vec =  fortran_str_array_2_vec(f_ins, *ins_str_len, *ins_array_len);
		vector<string> out_vec =  fortran_str_array_2_vec(f_out, *out_str_len, *out_array_len);
		string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
		string port =  fortran_str_2_string(f_port, *f_port_len);
		string info_filename =  fortran_str_2_string(f_port, *info_filename_len);
		fout_run_manager_log_file.open(info_filename);
		_run_manager_ptr_ = new RunManagerYAMR(comline_vec, tpl_vec, inp_vec, ins_vec, out_vec, storfile, 
			port, fout_run_manager_log_file, *n_max_fail);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int RMIF_CREATE_GENIE(char *f_comline, int  *comline_str_len, int *comline_array_len,
	char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
	char *f_inp, int  *inp_str_len, int *inp_array_len,
	char *f_ins, int  *ins_str_len, int *ins_array_len,
	char *f_out, int  *out_str_len, int *out_array_len,
	char *f_storfile, int *storfile_len,
	char *f_host, int *f_host_len,
	char *f_genie_tag, int *genie_tag_len)
{
	int err = 0;
	try {
		vector<string> comline_vec =  fortran_str_array_2_vec(f_comline, *comline_str_len, *comline_array_len);
		vector<string> tpl_vec =  fortran_str_array_2_vec(f_tpl, *tpl_str_len, *tpl_array_len);
		vector<string> inp_vec =  fortran_str_array_2_vec(f_inp, *inp_str_len, *inp_array_len);
		vector<string> ins_vec =  fortran_str_array_2_vec(f_ins, *ins_str_len, *ins_array_len);
		vector<string> out_vec =  fortran_str_array_2_vec(f_out, *out_str_len, *out_array_len);
		string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
		string host =  fortran_str_2_string(f_host, *f_host_len);
		string genie_tag =  fortran_str_2_string(f_genie_tag, *genie_tag_len);
		_run_manager_ptr_ = new RunManagerGenie(comline_vec, tpl_vec, inp_vec, ins_vec, out_vec, storfile, genie_tag);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}


int RMIF_ADD_RUN(double *parameter_data, int *npar, int *id)
{
	int err = 0;
	try {
		vector<double> data(parameter_data, parameter_data+*npar);
		*id = _run_manager_ptr_->add_run(data);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int RMIF_INITIALIZE(char *f_pname, int  *pname_str_len, int *pname_array_len,
				 char *f_oname, int  *oname_str_len, int *oname_array_len)

{
	int err = 0;
	try {
		vector<string> pname_vec =  fortran_str_array_2_vec(f_pname, *pname_str_len, *pname_array_len);
		vector<string> oname_vec =  fortran_str_array_2_vec(f_oname, *oname_str_len, *oname_array_len);
		_run_manager_ptr_->initialize(pname_vec, oname_vec);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int RMIF_INITIALIZE_RESTART(char *f_storfile, int *storfile_len)
{
	int err = 0;
	try {
		string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
		_run_manager_ptr_->initialize_restart(storfile);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int RMIF_REINITIALIZE()
{
    int err = 0;
	try {
	    _run_manager_ptr_->reinitialize();
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int RMIF_RUN()
{
	int err = 0;
	try {
		_run_manager_ptr_->run();
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int RMIF_GET_RUN(int *run_id, double *parameter_data, int *npar, double *obs_data, int *nobs)
{
	int err = 1;
	bool success = false;
	size_t n_par = *npar;
	size_t n_obs = *nobs;
	try {
		success = _run_manager_ptr_->get_run(*run_id, parameter_data, n_par, obs_data, n_obs);
		if (success) err = 0;
	}
	catch(PestIndexError ex) {
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int RMFI_DELETE()
{
	int err = 0;
	try {
		delete _run_manager_ptr_;
	}
		catch(...)
	{
		err = 1;
	}
	return err;
}

int RMIF_GET_NUM_FAILED_RUNS(int *nfail)
{
    int err = 0;
	*nfail = -999;
	try
	{
        const std::set<int> &fail_set = _run_manager_ptr_->get_failed_run_ids();
        *nfail = fail_set.size();
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int RMIF_GET_FAILED_RUN_IDS(int *run_id_array, int *len_run_id_array)
{
    int err = 0;
    try
	{
		std::fill_n(run_id_array, *len_run_id_array, -999);
        const std::set<int> &fail_set = _run_manager_ptr_->get_failed_run_ids();
        int n_failed_tmp = min(fail_set.size(), *len_run_id_array);
        std::copy_n(fail_set.begin(), n_failed_tmp, run_id_array);
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int RMIF_GET_NUM_TOTAL_RUNS(int *nruns)
{
    int err = 0;
    try {
        *nruns = _run_manager_ptr_->get_total_runs();
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

}