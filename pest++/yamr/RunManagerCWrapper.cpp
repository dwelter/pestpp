#include "RunManagerCWrapper.h"
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

RunManager* rmic_create_serial(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *rundir)
{
	RunManager *run_manager_ptr = nullptr;
	vector<string> comline_vec(comline, comline+comline_array_len);
	vector<string> tpl_vec(tpl, tpl+tpl_array_len);
	vector<string> inp_vec(inp, inp+inp_array_len);
	vector<string> ins_vec(ins, ins+ins_array_len);
	vector<string> out_vec(out, out+out_array_len);
	run_manager_ptr = new RunManagerSerial(comline_vec, tpl_vec, inp_vec, ins_vec, out_vec, storfile, rundir);	
	return run_manager_ptr;
}


RunManager* rmic_create_yamr(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *port,
	char *info_filename)
{
	RunManager *run_manager_ptr = nullptr;
	vector<string> comline_vec(comline, comline+comline_array_len);
	vector<string> tpl_vec(tpl, tpl+tpl_array_len);
	vector<string> inp_vec(inp, inp+inp_array_len);
	vector<string> ins_vec(ins, ins+ins_array_len);
	vector<string> out_vec(out, out+out_array_len);
	fout_run_manager_log_file.open(info_filename);
	run_manager_ptr = new RunManagerYAMR(comline_vec, tpl_vec, inp_vec, ins_vec, out_vec, storfile, port, fout_run_manager_log_file);
	return run_manager_ptr;
}

RunManager* rmic_create_genie(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *host,
	char *genie_tag)
{
	RunManager *run_manager_ptr = nullptr;
	vector<string> comline_vec(comline, comline+comline_array_len);
	vector<string> tpl_vec(tpl, tpl+tpl_array_len);
	vector<string> inp_vec(inp, inp+inp_array_len);
	vector<string> ins_vec(ins, ins+ins_array_len);
	vector<string> out_vec(out, out+out_array_len);
	run_manager_ptr = new RunManagerGenie(comline_vec, tpl_vec, inp_vec, ins_vec, out_vec, storfile, genie_tag);
	return run_manager_ptr;
}

int rmic_initialize(RunManager *run_manager_ptr, 
	char **pname, int pname_array_len,
	char **oname, int oname_array_len)

{
	int err = 0;
	try {
		vector<string> pname_vec(pname, pname+pname_array_len);
		vector<string> oname_vec(oname, oname+oname_array_len);
		run_manager_ptr->initialize(pname_vec, oname_vec);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmic_add_run(RunManager *run_manager_ptr, double *parameter_data, int npar, int *id)
{
	int err = 0;
	try {
		vector<double> data(parameter_data, parameter_data+npar);
		*id = run_manager_ptr->add_run(data);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}



int rmic_run(RunManager *run_manager_ptr)
{
	int err = 0;
	try {
		run_manager_ptr->run();
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmic_get_run(RunManager *run_manager_ptr, int run_id, double *parameter_data, int npar, double *obs_data, int nobs)
{
	int err = 0;
	try {
	err = run_manager_ptr->get_run(run_id, parameter_data, npar, obs_data, nobs);
	}
	catch(PestIndexError ex) {
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int rmic_delete(RunManager *run_manager_ptr)
{
	int err = 0;
	try {
		delete run_manager_ptr;
	}
		catch(...)
	{
		err = 1;
	}
	return err;
}

}