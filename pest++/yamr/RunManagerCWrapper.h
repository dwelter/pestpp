#include "RunManagerAbstract.h"

#ifndef RUNMANAGER_C_WRAP_H_
#define RUNMANAGER_C_WRAP_H_

typedef struct RunManagerAbstract RunManager;
extern "C"
{
RunManager* rmic_create_serial(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *rundir);

RunManager* rmic_create_yamr(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *port,
	char *genie_tag,
	char *info_filename);

RunManager* rmic_create_genie(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *host,
	char *genie_tag);

int rmic_initialize(RunManager *run_manager_ptr, 
	char **pname, int pname_array_len,
	char **oname, int oname_array_len);

int rmic_add_run(RunManager *run_manager_ptr, double *parameter_data, int npar, int id);

int rmic_run(RunManager *run_manager_ptr);

int rmic_get_run(RunManager *run_manager_ptr, int run_id, double *parameter_data, int npar, double *obs_data, int nobs);

int rmic_delete(RunManager *run_manager_ptr);

					 
}
#endif //RUNMANAGER_C_WRAP_H_