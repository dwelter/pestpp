#include "RunManagerAbstract.h"

#ifndef RUNMANAGER_C_WRAP_H_
#define RUNMANAGER_C_WRAP_H_

typedef struct RunManagerAbstract RunManager;
extern "C"
{
extern __declspec(dllexport)
RunManager* rmic_create_serial(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *rundir);

extern __declspec(dllexport)
RunManager* rmic_create_yamr(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *port,
	char *info_filename);

extern __declspec(dllexport)
RunManager* rmic_create_genie(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *host,
	char *genie_tag);

extern __declspec(dllexport)
int rmic_initialize(RunManager *run_manager_ptr, 
	char **pname, int pname_array_len,
	char **oname, int oname_array_len);

extern __declspec(dllexport)
int rmic_add_run(RunManager *run_manager_ptr, double *parameter_data, int npar, int *id);

extern __declspec(dllexport)
int rmic_run(RunManager *run_manager_ptr);

extern __declspec(dllexport)
int rmic_get_run(RunManager *run_manager_ptr, int run_id, double *parameter_data, int npar, double *obs_data, int nobs);

extern __declspec(dllexport)
int rmic_delete(RunManager *run_manager_ptr);

					 
}
#endif //RUNMANAGER_C_WRAP_H_