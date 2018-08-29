#ifndef RUNMANAGER_C_WRAP_H_
#define RUNMANAGER_C_WRAP_H_

typedef struct RunManagerAbstract RunManager;

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		RunManager* rmic_create_serial(char **comline, int comline_array_len,
			char **tpl, int tpl_array_len,
			char **inp, int inp_array_len,
			char **ins, int ins_array_len,
			char **out, int out_array_len,
			char *storfile,
			char *rundir,
			int n_max_fail);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		RunManager* rmic_create_panther(char *storfile,
			char *port,
			char *info_filename,
			int n_max_fail,
			double overdue_reched_fac, double overdue_giveup_fac);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		const char* rmic_err_msg();

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_initialize(RunManager *run_manager_ptr,
			char **pname, int pname_array_len,
			char **oname, int oname_array_len);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_reinitialize(RunManager *run_manager_ptr);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_initialize_restart_(RunManager *run_manager_ptr, char *store_filename);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_add_run(RunManager *run_manager_ptr, double *parameter_data, int npar, int model_exe_index, int *id);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_add_run_with_info(RunManager *run_manager_ptr, double *parameter_data, int npar, 
			int model_exe_index, char *info_txt, double info_value, int *id);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_run(RunManager *run_manager_ptr);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_run_until(RunManager *run_manager_ptr, int condition, int no_ops, double time_sec, int * return_cond);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_cancel_run(RunManager *run_manager_ptr, int run_id);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_get_run_status_info(RunManager *run_manager_ptr, int run_id, int *run_status, double *max_runtime, int *n_concurrent_runs);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_get_run(RunManager *run_manager_ptr, int run_id, double *parameter_data, int npar, double *obs_data, int nobs);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_get_run_with_info(RunManager *run_manager_ptr, int run_id, double *parameter_data, int npar, double *obs_data, int nobs, char *info_txt, int info_txt_len, double info_value);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_get_n_failed_runs(RunManager *run_manager_ptr, int *nfail);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_get_failed_run_ids_alloc(RunManager *run_manager_ptr, int *run_id_array, int *nfail);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_get_failed_run_ids(RunManager *run_manager_ptr, int *run_id_array, int nfail);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_get_n_cur_runs(RunManager *run_manager_ptr, int *nruns);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_get_n_total_runs(RunManager *run_manager_ptr, int *total_runs);

#ifdef OS_WIN
	extern __declspec(dllexport)
#endif
		int rmic_delete(RunManager *run_manager_ptr);

#ifdef __cplusplus
}
#endif					 

#endif //RUNMANAGER_C_WRAP_H_
