#include "RunManagerAbstract.h"

#ifndef RUNMANAGER_FORTRAN_WRAP_H_
#define RUNMANAGER_FORTRAN_WRAP_H_

extern "C"
{
void RMIF_CREATE_SERIAL(char *f_comline, int  *comline_str_len, int *comline_array_len,
					 char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
					 char *f_inp, int  *inp_str_len, int *inp_array_len,
					 char *f_ins, int  *ins_str_len, int *ins_array_len,
					 char *f_out, int  *out_str_len, int *out_array_len,
					 char *f_storfile, int *storfile_len,
					 char *f_rundir, int *rundir_len);

void RMIF_CREATE_YAMR(char *f_comline, int  *comline_str_len, int *comline_array_len,
					 char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
					 char *f_inp, int  *inp_str_len, int *inp_array_len,
					 char *f_ins, int  *ins_str_len, int *ins_array_len,
					 char *f_out, int  *out_str_len, int *out_array_len,
					 char *f_storfile, int *storfile_len,
					 char *f_port, int *f_port_len,
					 char *f_info_filename, int *info_filename_len );

void RMIF_CREATE_GENIE(char *f_comline, int  *comline_str_len, int *comline_array_len,
					 char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
					 char *f_inp, int  *inp_str_len, int *inp_array_len,
					 char *f_ins, int  *ins_str_len, int *ins_array_len,
					 char *f_out, int  *out_str_len, int *out_array_len,
					 char *f_storfile, int *storfile_len,
					 char *f_host, int *f_host_len,
					 char *f_genie_tag, int *genie_tag_len);

int RMIF_ADD_RUN(double *parameter_data, int *npar);

int RMIF_INITIALIZE(char *f_pname, int  *pname_str_len, int *pname_array_len,
				 char *f_oname, int  *oname_str_len, int *oname_array_len);

void RMIF_RUN();

int RMIF_GET_RUN(int *run_id, double *parameter_data, int *npar, double *obs_data, int *nobs);

void RMFI_DELETE();

					 
}
#endif //RUNMANAGER_FORTRAN_WRAP_H_