#include <stdio.h>
#include "RunManagerCWrapper.h"

main()
{
	printf("Hello World\n");
	RunManager *run_mng;
	char port[21];
	/* Body of run_manager_fortran_test */
	char storfile[21] = "tmp_run_data.bin";
	/* PANTHER parameters */
	char rmi_info_file[21] = "run_manager_info.txt";

	int nruns = 5;
	int npar = 3;
	int nobs = 16;
	char **p_names = (char*[]) { "recharge", "cond", "scoeff" };
	char **o_names = (char*[]) {
		"head1", "head2", "head3", "head4", "head5",
			"head6", "head7", "head8", "head9", "head10", "head11", "head12",
			"head13", "head14", "head15", "head16"
	};
	double pars[3] = { 0.1, .005, .05 };
	double bad_pars[3] = { 0.1, -9.0e9, .05 };

	double obs[16];
	int err = 0;
	/* instantiate PANTHER run manager */
	printf("please enter port:\n");
	scanf("%s", &port);
	printf("%s\n", port);
	/* run_mng = rmic_create_panther(storfile, port, rmi_info_file, 2, 1.15, 100.0); */
	char **command_line = (char *[]){ ".\\exe\\storage1.exe" };
	char **tpl = (char *[]){ ".\\misc\\input.tpl" };
	char **inp = (char *[]){ ".\\input.dat" };
	char **ins = (char *[]){ ".\\misc\\output.ins" };
	char **out = (char *[]){ ".\\output.dat" };
	char *run_dir = "C:\\Users\\Dave\\Documents\\GitHub\\pestpp\\benchmarks\\stor\\test";
	run_mng = rmic_create_serial(command_line, 1, tpl, 1, inp, 1, ins, 1, out, 1, storfile,run_dir, 1);

	/* intitialize run manager - allocate memory initialize parameter and observation names */
	err = rmic_initialize(run_mng, p_names, npar, o_names, nobs);

	/* add model runs to the queue */
	int run_id = -999;
	err = rmic_add_run(run_mng, bad_pars, npar, 1, &run_id);
	/* Invalid model exe index - this run should fail*/
	err = rmic_add_run(run_mng, pars, npar, 4, &run_id);
	for (int i = 0; i < nruns - 2; ++i)
	{
		pars[0] = pars[0] * 0.2;
		err = rmic_add_run(run_mng, pars, npar, 1, &run_id);
	}

	printf("  ------------------- Status of Model Runs -------------------\n");
	for (int irun = 0; irun < nruns; ++irun)
	{
		int status;
		double runtime;
		int n_concurrent;
		err = rmic_get_run_status(run_mng, irun, &status, &runtime, &n_concurrent);
		printf("    run_id = %d;   status = %d;   runtime = %f; concurrent = %d\n", irun, status, runtime, n_concurrent);
	}
	printf("  ------------------------------------------------------------ \n");

	/* perform model runs */
	printf("Performing model runs...\n");
	/*set flag to return after 15sec */
	int run_input_flag = 2;
	int run_until_no_ops = 0;
	double run_until_time_sec = .05;
	int run_output_flag = -1;
	/*
	run_input_flag parameter is used specifiy when the run_until function should return
		run_input_flag = 0  -> return after all runs are complete
		run_input_flag = 1  -> return after "run_until_no_ops" consecutive loops in the internal run manager code have not had any tcp / ip action
		run_input_flag = 2  -> return after "run_until_time_sec" seconds
		run_input_flag = 3 ->  return after "run_until_no_ops" consecutive loops in the internal run manager code have not had any tcp / ip action
		or "run_until_time_sec" seconds.Which ever comes first

		run_output_flag returns the reason the call to rmif_run_until is returning
		run_output_flag = 0->all runs are complete
		run_output_flag = 1 -> "run_until_no_ops" was active and this condition was satifisfied
		run_output_flag = 2-> "run_until_time_sec" was active and this condition occured
	*/
	while (run_output_flag != 0)
	{
		err = rmic_run_until(run_mng, run_input_flag, run_until_no_ops, run_until_time_sec, &run_output_flag);
		/* check and print status of model runs */
		printf("  ------------------- Status of Model Runs -------------------\n");
		for (int irun = 0; irun < nruns; ++irun)
		{
			int status;
			double runtime;
			int n_concurrent;
			err = rmic_get_run_status(run_mng, irun, &status, &runtime, &n_concurrent);
			printf("    run_id = %d;   status = %d;   runtime = %f; concurrent = %d\n", irun, status, runtime, n_concurrent);
		}
		printf("  ------------------------------------------------------------ \n");
		printf("  Doing useful stuff while program runs\n\n");
		/* cancel the last run */
		run_id = nruns - 1;
		err = rmic_cancel_run(run_mng, run_id); /* ony gets canceled on the first call*/
	}
	/* read results */
	for (int irun = 0; irun < nruns; ++irun)
	{
		err = rmic_get_run(run_mng, irun, pars, npar, obs, nobs);
		printf("\n\nResults for model run %d\n", irun);
		printf("Parameter Values:\n");
		for (int ipar = 0; ipar < npar; ++ipar)
		{
			printf("    %20s = %20f\n", p_names[ipar], pars[ipar]);
		}
		if (err == 0)
		{
			printf("Observation Values:\n");
			for (int iobs = 0; iobs < nobs; ++iobs)
				printf("    %20s = %20f\n", o_names[iobs], obs[iobs]);
		}
		else
		{
			printf("run failed...\n");
		}
	}

}