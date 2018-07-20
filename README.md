# PEST++
Object Oriented Inverse Modeling Software
<br><br><br>
[![Travis Status](https://travis-ci.org/dwelter/pestpp.svg?branch=develop)](https://travis-ci.org/dwelter/pestpp)
[![Build status](https://ci.appveyor.com/api/projects/status/6tktqgy67d47wkjf?svg=true)](https://ci.appveyor.com/project/jtwhite79/pestpp)
## Overview
The PEST++ software suite includes several stand-alone tools for model-independent (non-intrusive) computer model parameter estimation and uncertainty analysis.  Codes include:

* ``pestpp``: deterministic GLM parameter estimation using "on-the-fly" subspace reparameterization, effectively reproducing the SVD-Assist methodology of PEST without any user intervention

* ``pestpp-gsa``: Global senitivity analysis using either Morris or Sobol

* ``pestpp-swp``: a generic parallel run utility driven by a CSV file of parameter values

* ``pestpp-opt``: chance-constrainted linear programming

* ``pestpp-ies``: iterative ensemble smoother implementation of GLM.

All members of the software suite can be compiled for PC, MAC, or Linux and have several run managers to support parallelization.  precompiled binaries are available in the "exe" folder.  Windows users should use the ``intel_c_windows`` branch binaries to avoid the dreaded MSVC missing runtime DLL issue

## Recent Updates
<b> update 4 July 2018 </b>: PESTPP++ version 4.0.0 has been released to support the newly-developed ``pestpp-ies``. A manuscript documenting ``pestpp-ies`` is available here: [https://www.sciencedirect.com/science/article/pii/S1364815218302676](https://www.sciencedirect.com/science/article/pii/S1364815218302676).  Stay tuned for an actual manual to accompany version 4!

<b> update 2 May 2018 </b>: some refactoring is underway.  ``sweep`` has been renamed ``pestpp-swp`` and ``gsa`` has been renamed ``pestpp-gsa``.  Also, the initial version of the new iterative ensemble smoother is avaiable as ``pestpp-ies``.  The basic ``++`` options needed for fine-grained control of ``pestpp-ies`` are listed below.   

<b> update 09/20/2017</b>: the new optimization under uncertainty tool is ready!  A supporting publication is in the works and should be available soon (a link will be posted once it is accepted).  This new tool uses the same control file/template file/instruction file approach as other PEST(++) applications, so applying this tool to your problem should be seamless.  Optional "++" args for tool are available further done this page.

<b>update 01/25/2017</b>: intel C++ builds are avaiable for mac and for windows.  For mac users, these are statically-linked so they do not require compilers to be installed.  For windows users, the intel build circumvents the "missing VCOMP140.DLL" error.  Note the intel windows builds are currently in the ``intel_c_windows`` branch.

<b>update 11/25/2016</b>: PEST++ version 3.6 is now available. Some of the many enhancements available in 3.6 include:

* a new approach to implementing regularization. Rather than using the standard pest control file parameters such as ``phimlim``, ``fracphim``, etc, we now offer a single pest++ argument, ``++reg_frac()``, that allows users to specify what fraction  of the composite objective function should be regularization penalty. For example, ``++reg_frac(0.5)`` would result in equal parts data misfit and regularization penalty, which results in the *maximum a posteriori* (MAP) parameter estimate. Using ``++reg_frac()`` will result in substantial speed ups during the lambda calculation process

* a new program for sequential linear programming under uncertainty.  ``pestpp-opt`` is a new executable in the PEST++ suite that uses the standard PEST model independent interface to solve a (sequential) linear programming (LP) problem.  ``pestpp-opt`` relies on the COIN-OR Linear Programming (CLP) solver [https://projects.coin-or.org/Clp]((https://projects.coin-or.org/Clp)).  Also, users have the option to use FOSM-based uncertainty estimation in the evaluation of model-based constraints (such as water levels, stream flows, stream flow depletion, etc) so that a risk-based optimal solution can be found.  See below for the required and optional ``++`` arguments needed to apply ``pestpp-opt``.  Two example problems using the ``pestpp-opt`` tool have been added to the ``benchmarks`` dir.  A publication about this tool is in the works.

* global optimization with differential evolution.  We now have a fully-parallel global solver that implements the differential evolution algorithm (DE) integrated into the pest++ executable.  See below for required and optional ``++`` arguments needed to use the DE solver.

* a new randomization-based SVD solver using the implementation of [https://github.com/ntessore/redsvd-h](https://github.com/ntessore/redsvd-h).  This solver is activated using ``++svd_pack(redsvd)``.  Testing shows it to be very efficient for problems for a wide range of problem sizes, especially with judicious use of ``max_sing``.

* upgrade parameter covariance scaling.  Through the ``++parcov_scale_fac()``, pest++ can now scale the normal matrix (J^tQJ) by a user specified parameter covariance matrix.  If no ``++parcov_filename()`` is provided, pest++ will construct a diagonal parameter covariance matrix from the parameter bounds.  This is a relatively new option and needs more testing, but limited testing to date shows that upgrade vectors resulting from a covariance-scaled normal matrix are more in harmony with expert knowledge.

<b>Update 05/26/2016</b>: PEST++ V3 has been officially released.  It supports a number of really cool features, including global sensitivity analyses, and automatic Bayes linear (first-order, second-moment) parameter and forecast uncertainty estimates.  We also have a utility for fully-parallel parametric sweeps from csv-based parameter files, which is useful for Monte Carlo, design of experiments, surrogate construction, etc.  All of these tools are based on the model-independent communication framework of PEST, so if you have a problem already setup, these tools are ready for you!

<b>Update 10/1/2014</b>: recent stable versions of PEST++ implement dynamic regularization, full restart capabilities, additional options for formulating the normal equations, and an iterative SVD algorithm for very-large problems.  Additionally the YAMR run manager has been improved to use threaded workers so that the master can more easily load balance.  
##Latest Report and Documentation
Welter, D.E., White, J.T., Hunt, R.J., and Doherty, J.E., 2015, Approaches in highly parameterized inversion— PEST++ Version 3, a Parameter ESTimation and uncertainty analysis software suite optimized for large environmental models: U.S. Geological Survey Techniques and Methods, book 7, chap. C12, 54 p., <a ref="http://dx.doi.org/10.3133/tm7C12">http://dx.doi.org/10.3133/tm7C12</a>.


## PEST++ References:

Morris, M.D. 1991. "Factorial Sampling Plans for Preliminary Computational Experiments".  Technometrics 33(2)"161-174

Sobol, I.M. (1993), “Sensitivity Estimates for Nonlinear Mathematical Models,” Mathematical 
Modeling and Computation, 1(4):407-414. 

Welter, D.E., Doherty, J.E., Hunt, R.J., Muffels, C.T., Tonkin, M.J., and Schreüder, W.A., 2012, Approaches in highly parameterized inversion—PEST++, a Parameter ESTimation code optimized for large environmental models: U.S. Geological Survey Techniques and Methods, book 7, section C5, 47 p., available only at <a ref="http://pubs.usgs.gov/tm/tm7c5">http://pubs.usgs.gov/tm/tm7c5</a>.

### Related Links:

* <a ref="http://www.pesthomepage.org">http://www.pesthomepage.org </a>
* <a ref="http://wi.water.usgs.gov/models/pestplusplus/">http://wi.water.usgs.gov/models/pestplusplus</a>
* <a ref="http://wi.water.usgs.gov/models/genie/">http://wi.water.usgs.gov/models/genie/ </a>
* <a ref="https://github.com/jtwhite79/pyemu">https://github.com/jtwhite79/pyemu </a>

## Compiling
The master branch includes a Visual Studio 2015 project, as well as makefiles for linux and mac.

## Testing
The benchmarks/ folder contain several test problems of varying problem size which are used to evaluate the performance of various aspects of the PEST++ algorithm and implementation.  

## Dependencies
Much work has been done to avoid additional external dependencies in PEST++.  As currently designed, the project is fully self-contained and statically linked.  

## pestpp arguments
Here is a (more or less) complete list of ``++`` arguments that can be added to the control file
* ``++overdue_resched_fac(1.2)``:YAMR only, if a run is more than <``overdue_resched_fac``> X average run time, reschedule it on available resources
* ``++overdue_giveup_fac(2.0)``:YAMR only, if a run is more than <``overdue_giveup_fac``> X average run time, mark it as failed
* ``++max_n_super(20)``: maximum number of super parameters to use

* ``++super_eigthres(1.0e-8)`` ratio of max to min singular values used to truncate the singular components when forming the super parameter problem

* ``++n_iter_base(1)``:number of base (full) parameter iterations to complete as part of the on-the-fly combined base-parameter/super-parameter iteration process.  A value of -1 results in calculation of the base jacobian and formation of the super parameter problem without any base parameter upgrades, replicating the behavior of the "svd-assist" methodology of PEST

* ``++n_iter_super(4)``: number of super (reduced dimension) parameter iterations to complete as part of the on-the-fly combined base-parameter/super-parameter iteration process

* ``++svd_pack(propack)``: which SVD solver to use.  valid arguments are ``eigen``(jacobi solution), ``propack``(iterative Lanczos solution) and ``redsvd`` (randomized solution).

* ``++lambdas(0.1,1,10,100,1000)``: the values of lambda to test in the upgrade part of the solution process. Note that this base list is augmented with values bracketing the previous iterations best lambda.  However, if a single value is specified, only one lambda will be used.

* ``++lambda_scale_fac(0.9,0.8,0.7,0.5)``: the values to scale each lambda upgrade vector.  This results in a line search along each upgrade vector direction, so that the number of upgrade vectors = len(lambdas) * len(lambda_scale_fac).  To disable, set = 1.0.

* ``++reg_frac(0.1)``: the portion of the composite phi that will be regularization. If this argument is specified, the ``* regularization`` section of the control file is ignored.  For limited testing, values ranging from 0.05 to 0.25 seem to work well.

* ``++base_jacobian(filename)``: an existing binary jacobian file to use for the first iteration

* ``++parcov(filename)``: an ASCII PEST-style matrix file or uncertainty file to use as the prior parameter covariance matrix in FOSM uncertainty calculations and/or normal matrix scaling.  If not specified and a prior is needed, one is constructed on-the-fly from parameter bounds

* ``++uncertainty(true)``:flag to activate or deactivate FOSM-based parameter and (optionally) forecast uncertainty estimation

* ``++forecasts(fore1,fore2...)``:comma separated list of observations to treat as forecasts in the FOSM-based uncertainty estimation

* ``++iteration_summary(true)``:flag to activate or deactivate writing iteration-based CSV files summarizing parameters (<base_case>.ipar), objective function (<base_case>.iobj) and sensitivities (<base_case>.isen), as well as upgrade summary (<base_case>.upg.csv) and a jacobian parameter-to-run_id mapping (<base_case>.rid). 
* ``++jac_scale(true)``: use PEST-style jacobian scaling. Important, but can be costly because it densifies the normal matrix, making SVD take longer.

* ``++upgrade_augment(true)``: augment the values of lambda to test by including the best lambda from the previous iteration, as well as best lambda * 2.0 and best lambda / 2.0.  If ``true``, then additional lambdas will be included by attempting to extend each upgrade vector along the region of parameter space defined by parameter bounds.  If ``false``, then only the vectors listed in the ``++lambda()`` arg will be tested and no extended upgrade will be included.  

* ``++upgrade_bounds(true)``: use additional tricks and upgrades to deal with upgrades that are going out bounds.  If ``true``, can result in substantial phi improvement, but in some cases can produce NaNs in the upgrade vectors.

* ``++hotstart_resfile(mycase.res)``: use an exising residual file to restart with an existing jacobian to forego the initial, base run and jump straight to upgrade calculations (++base_jacobian arg required).

* ``++max_run_fail(4)``:maximum number of runs that can fail before the run manager emits an error.

* ``++mat_inv(jtqj)``: the form of the normal matrix to use in the solution process. Valid values are "jtqj" and "q1/2j".

* ``++der_forgive(true)``: a flag to tolerate run failures during the derivative calculation process

* ``++parcov_scale_fac(0.01)``: scaling factor to scale the prior parameter covariance matrix by when scaling the normal matrix by the inverse of the prior parameter covariance matrix.  If not specified, no scaling is undertaken; if specified, ``++mat_inv`` must be "jtqj".

* ``++condor_submit_file(pest.sub)``: a HTCondor submit file.  Setting this arg results in use of a specialized version of the YAMR run manager where the ``condor_submit()`` command is issued before the run manager starts, and, once a set of runs are complete, the workers are released and the ``condor_rm()`` command is issued.  This specialized run manager is useful for those sharing an HTCondor pool so that during the upgrade calculation process, all workers are released and during upgrade testing, only the required number workers are queued.  As with all things PEST and PEST++, it is up to the user to make sure the relative paths between the location of the submit file, the control file and the instance of PEST++ are in sync.

### pestpp-swp ``++`` arguments
``sweep`` is a utility to run a parametric sweep for a series of parameter values.  Useful for things like monte carlo, design of experiment, etc. Designed to be used with ``pyemu`` and the python pandas library.

* ``++sweep_parameter_csv_file(filename)``: the CSV file that lists the runs to be evaluated. "sweep_in.csv" is the default
* ``++sweep_output_csv_file(filename)``: the output CSV file from the parametric sweep.  If not passed, output is written to "sweep_out.csv"
* ``++sweep_chunk(500)``: number of runs to batch queue for the run manager.  Each chunk is read, run and written as a single batch
* ``++sweep_forgive(false)``: a flag to forgive missing parameters in the input csv file.  If ``true``, then missing parameters are filled with the initial parameter value in the control file.

### pestpp-opt ``++`` arguments
``pestpp-opt`` is an implementation of sequential linear programming under uncertainty for the PEST-style model-independent interface

* ``++opt_dec_var_groups(<group names>)``: comma-separated string identifying which parameter groups are to be treated as decision variables.  If not passed, all adjustable parameters are treated as decision variables

* ``++opt_external_dec_var_groups(<group_names>)``: comma-separated string identifying which parameter groups are to be treated as "external" decision variables, that is decision variables that do not influence model-based constraints and therefore do not require a perturbation run of the model to fill columns in the response matrix.

* ``++opt_constraint_groups(<group names>)``: comma-separated string identifying which observation and prior information groups are to be treated as constraints.  Groups for "less than" constraints must start with either "l_" or "less"; groups for "greater than" constraints must start with "g_" or "greater".  If not passed, all observation and prior information groups that meet the naming rules will be treated as constraints

* ``++opt_obj_func(<obj func name >)``: string identifying the prior information equation or two-column ASCII file that contains the objective function coefficients.  If not passed, then each decision variable is given a coefficient of 1.0 in the objective function.

* ``++opt_direction(<direction>)``: either "min" or "max", whether to minimize or maximize the objective function. 

* ``++opt_risk(<risk>)``: a float ranging from 0.0 to 1.0 that is the value to use in the FOSM uncertainty estimation for model-based constraints. a value of 0.5 is a "risk neutral" position and no FOSM measures are calculated.  A value of 0.95 will seek a 95% risk averse solution, while a value of 0.05 will seek a 5% risk tolerant solution. See Wagner and Gorelick, 1987, *Optimal groundwater quality management under parameter uncertainty* for more background on chance-constrained linear programming

* ``++opt_skip_final(<skip_final>)``: a flag to skip the final model run using optimal decision variable values.  If ``true`` and ``++base_jacobian()`` and ``++hotstart_resfile()`` are set and (lot of "ands") ``noptmax``=0, then this causes no model runs to happen, pestpp-opt simply solves the chance constrainted LP problem, report optimal decision variables and phi, then exits.  Default is ``false``.

* ``++opt_std_weights(<std_weights>)``:a flag to treat model-based constraints (listed in the ``* observation data`` section) as standard deviations for chance constraints.  This can result in substantial time savings since the FOSM calculation process can be skipped.  This can also be used to specific empirical constraint uncertainty (e.g. from ensemble methods).  Default is ``false``.

### pestpp-ies ``++`` arguments
``pestpp-ies`` is an implementation of the iterative ensemble smoother GLM algorithm of Chen and Oliver 2012. So far, this tool has performed very well across a range of problems.  It functions without any additional ``++`` arguments. However, several ``++`` arguments can be used to fine-tune the function of ``pestpp-ies``.  

* ``++ies_parameter_ensemble(<par_en>)``: file containing parameter ensemble.  File extension is used to determine file type: ``.csv`` and ``.jcb`` are supported.  If not passed, a parameter ensemble is generated from prior parameter distribution

* ``++ies_observation_ensemble(<obs_en>)``: file containing observation noise ensemble (obs vals + noise realizatons). File extension is used to determine file type: ``.csv`` and ``.jcb`` are supported. If not passed, an observation ensemble is generated using observation weights

* ``ies_restart_obs_en(<restart_obs_en>)``: file containting a restart observation ensemble (simulated outputs from a previous pestpp-ies run or from pestpp-swp).  File extension is used to determine file type: ``.csv`` and ``.jcb`` are supported.

* ``ies_num_reals(<num_reals>)``: number of realizations to use.  If ``par_en``, ``obs_en`` and/or ``restart_obs_en`` are passed and are larger than ``num_reals``, then only the first ``num_reals`` realizations are used. Default is 50.  

* ``ies_bad_phi(<bad_phi>)``: realizations yielding a phi greater than ``bad_phi`` are dropped from ``par_en`` and ``obs_en``.  Default is 1.0e+30

* ``ies_lambda_mults(<lambda_mults>)``: lambda multiplers to test during upgrade calculations.  Default is [0.1,0.5,1.0,2.0,5.0].

* ``ies_initial_lambda(<init_lambda>)``: initial lambda value to use with ``lambda_mults`` to get lambda values for testing during the first iteration.  If not passed, ``init_lambda`` is derived from the initial mean phi.  

* ``ies_subset_size(<subset_size>)``: the number of realizations to test for each upgrade parameter ensemble.  The total number of upgrade testing runs is ``len(lambda_mults)`` * ``len(lambda_scale_fac)`` * ``subset_size``.  Default is 5, meaning only the first 5 realizations are evaluated during testing; if a successful (subset) upgrade ensemble is found, the remaining ``num_reals`` - ``subset_size`` realizations are evaluated.

* ``ies_reg_factor(<reg_factor>)``: fraction of zeroth-order Tikhonov regularization penalty to add to the measurement phi to form composite phi.  Default is 0.0

* ``ies_use_approx(<use_approx>)``: flag to use the approximate upgrade solution process.  Default is ``true``.

* ``ies_use_prior_scaling(<use_prior_scaling>)``: flag to scale various quantities by the prior parameter covariance matrix during upgrade calculations.  Default is ``false``.

* ``ies_use_empirical_prior(<use_empirical_prior>)``: flag to calculate prior parameter covariance matrix from the parameter ensemble.  Requires passing of ``par_en``.  Default is ``false``.

* ``ies_verbose_level(<verbose_level>)``: integer to control how much pestpp-ies writes.  Can be 0, 1, 2, or 3.  Default is 1.

* ``ies_add_bases(<add_bases>)``: flag to add initial parameter values to ``par_en`` and add actual observation values (no noise) to ``obs_en``.  This results in an approximation to the minimum error variance parameter set being carried through the pestpp-ies analysis.  If ``true``, results in ``num_reals`` + 1 realizations.  Default is ``true``.

* ``ies_enforce_bounds(enforce_bounds>)``: flag to enforce parameter bounds during upgrade calculations.  Default is ``true``.

* ``ies_save_binary(<save_binary>)``: flag to save iteration parameter and observation ensembles to pest-compatible (jacobian format) binary files.  Default is ``false``

* ``ies_accept_phi_fac(<accept_phi_fac>)``: tolerance for accepting the results for a (subset) ensemble evaluation. If the resulting mean phi * ``accept_phi_fac`` is greater than the best mean phi from the last iteration, then the upgrade is rejected.  Default is 1.05 (5% tolerance).

* ``ies_lambda_inc_fac(<lambda_inc_fac>)``: factor increase current lambda by if current upgrade testing was not successful.  Default is 10.0

* ``ies_lambda_dec_fac(<lambda_dec_fac>)``: factor to decrease current lambd by if current upgrade testing was successful.  Default is 0.75.

* ``ies_save_lambda_en(<save_lambda_en>)``: flag to save lambda testing parameter ensembles.  Can be use for finding parameters vectors that are causing run failures.  Default is ``false``.

* ``parcov_filename(<parcov_filename>)``: (repeated argument for pestpp and pestpp-opt).  Name of existing prior parameter covariance matrix.  Can be ASCII (``.mat`` or ``.cov``), binary (``.jco`` or ``.jcb``) or an uncertainty file (``.unc``).  If not passed, the prior parameter covariance matrix is constructed from parameter bounds (or from the parameter ensemble if ``use_empirical_prior`` is passed).  

* ``par_sigma_range(<par_sigma_range>)``: the number of standard deviations implied by parameter bounds.  Used to construct prior parameter covariance matrix from parameter bounds.  Default is 4.0.

* ``lambda_scale_fac(<lambda_scale_fac>)``: line search lambda scaling factors.  Each lambda parameter upgrade ensemble is scaled by each ``lambda_scale_fac``.  Default is [0.5,0.75,0.9,1.0,1.1].  

* ``ies_subset_how(<subset_how>)``: choice for how the subset is selected.  Choices are "first" (use the first ``subset_size`` realizations),"last" (use the last ``subset_size`` realizations),"random" (randomly select ``subset_size`` realizations each iteration),"phi_based" (chose ``subset_size`` realizations spread across the composite phi from the last iteration). Default is "first".  

* ``ies_localizer(<localizer>)``: an optional localizer to use to localize spurious cross-correlations between observations and parameters.  The file format is determined by the extension: "mat","cov","csv","jco"/"jcb".  The row names of the matrix are observation names and/or observation group names and columns are parameter names and/or parameter group names (each obs/par can only be include once even through groups).  Currently only zero vs nonzero entires are used to localize.  

### USGS disclaimer

This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use
