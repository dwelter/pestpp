#PEST++
Object Oriented Inverse Modeling Software
<br><br><br>
##Overview
PEST++ can be compiled for PC, MAC, or Linux and has several run managers to support parallelization.  PEST++ implements "on-the-fly" subspace reparameterization, effectively reproducing the SVD-Assist methodology of PEST without any user intervention.  The user simply specifies how frequently a base Jacobian should be recalculated.  PEST++ also includes automatic Bayes linear parameter and forecast uncertainty estimation, which is completed for free at the end of the PEST++ run.  Also included in PEST++ are global sensitivity analysis codes that implement the method of Morris and the method of Sobol.  Both of these codes include a parallel run manager as well and are fully compatible with the PEST model-independent framework.

##Under Development
* optimization within the PEST model-independent framework, including chance constraints from Bayes Linear estimation
* global optimization and inversion with Differential Evolution
* new ways to enforce regularization as component of the composite objective function and in the linear algebra of the inversion problem

If any of these items are of interest to you, we are looking for contributors!

Visit the <a href="http://www.pestpp.org">homepage </a> to download the most recent pre-compiled version 

##Recent Updates
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


##PEST++ References:

Morris, M.D. 1991. "Factorial Sampling Plans for Preliminary Computational Experiments".  Technometrics 33(2)"161-174

Sobol, I.M. (1993), “Sensitivity Estimates for Nonlinear Mathematical Models,” Mathematical 
Modeling and Computation, 1(4):407-414. 

Welter, D.E., Doherty, J.E., Hunt, R.J., Muffels, C.T., Tonkin, M.J., and Schreüder, W.A., 2012, Approaches in highly parameterized inversion—PEST++, a Parameter ESTimation code optimized for large environmental models: U.S. Geological Survey Techniques and Methods, book 7, section C5, 47 p., available only at <a ref="http://pubs.usgs.gov/tm/tm7c5">http://pubs.usgs.gov/tm/tm7c5</a>.

###Related Links:

* <a ref="http://www.pesthomepage.org">http://www.pesthomepage.org </a>
* <a ref="http://wi.water.usgs.gov/models/pestplusplus/">http://wi.water.usgs.gov/models/pestplusplus</a>
* <a ref="http://wi.water.usgs.gov/models/genie/">http://wi.water.usgs.gov/models/genie/ </a>
* <a ref'"https://github.com/jtwhite79/pyemu">https://github.com/jtwhite79/pyemu </a>

##Compiling
The master branch includes a Visual Studio 2015 project, as well as makefiles for linux and mac.

##Testing
The benchmarks/ folder contain several test problems of varying problem size which are used to evaluate the performance of various aspects of the PEST++ algorithm and implementation.  

##Dependencies
Much work has been done to avoid additional external dependencies in PEST++.  As currently designed, the project is fully self-contained and statically linked.  

##PEST++ arguments
Here is a (more or less) complete list of ``++`` arguments that can be added to the control file
* ``++max_n_super(20)``: maximum number of super parameters to use

* ``++super_eigthres(1.0e-8)`` ratio of max to min singular values used to truncate the singular components when forming the super parameter problem

* ``++n_iter_base(1)``:number of base (full) parameter iterations to complete as part of the on-the-fly combined base-parameter/super-parameter iteration process.  A value of -1 results in calculation of the base jacobian and formation of the super parameter problem without any base parameter upgrades, replicating the behavior of the "svd-assist" methodology of PEST

* ``++n_iter_super(4)``: number of base (full) parameter iterations to complete as part of the on-the-fly combined base-parameter/super-parameter iteration process

* ``++svd_pack(propack)``: which SVD solver to use.  valid arguments are ``eigen``(jacobi solution), ``propack``(iterative Lanczos solution) and ``redsvd`` (randomized solution).

* ``++lambdas(0.1,1,10,100,1000)``: the values of lambda to test in the upgrade part of the solution process. Note that this base list is augmented with values bracketing the previous iterations best lambda.  However, if a single value is specified, only one lambda will be used.

* ``++reg_frac(0.1)``: the portion of the composite phi that will be regularization. If this argument is specified, the ``* regularization`` section of the control file is ignored.  For limited testing, values ranging from 0.05 to 0.25 seem to work well.

* ``++base_jacobian(filename)``: an existing binary jacobian file to use for the first iteration

* ``++parcov(filename)``: an ASCII PEST-style matrix file or uncertainty file to use as the prior parameter covariance matrix in FOSM uncertainty calculations and/or normal matrix scaling.  If not specified and a prior is needed, one is constructed on-the-fly from parameter bounds

* ``++uncertainty(true)``:flag to activate or deactivate FOSM-based parameter and (optionally) forecast uncertainty estimation

* ``++forecasts(fore1,fore2...)``:comma separated list of observations to treat as forecasts in the FOSM-based uncertainty estimation

* ``++iteration_summary(true)``:flag to activate or deactivate writing iteration-based CSV files summarizing parameters (<base_case>.ipar), objective function (<base_case>.iobj) and sensitivities (<base_case>.isen)

* ``++max_run_fail(4)``:maximum number of runs that can fail before the run manager emits an error.

* ``++mat_inv(jtqj)``: the form of the normal matrix to use in the solution process. Valid values are "jtqj" and "q1/2j".

* ``++der_forgive(true)``: a flag to tolerate run failures during the derivative calculation process

* ``++parcov_scale_fac(0.01)``: scaling factor to scale the prior parameter covariance matrix by when scaling the normal matrix by the inverse of the prior parameter covariance matrix.  If not specified, no scaling is undertaken; if specified, ``++mat_inv`` must be "jtqj".

###sweep ``++`` arguments
``sweep`` is a utility to run a parametric sweep for a series of parameter values.  Useful for things like monte carlo, design of experiment, etc. Designed to be used with ``pyemu`` and the python pandas library.

* ``++sweep_parameter_csv_file(filename)``: the CSV file that lists the runs to be evaluated. REQUIRED
* ``++sweep_output_csv_file(filename)``: the output CSV file from the parametric sweep.  If not passed, output is written to "output.csv"
* ``++sweep_chunk(500)``: number of runs to batch queue for the run manager.  Each chunk is read, run and written as a single batch
* ``++sweep_forgive(true)``: a flag to allow the ``sweep_parameter_csv_file`` to only include a subset of parameters listed in the control file.  If ``true``, then parameters not listed in the ``sweep_parrameter_csv_file`` are given the corresponding ``parval1`` value in the control file
* ``sweep_base_run(true)``: flag to include a "base" run of control file parameter values in the parametric sweep

###pestpp-opt ``++`` arguments
``pestpp-opt`` is a implementation of sequential linear programming under uncertainty for the PEST-style model-independent interface

* ``++opt dec var groups(<group names>)``: comma-separated string identifying which parameter groups are to be treated as decision variables.  If not passed, all adjustable parameters are treated as decision variables

* ``++opt_external_dec_var_groups(<group_names>)``: comma-separated string identifying which parameter groups are to be treated as "external" decision variables, that is decision variables that do not influence model-based constraints and therefore do not require a perturbation run of the model to fill columns in the response matrix.

* ``++opt constraint groups(<group names>)``: comma-separated string identifying which observation and prior information groups are to be treated as constraints.  Groups for "less than" constraints must start with either "l_" or "less"; groups for "greater than" constraints must start with "g_" or "greater".  If not passed, all observation and prior information groups that meet the naming rules will be treated as constraints

* ``++opt obj func(<obj func name >)``: string identifying the prior information equation or two-column ASCII file that contains the objective function coefficients.  If not passed, then each decision variable is given a coefficient of 1.0 in the objective function.

* ``++opt_direction(<direction>)``: either "min" or "max", weather to minimize or maximize the objective function. 

* ``++opt_risk(<risk>)``: a float ranging from 0.0 to 1.0 that is the value to use in the FOSM uncertainty estimation for model-based constraints. a value of 0.5 is a "risk neutral" position and no FOSM measures are calculated.  A value of 0.95 will seek a 95% risk averse solution, while a value of 0.05 will seek a 5% risk tolerant solution. See Wagner and Gorelick, 1987, *Optimal groundwater quality management under parameter uncertainty* for more background on chance-constrained linear programming
