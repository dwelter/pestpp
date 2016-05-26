#PEST++
Object Oriented Inverse Modeling Software
<br><br><br>
##Overview
PEST++ can be compiled for PC, MAC, or Linux and has several run managers to support parallelization.  PEST++ implements "on-the-fly" subspace reparameterization, effectively reproducing the SVD-Assist methodology of PEST without any user intervention.  The user simply specifies how frequently a base Jacobian should be recalculated.  PEST++ also includes automatic Bayes linear parameter and forecast uncertainty estimation, which is completed for free at the end of the PEST++ run.  Also included in PEST++ are global sensitivity analysis codes that implement the method of Morris and the method of Sobol.  Both of these codes include a parallel run manager as well and are fully compatible with the PEST model-independent framework.

##Under Developement
* optimization within the PEST model-independent framework, including chance constraints from Bayes Linear estimation
* global optimization and inversion with Differential Evolution
* new ways to enforce regularization as component of the composite objective function and in the linear algebra of the inversion problem

If any of these items are of interest to you, we are looking for contributors!

Visit the <a href="http://www.pestpp.org">homepage </a> to download the most recent pre-compiled version 

##Recent Updates
<b>Update 05/26/2016</b>: PEST++ V3 has been officially released.  It supports a number of really cool features, including global sensitivity analyses, and automatic Bayes linear (first-order, second-moment) parameter and forecast uncertainty estimates.  We also have a utility for fully-parallel parameteric sweeps from csv-based parameter files, which is useful for Monte Carlo, design of experiments, surrogate construction, etc.  All of these tools are based on the model-independent communication framework of PEST, so if you have a problem already setup, these tools are ready for you!

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
The master branch includes a Visual Studio 2013 project, as well as a Makefile.

##Testing
The benchmarks/ folder contain several test problems of varying problem size which are used to evaluate the performance of various aspects of the PEST++ algorithm and implementation.  

##Dependencies
Much work has been done to avoid additional external dependencies in PEST++.  As currently designed, the project is fully self-contained and statically linked.  
