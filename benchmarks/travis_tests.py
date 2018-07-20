import ies_test
import opt_test

opt_test.std_weights_test()
ot_test.run_supply2_test()

ies_test.setup_suite_dir("ies_freyberg")
ies_test.setup_suite_dir("ies_10par_xsec")
#ies_test.run_suite("ies_10par_xsec")

ies_test.test_freyberg_full_cov_reorder()
ies_test.test_freyberg_full_cov_reorder_run()
ies_test.test_freyberg_full_cov()
ies_test.test_freyberg_ineq()
#ies_test.freyberg_localizer_test1()
#ies_test.freyberg_localizer_test2()
ies_test.freyberg_localizer_test3()

ies_test.tenpar_subset_test()
ies_test.tenpar_full_cov_test()
ies_test.tenpar_tight_tol_test()
ies_test.tenpar_narrow_range_test()
ies_test.tenpar_fixed_test()
ies_test.tenpar_fixed_test2()
ies_test.tenpar_subset_how_test()
#ies_test.tenpar_localizer_test1()
#ies_test.tenpar_localizer_test2()
ies_test.tenpar_localizer_test3()

ies_test.test_chenoliver()
