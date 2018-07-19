import ies_test
ies_test.setup_suite_dir("ies_10par_xsec")
ies_test.run_suite("ies_10par_xsec")
ies_test.test_chenoliver()
ies_test.tenpar_subset_test()
