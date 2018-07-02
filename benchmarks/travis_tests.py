import ies_test
#ies_test.setup_suite_dir("ies_10par_xsec")
ies_test.run_suite("ies_10par_xsec")
ies_test.test_10par_xsec(silent_master=True)
ies_test.test_chenoliver()
