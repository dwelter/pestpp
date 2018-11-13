from panther_pi import *

def test_serial():
    comline = [r'.\exe\storage1.exe']
    tpl = [r'.\misc\input.tpl']
    inp = [r'input.dat']
    ins = [r'.\misc\output.ins']
    out = [r'output.dat']
    stor_file = 'tmp_run_data.bin'
    #rundir = r'..\..\benchmarks\stor\test'
    rundir = r'C:\Users\Dave\Documents\GitHub\pestpp\benchmarks\stor\test'

    par_names = ['recharge', 'cond', 'scoeff']
    obs_names = ['head1', 'head2', 'head3', 'head4',
                 'head5', 'head6', 'head7', 'head8',
                 'head9', 'head10', 'head11', 'head12',
                 'head13', 'head14', 'head15', 'head16']

    run_mngr = panther()
    run_mngr.create_serial(comline, tpl, inp, ins, out, stor_file, rundir, 1)
    run_mngr.initialize(par_names, obs_names)
    #add 3 runs
    nruns = 0
    for i in range(0, 3):
        data = np.array([0.1 + .02*i, 0.005, 0.05], dtype=np.double)
        run_mngr.add_run(data, 1)
        nruns += 1
    data = np.array([-9.0e99, 0.005, 0.05], dtype=np.double)
    run_mngr.add_run(data, 1)
    nruns += 1

    err = run_mngr.run()

    failed_runs = run_mngr.get_failed_run_ids()
    print('failed runs:', failed_runs)
    for i in range(0, nruns):
        print('Results for run: %d' % i)
        par, obs = run_mngr.get_run(i)
        print(par)
        print(obs)
        print('')


def test_panther():
    stor_file = 'tmp_run_data.bin'
    port = '4005'
    print('Please start workers using port %s' % port)
    out_file = 'panther_test.out'
    par_names = ['recharge', 'cond', 'scoeff']
    obs_names = ['head1', 'head2', 'head3', 'head4',
                 'head5', 'head6', 'head7', 'head8',
                 'head9', 'head10', 'head11', 'head12',
                 'head13', 'head14', 'head15', 'head16']

    run_mngr = panther()
    run_mngr.create_panther(stor_file, port, out_file, 3, 1.15, 100.0)
    run_mngr.initialize(par_names, obs_names)
    #add 3 runs
    nruns = 3
    for i in range(0, nruns):
        data = np.array([0.1 * .02*nruns, 0.005, 0.05], dtype=np.double)
        run_mngr.add_run(data,1)

    run_mngr.run()

    for i in range(0, nruns):
        print('Results for run: %d' % i)
        par, obs = run_mngr.get_run(i)
        print(par)
        print(obs)
        print('')

print('Starting serial test')
test_serial() 
print('Starting PANTHER test')
test_panther()
    