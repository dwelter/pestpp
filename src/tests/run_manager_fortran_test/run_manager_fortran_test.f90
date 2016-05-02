!  run_manager_fortran_test.f90 
!
!  FUNCTIONS:
!  run_manager_fortran_test - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: run_manager_fortran_test
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program run_manager_fortran_test

    implicit none
     ! Variables
     
    external rmif_create_serial
    integer rmif_create_serial
    external rmif_create_yamr
    integer rmif_create_yamr
    external rmif_create_genie
    integer rmif_create_genie
    external rmif_initialize
    integer rmif_initialize
    external rmif_reinitialize
    integer rmif_reinitialize
    external rmif_add_run
    integer rmif_add_run
    external rmif_run
    integer rmif_run
    external rmif_get_run
    integer rmif_get_run
    external rmif_get_num_failed_runs
    integer rmif_get_num_failed_runs
    external rmif_get_failed_run_ids
    integer rmif_get_failed_run_ids
    external rmif_get_num_total_runs
    integer rmif_get_num_total_runs
    
    
    external rmif_delete
    integer rmif_delete
    
    character*20 comline(1)
    character*20 tpl(1)
    character*20 inp(1)
    character*20 ins(1)
    character*20 out(1)
    character*20 storfile
    character*20 port
    character*20 rmi_info_file
    character*80 rundir
    character*20 genie_host
    character*80 genie_tag
    character*20 p_names(3)
    character*50 o_names(16)
    character buf
    DOUBLE PRECISION pars(3)
    DOUBLE PRECISION bad_pars(3)
    DOUBLE PRECISION obs(16)
    integer nruns, npar, nobs
    integer itype
    integer err
    integer ipar, iobs, irun
    integer nfail
    integer failed_run_ids(100)
    integer n_total_runs
    
    
    ! Body of run_manager_fortran_test
    data comline   /'storage1_win.exe  '/
    data tpl       /'input.tpl           '/
    data inp       /'input.dat           '/
    data ins       /'output.ins          '/
    data out       /'output.dat          '/
    data storfile  /'tmp_run_data.bin    '/
    !YAMR parameters
    data rmi_info_file  /'run_manager_info.txt'/
    !serial run manager parameters
    data rundir    /'.\                  '/
    !genie run manager parameters
    data genie_host      /'localhost:4005      '/
    data genie_tag           /'rmif test           '/
    
    data p_names   /'recharge            ',&
                    'cond                ',&
                    'scoeff              '/
    
    data o_names   /'head1               ',&
                    'head2               ',&
                    'head3               ',&
                    'head4               ',&
                    'head5               ',&
                    'head6               ',&
                    'head7               ',&
                    'head8               ',&
                    'head9               ',&
                    'head10              ',&
                    'head11              ',&
                    'head12              ',&
                    'head13              ',&
                    'head14              ',&
                    'head15              ',&
                    'head16              '/
    
    data pars / 0.1, .005, .05/
    data bad_pars / 0.1, -9e9, .05/
    
    ! instantiate run manager
    itype = -1
    do while (itype < 1)
        write(*,*) "Select a run manager"
        write(*,*) "  1 - serial run manager"
        write(*,*) "  2 - YAMR run manager"
        write(*,*) "  3- GENIE run manager"
        write(*,*) ""
        read(*,*) itype
        if (itype == 1) then
             write(*,*) 'please enter model run directory: '
             read(*,*) rundir
             err = rmif_create_serial(comline, 20, 1,&
                        tpl, 20, 1,&
                        inp, 20, 1,&
                        ins, 20, 1,&
                        out, 20, 1,&
                        storfile, 20,&
                        rundir, 80, 3)
        else if (itype == 2) then
            write(*,*) 'please enter port:'
            read(*,*) port
            err = rmif_create_yamr(storfile, 20,&
                        port, 20,&
                        rmi_info_file, 20, 2, 1.15, 100.0)
        else if (itype == 3) then
            write(*,*) 'please enter GMAN socket:'
            read(*,*) genie_host
            write(*,*) 'please enter GENIE id tag:'
            read(*,*) genie_tag
            err = rmif_create_genie(comline, 20, 1,&
                        tpl, 20, 1,&
                        inp, 20, 1,&
                        ins, 20, 1,&
                        out, 20, 1,&
                        storfile, 20,&
                        genie_host, 20,&
                        genie_tag, 80 )
        else
                itype = -1
        end if
    end do 
            
    
    
    nruns = 4
    npar = 3
    nobs = 16
    ! intitialize run manager - allocate memory initialize parameter and observation names
    err = rmif_initialize(p_names, 20, npar, o_names, 50, nobs)
    
    ! add model runs to the queue
    err = rmif_add_run(bad_pars, npar, irun)
    do irun = 1, nruns - 1
        pars(1) = pars(1) + pars(1) * 0.2
        err = rmif_add_run(pars, npar, irun)
    end do
    
    ! perform mdoel runs
    write(*,*) 'Performing model runs...'
    err = rmif_run()
    
    ! get number of failed model runs
    err = rmif_get_num_failed_runs(nfail)
    write(*,*) 'Number of failed runs: ', nfail
    err = rmif_get_failed_run_ids(failed_run_ids, 100)
    do irun = 1, nfail
         write(*,*) '   failed run id = ', failed_run_ids(irun)
    end do    
    
    ! read results
    do irun = 0, nruns-1
        err = rmif_get_run(irun, pars,npar, obs, nobs)
        
        write(*,*) ""
        write(*,*) ""
        write(*,*) "Results for model run", irun
        write(*,*) "Parameter Values:"
        do ipar = 1, npar
            write(*,*) "  parameter ", p_names(ipar), " = ", pars(ipar)
        end do
        if (err ==0 ) then
            write(*,*) "Observation Values:"
            do iobs = 1, nobs
                write(*,*) "  observation ", o_names(iobs), " = ", obs(iobs)
            end do
        else
            write(*,*) "run failed..."
        endif
    end do
    
    err = rmif_get_num_total_runs(n_total_runs)
    write(*,*) ''
    write(*,*) 'Total number of successful model runs:', n_total_runs
    
    ! reinitialize run manager and make another set of runs
    err = rmif_reinitialize()
    err = rmif_add_run(pars, npar, irun)
    pars(1) = pars(1) + pars(1) * 0.2
    err = rmif_add_run(pars, npar, irun)
    err = rmif_run()
    
    err = rmif_get_run(0, pars,npar, obs, nobs)
    
    err = rmif_get_num_total_runs(n_total_runs)
    write(*,*) ''
    write(*,*) 'Total number of successful model runs:', n_total_runs
    
    !clean up
    err = rmif_delete()
   
    write (*,*) ""
    write (*,*) "Press ENTER to continue"
    read (*,*)
    end program run_manager_fortran_test

