program neb

#ifdef USE_MPI
  use mpi
  use slave
#endif
  use utility
  use geometry
  use idpp
  use pes
  use elastic
  use climbing
  use optimization
  use output
  use input
  use computation_info

  implicit none

  integer       :: i
  integer       :: cmdargcount
  integer(LONG) :: end_clock
  character(80) :: arg
  character(80) :: cmd_name
  character(80) :: fname_in
  logical       :: fmain_converged
#ifdef USE_MPI
  integer       :: priv_comm_sz
  integer       :: priv_proc_id
  integer       :: err_n
#endif

  !==================================================================
  ! Get compilation info
  !==================================================================

  ! if compiled with MPI ----------------------------------
#ifdef USE_MPI
  call MPI_Init(err_n)
  call MPI_Comm_size(MPI_COMM_WORLD,priv_comm_sz,err_n)
  call MPI_Comm_rank(MPI_COMM_WORLD,priv_proc_id,err_n)

  call update_comm_sz(priv_comm_sz)
  call update_proc_id(priv_proc_id)
#endif

  !==================================================================
  ! Master/Slaves branching
  !==================================================================

  ! all slaves await commands -----------------------------
#ifdef USE_MPI
  if (proc_id/=0) then
    call slave_on_idle()
  end if
#endif

  ! only the master continue with the main execution ------

  !==================================================================
  ! Arguments Parsing
  !==================================================================

  ! get program name --------------------------------------
  call get_command_argument(0,cmd_name)
  call set_main_program_name(cmd_name)
  
  ! get number of arguments -------------------------------
  cmdargcount = command_argument_count()

  if (cmdargcount<1) then
    call write_help(cmd_name)
    call end_main_exec() ! no arguments
  end if

  ! parse arguments ---------------------------------------
  do i=1, cmdargcount
    call get_command_argument(i,arg)

    select case (arg)
    case ("-h","--help")
      call write_help(cmd_name)
      call end_main_exec() ! help argument
    case ("-v","--version")
      call write_build_version()
      call end_main_exec() ! version argument
    case ("-t","--template")
      call write_neb_input_template(.false.)
      call end_main_exec() ! template argument
    case ("--verbose-template")
      call write_neb_input_template(.true.)
      call end_main_exec() ! verbose template argument
    case default
      if (i==cmdargcount) then
        fname_in = arg
      else
        call error("main: wrong argument: "//arg)
      end if
    end select
  end do

  !==================================================================
  ! Open output stream
  !==================================================================

  call set_fileout(fname_in)

#ifdef DBG0
  !DEBUG --------------------------------------------------
  ! This section is intended to check the results
  ! of compute_total_forces subroutine.
  !--------------------------------------------------------

  write(FILEOUT,'(" *** main: ",A," Debug Version")') trim(cmd_name)
  call write_build_version()
  write(FILEOUT,'(5X,"Intended to check the results of compute_total_forces")')

  ! start time
  call set_start_clock()
  call write_date(trim(cmd_name)//" Launched on")

  write(FILEOUT,*) "*** main: input file: ", trim(fname_in)

  write(FILEOUT,*) "*** main: input_dumping"
  call input_dumping(fname_in)

  write(FILEOUT,*) "*** main: read_input"
  call read_input(fname_in)

  write(FILEOUT,*) "*** main: init_images"
  call init_images()
  call last_geom_bkp(.false.,"linear-interpolation.xyz")

  write(FILEOUT,*) "*** main: init_pes_module"
  call init_pes_module()

  write(FILEOUT,*) "*** main: init_elastic_module"
  call init_elastic_module()

  write(FILEOUT,*) "*** main: write_parallelization_info"
  call write_parallelization_info()

  write(FILEOUT,*) "*** main: write_pes_info"
  call write_pes_info()

  if (flag_idpp) then
    write(FILEOUT,*) "*** main: init_idpp"
    call init_idpp()

    write(FILEOUT,*) "*** main: optmz_idpp"
    call optmz_idpp(fmain_converged)
    if (fmain_converged.eqv..false.) then
      call error("main: optmz_idpp did not converge")
    end if
    call last_geom_bkp(.false.,"idpp-interpolation.xyz")
  end if

  write(FILEOUT,*) "*** main: write_all_images"
  call write_all_images(FILEOUT,.true.)
  
  write(FILEOUT,*) "*** main: compute_total_forces"
  call compute_total_forces(PES_MODE,.true.)

  write(FILEOUT,*) "*** main: write_pes_energy"
  call write_pes_energy()

  write(FILEOUT,*) "*** main: write_pes_forces"
  call write_pes_forces()

  write(FILEOUT,*) "*** main: write_perpen_pes_forces"
  call write_perpen_pes_forces()

  write(FILEOUT,*) "*** main: write_parall_elastic_forces"
  call write_parall_elastic_forces()

  write(FILEOUT,*) "*** main: write_total_forces"
  call write_total_forces()

  write(FILEOUT,*) "*** main: Successful Termination"

  ! end time
  call system_clock(end_clock)
  call human_time("Total Execution Time",start_clock,end_clock,clock_rate)
  call write_date(trim(cmd_name)//" Terminated on")

  call close_fileout()
  call end_main_exec() ! DBG0
#endif

  !==================================================================
  ! Working Section
  !==================================================================

  write(FILEOUT,'(" *** main: ",A," Standard Version")') trim(cmd_name)
  call write_build_version()

  ! start time
  call set_start_clock()
  call write_date(trim(cmd_name)//" Launched on")

  write(FILEOUT,*) "*** main: input file: ", trim(fname_in)

  write(FILEOUT,*) "*** main: input_dumping"
  call input_dumping(fname_in)

  write(FILEOUT,*) "*** main: read_input"
  call read_input(fname_in)

  write(FILEOUT,*) "*** main: init_images"
  call init_images()
  call last_geom_bkp(.false.,"linear-interpolation.xyz")

  write(FILEOUT,*) "*** main: init_pes_module"
  call init_pes_module()

  write(FILEOUT,*) "*** main: init_elastic_module"
  call init_elastic_module()

  write(FILEOUT,*) "*** main: write_parallelization_info"
  call write_parallelization_info()

  write(FILEOUT,*) "*** main: write_pes_info"
  call write_pes_info()

  if (flag_idpp) then
    write(FILEOUT,*) "*** main: init_idpp"
    call init_idpp()
    
    write(FILEOUT,*) "*** main: optmz_idpp"
    call optmz_idpp(fmain_converged)
    if (fmain_converged.eqv..false.) then
      call error("main: optmz_idpp did not converge")
    end if
    call last_geom_bkp(.false.,"idpp-interpolation.xyz")
  end if

  write(FILEOUT,*) "*** main: write_all_images"
  write(FILEOUT,*) "*** main: Initial Geometries"
  call write_all_images(FILEOUT,.true.)

  if (flag_only_interpolation) then
    write(FILEOUT,*) "*** main: Specified: Only Interpolation"
    write(FILEOUT,*) "*** main: Successful Termination"
  else
    write(FILEOUT,*) "*** main: optmz_pes"
    call optmz_pes(fmain_converged)

    write(FILEOUT,*) "*** main: write_all_images"
    write(FILEOUT,*) "*** main: Final Geometries"
    call write_all_images(FILEOUT,.true.)

    if (fmain_converged) then
      write(FILEOUT,*) "*** main: write_transition_state"
      call write_transition_state()

      call last_geom_bkp(.false.,"final_geometries.xyz")

      write(FILEOUT,*) "*** main: Successful Termination :: Convergence Achieved"
    else
      write(FILEOUT,*) "*** main: Successful Termination :: Convergence NOT Achieved"
    end if
  end if

  ! end time
  call system_clock(end_clock)
  call human_time("Total Execution Time",start_clock,end_clock,clock_rate)
  call write_date(trim(cmd_name)//" Terminated on")

  !==================================================================
  ! Close output stream
  !==================================================================

  call close_fileout()

  !==================================================================
  ! Finalize
  !==================================================================

  call end_main_exec() ! main

contains

!====================================================================

subroutine write_help(cmd_name)

  character(*), intent(IN) :: cmd_name
  
  write(STDOUT,'(A)') "Usage: "//trim(cmd_name)//" [OPTION] FILE"
  write(STDOUT,'(A)')
  write(STDOUT,'(A)') "OPTION list"
  write(STDOUT,'(A)') "  -h, --help                Print this help and exit."
  write(STDOUT,'(A)') "  -t, --template            Write an input template and exit."
  write(STDOUT,'(A)') "      --verbose-template    Write a commented input template and exit."
  write(STDOUT,'(A)') "                              Very useful if this is your first time using the program."
  write(STDOUT,'(A)') "  -v, --version             Print version info and exit."

end subroutine write_help

!====================================================================

subroutine write_build_version()
  
  integer :: out_stream
#ifdef USE_MPI
  integer :: version
  integer :: subversion
#endif

  if (flag_fileout) then
    out_stream=FILEOUT
  else
    out_stream=STDOUT
  end if

  if (flag_mpi) then
    write(out_stream,'("Parallel Version")')
#ifdef USE_MPI
    call mpi_get_version(version,subversion,err_n)
    write(out_stream,'("Distributed Memory via MPI (version ",&
      &I1,".",I1,")")') version, subversion
#endif
  else
    write(out_stream,'("Serial Version")')
  end if

end subroutine write_build_version

!====================================================================

end program neb

