! Copyright (C) 2020-2021  Marco Bertini
!
! This file is part of Empathes.
!
! Empathes is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Empathes is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Empathes.  If not, see <https://www.gnu.org/licenses/>.

module pes

  !==================================================================
  !   PES Module
  !==================================================================
  ! This module executes external programs
  ! from which gets energy and force values.
  !==================================================================

#ifdef USE_MPI
  use mpi
#endif
  use utility
  use c_utility_wrapper
  use geometry
  use pes_data
  use pes_input_template

  implicit none
  save
  private

  ! protected variables -----------------------------------
  public    :: pes_program,               &
               pes_exec,                  &
               pes_proc,                  &
               pes_forces,                &
               pes_energy,                &
               start_energy,              &
               end_energy,                &
               flag_init_pes_module,      &
               flag_new_pes_program
  protected :: pes_program,               &
               pes_exec,                  &
               pes_proc,                  &
               pes_forces,                &
               pes_energy,                &
               start_energy,              &
               end_energy,                &
               flag_init_pes_module,      &
               flag_new_pes_program
  ! public procedures -------------------------------------
  public    :: init_pes_module,           &
               set_new_pes_program,       &
               set_start_energy,          &
               set_end_energy,            &
               set_pes_program,           &
               set_pes_program_with_mpi,  &
               set_pes_exec,              &
               set_pes_proc,              &
               set_pes_mem,               &
               compute_pes_forces,        &
               get_pes_forces,            &
               update_pes_energy,         &
               update_pes_forces,         &
               get_scfconv,               &
               get_scfcycle

  !--------------------------------------------------------
  character(*), parameter                   :: base_name                  = "neb"
  character(*), parameter                   :: base_dirname               = "dir_neb"
  character(*), parameter                   :: base_auxdirname            = "dir_aux"
  logical                                   :: flag_init_pes_module       = .false.
  logical                                   :: flag_new_pes_program       = .false.
  logical                                   :: flag_pes_program           = .false.
  logical                                   :: flag_pes_program_with_mpi  = .false.
  logical                                   :: flag_pes_exec              = .false.
  logical                                   :: flag_pes_proc              = .false.
  logical                                   :: flag_pes_mem               = .false.
  logical                                   :: flag_start_energy          = .false.
  logical                                   :: flag_end_energy            = .false.
  character(30)                             :: pes_program                = "null"
  character(120)                            :: pes_exec
  integer                                   :: pes_proc                   = 1
  integer                                   :: pes_mem
  character                                 :: pes_mem_scale
  real(DBL), allocatable, dimension(:)      :: pes_energy
  real(DBL), allocatable, dimension(:,:)    :: pes_forces
  real(DBL)                                 :: start_energy
  real(DBL)                                 :: end_energy

contains

!====================================================================
! Public
!====================================================================

subroutine set_new_pes_program(flag)

  logical,     intent(IN) :: flag

  character(*), parameter :: my_name    = "set_new_pes_program"
  logical, save           :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  flag_new_pes_program = flag

  first_call = .false.

end subroutine set_new_pes_program

!====================================================================

subroutine set_start_energy(str)

  character(*), intent(IN) :: str

  if (flag_start_energy) then
    call error("set_start_energy: start energy already setted")
  end if

  if (.not.isreal(trim(adjustl(str)))) then
    call error("set_start_energy: argument """//str//""" is not valid")
  end if

  read(str,*) start_energy

  flag_start_energy = .true.

end subroutine set_start_energy

!====================================================================

subroutine set_end_energy(str)

  character(*), intent(IN) :: str

  if (flag_end_energy) then
    call error("set_end_energy: end energy already setted")
  end if

  if (.not.isreal(trim(adjustl(str)))) then
    call error("set_end_energy: argument """//str//""" is not valid")
  end if

  read(str,*) end_energy

  flag_end_energy = .true.

end subroutine set_end_energy

!====================================================================

subroutine set_pes_program(str)

  character(*), intent(INOUT) :: str

#ifdef USE_MPI
  character(8)                :: istr
  logical                     :: flag_warning
#endif

  ! preliminary checks ------------------------------------
  if (flag_pes_program) then
    call error("set_pes_program: pes program already setted")
  end if

  if (len_trim(str)==0) then
    call error("set_pes_program: pes program not specified")
  end if

#ifdef USE_MPI
  flag_warning = .false.
#endif
  call tolower(str)
  pes_program=str

  ! check that pes_program is valid -----------------------
  select case (pes_program)
  case ("null")
    call error("set_pes_program: pes program not specified, use #PESPROGRAM")
  case ("gaussian")
  case ("siesta")
#ifdef USE_MPI
    if (flag_pes_program_with_mpi) then
      flag_warning = .true.
    end if
#endif
  case default
    if (flag_new_pes_program.eqv..false.) then
      call error("set_pes_program: invalid option """//trim(pes_program)//"""")
    end if
  end select

  ! final checks ------------------------------------------
#ifdef USE_MPI
  if ((flag_warning).and.(comm_sz>1)) then
    write(istr,'(I8)') comm_sz
    istr=adjustl(istr)
    call error("set_pes_program: if you want to run the program """//trim(pes_program)//&
      &""" with MPI, "//trim(main_program_name)//&
      &" must be executed with 1 process (executed with "//trim(istr)//" processes)")
  end if
#endif

  flag_pes_program = .true.

end subroutine set_pes_program

!====================================================================

subroutine set_pes_program_with_mpi(flag)

  logical,     intent(IN) :: flag

  character(*), parameter :: my_name = "set_pes_program_with_mpi"
  logical, save           :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  flag_pes_program_with_mpi = flag

  first_call = .false.

end subroutine set_pes_program_with_mpi

!====================================================================

subroutine set_pes_exec(str)

  character(*), intent(IN) :: str

  if (flag_pes_exec) then
    call error("set_pes_exec: pes exec already setted")
  end if

  if (len_trim(str)==0) then
    call error("set_pes_exec: pes exec not specified")
  else
    pes_exec = str
  end if

  flag_pes_exec = .true.

end subroutine set_pes_exec

!====================================================================

subroutine set_pes_proc(str)

  character(*), intent(IN) :: str

  if (flag_pes_proc) then
    call error("set_pes_proc: pes processors already setted")
  end if

  if (len_trim(str)==0) then
    call error("set_pes_proc: pes processors not specified")
  end if

  if (isinteger(trim(adjustl(str)))) then
    read(str,*) pes_proc
  else
    call error("set_pes_proc: wrong argument """//trim(str)//"""")
  end if

  if (pes_proc<=0) then
    call error("set_pes_proc: argument must be a non-zero positive integer")
  end if

  flag_pes_proc = .true.

end subroutine set_pes_proc

!====================================================================

subroutine set_pes_mem(str)

  character(*), intent(IN) :: str

  integer                  :: i
  character                :: c

  if (flag_pes_mem) then
    call error("set_pes_mem: pes memory already setted")
  end if

  if (len_trim(str)==0) then
    call error("set_pes_mem: pes memory not specified")
  end if

  i=len_trim(str)

  if (isalpha(str(i:i))) then
    c = str(i:i)
    call tolower(c)
    if (.not.((c=="k").or.(c=="m").or.(c=="g"))) then
      call error("set_pes_mem: scale """//c//""" is not valid")
    end if
    pes_mem_scale = c

    if (isinteger(trim(adjustl(str(:i-1))))) then
      read(str(:i-1),*) pes_mem
    else
      call error("set_pes_mem: wrong argument """//trim(str)//"""")
    end if
  else
    pes_mem_scale = "0"

    if (isinteger(trim(adjustl(str)))) then
      read(str,*) pes_mem
    else
      call error("set_pes_mem: wrong argument """//trim(str)//"""")
    end if
  end if

  flag_pes_mem = .true.

end subroutine set_pes_mem

!====================================================================

subroutine init_pes_module()

  !--------------------------------------------------------
  ! init_pes_module requires init_images be called.
  !
  ! init_pes_module performs all the preliminary checks
  ! required to compute energy and forces by an external program.
  ! For this reason, all the other procedures inside this module
  ! don't need to check if their arguments are valid.
  !--------------------------------------------------------

  integer        :: err_n
  character(120) :: err_msg

  ! if I'm a slave, I go to... ----------------------------
#ifdef USE_MPI
  if (proc_id/=0) then
    call mmpi_init_pes_module()
    return
  end if
#endif

  ! else, if I'm the master... ----------------------------

  ! preliminary checks ------------------------------------
  if (flag_init_pes_module) then
    call error("init_pes_module: module pes already initialized")
  end if

  if (flag_init_images.eqv..false.) then
    call error("init_pes_module: module geometry not initialized")
  end if

  ! allocation section ------------------------------------
  allocate(pes_energy(0:image_n+1),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_pes_module: "//trim(err_msg))
  end if

  allocate(pes_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_pes_module: "//trim(err_msg))
  end if

  ! initialize pes parameters -----------------------------
  pes_energy            = 0.0_DBL
  pes_energy(0)         = start_energy
  pes_energy(image_n+1) = end_energy

  flag_init_pes_module = .true.

  ! Master finished its initialization. -------------------
  ! If slave processes exist,
  ! master initializes them.
#ifdef USE_MPI
  if (comm_sz>1) then
    call mmpi_init_pes_module()
  end if
#endif

end subroutine init_pes_module

!====================================================================

subroutine compute_pes_forces()

  !--------------------------------------------------------
  ! WARNING: compute_pes_forces *must* be called
  !          before init_tangents
  !--------------------------------------------------------

  integer   :: i
  integer   :: tid
  real(DBL) :: conv_threshold
  logical   :: converged

  ! use MPI variant? --------------------------------------
#ifdef USE_MPI
  if (comm_sz>1) then
    call mmpi_compute_pes_forces()
    return
  end if
#endif

  ! preliminary checks ------------------------------------
  if (flag_init_pes_module.eqv..false.) then
    call error("compute_pes_forces: module pes not initialized")
  end if

  ! working section ---------------------------------------
  tid = 0
  do i=1, image_n
    conv_threshold = get_scfconv()

    do
      call get_pes_forces(i,tid,conv_threshold,converged)

      if (converged) then
        exit
      else
        conv_threshold = conv_threshold*10.0_DBL
        write(FILEOUT,'(5X,"compute_pes_forces: process ",I3,&
          &": convergence threshold on image ",I3,&
          &" reduced to ",ES8.1)') tid,i,conv_threshold
      end if
    end do
  end do

end subroutine compute_pes_forces

!====================================================================

subroutine get_pes_forces(i,tid,conv_threshold,flag_conv,ig,pesf,pesg)

  !--------------------------------------------------------
  ! input:  an integer, the index of image_geom array
  !
  ! get_pes_forces prepares an input file, execute the external
  ! program, gets forces and energy from the latter's output,
  ! and removes all the intermediary files
  !--------------------------------------------------------
  
  integer,                           intent(IN)    :: i
  integer,                           intent(IN)    :: tid
  real(DBL),                         intent(IN)    :: conv_threshold
  logical,                           intent(INOUT) :: flag_conv
  real(DBL), dimension(:), optional, intent(IN)    :: ig
  real(DBL),               optional, intent(OUT)   :: pesf
  real(DBL), dimension(:), optional, intent(OUT)   :: pesg

  ! @end_user: add a new real(DBL) "max_programname_threshold" parameter
  real(DBL), parameter                             :: max_gaussian_threshold = 1.0E-5
  real(DBL), parameter                             :: max_siesta_threshold   = 1.0E-1
  integer,   parameter                             :: base_fnumb_in          = 1000
  integer,   parameter                             :: base_fnumb_out         = 2000
  integer,   parameter                             :: tid_lim                = base_fnumb_out-base_fnumb_in-1
  integer,   parameter                             :: opt_arg                = 3
  integer                                          :: fnumb_in
  integer                                          :: fnumb_out
  character(8)                                     :: tid_lim_str
  character(8)                                     :: real_str
  character(120)                                   :: fname_in
  character(120)                                   :: fname_out
  character(120)                                   :: dirname
  character(120)                                   :: auxdirname
  logical, dimension(opt_arg)                      :: arg_presence

  ! checking argument presence ----------------------------
  arg_presence = .false.

  if (present(ig)) then
    arg_presence(1) = .true.
  end if

  if (present(pesf)) then
    arg_presence(2) = .true.
  end if

  if (present(pesg)) then
    arg_presence(3) = .true.
  end if

  if (alltrue(arg_presence)) then
    if (flag_init_pes_module.eqv..false.) then
      call error("get_pes_forces: module pes not initialized")
    end if
  else if (.not.allfalse(arg_presence)) then
    call error("get_pes_forces: optional arguments were only partially supplied")
  end if

  ! check tid ---------------------------------------------
  if (tid>tid_lim) then
    write(tid_lim_str,'(I8)') tid_lim
    tid_lim_str=adjustl(tid_lim_str)
    call error("get_pes_forces: max processes allowed: "//trim(tid_lim_str))
  end if

  ! SCF Convergence threshold check -----------------------
  ! @end_user: add a new "case ("programname")"
  !            with the proper check for SCF convergence threshold
  !            yust before "case default"
  select case (pes_program)
  case ("gaussian")
    if (nint(log10(conv_threshold))>nint(log10(max_gaussian_threshold))) then
      write(real_str,'(ES8.1)') max_gaussian_threshold
      real_str=adjustl(real_str)
      call error("get_pes_forces: convergence threshold above"//&
        &" the default value of "//trim(real_str))
    end if
  case ("siesta")
    if (conv_threshold>max_siesta_threshold) then
      write(real_str,'(ES8.1)') max_siesta_threshold
      real_str = adjustl(real_str)
      call error("get_pes_forces: convergence threshold above"//&
        &" the default value of "//trim(real_str))
    end if
  case default
    call error("get_pes_forces: unknown program """//trim(pes_program)//"""")
  end select

  ! Preparing and executing the external program ----------
  fnumb_in  = base_fnumb_in + tid
  fnumb_out = base_fnumb_out + tid
  call set_dir(i,dirname,auxdirname)
  call f_chdir(dirname)

  ! @end_user: add a new "case ("programname")"
  !            with the calls to the 3 new subroutines
  !            yust before "case default"
  select case (pes_program)
  case ("gaussian")
    if (arg_presence(1)) then
      call write_gaussian_input(i,conv_threshold,fnumb_in,fname_in,fname_out,ig)
    else
      call write_gaussian_input(i,conv_threshold,fnumb_in,fname_in,fname_out)
    end if
    call exec_gaussian(fname_in,flag_conv)
    if (flag_conv) then
      if (arg_presence(1)) then
        call get_gaussian_output(i,fnumb_out,fname_out,pesf,pesg)
      else
        call get_gaussian_output(i,fnumb_out,fname_out)
      end if
    end if
  case ("siesta")
    call write_siesta_input(i,conv_threshold,fnumb_in,fname_in,fname_out)
    call exec_siesta(fnumb_in,fname_in,fname_out)
    call get_siesta_output(i,fnumb_out,fname_out,flag_conv)
  case default
    call error("get_pes_forces: unknown program """//trim(pes_program)//"""")
  end select

  call f_chdir("..")
  call get_auxiliary_files(dirname,auxdirname)
  call remove_dir(dirname)

end subroutine get_pes_forces

!====================================================================

subroutine update_pes_energy(arr)

  real(DBL), dimension(:), intent(IN) :: arr

  ! preliminary checks ------------------------------------
  if (flag_init_pes_module.eqv..false.) then
    call error("update_pes_energy: module pes not initialized")
  end if

  if (size(arr,1)/=image_n) then
    call error("update_pes_energy: wrong argument size")
  end if

  ! update ------------------------------------------------
  pes_energy(1:image_n) = arr

end subroutine update_pes_energy

!====================================================================

subroutine update_pes_forces(arr)

  real(DBL), dimension(:,:), intent(IN) :: arr

  ! preliminary checks ------------------------------------
  if (flag_init_pes_module.eqv..false.) then
    call error("update_pes_forces: module pes not initialized")
  end if

  if ((size(arr,1)/=image_n).or.(size(arr,2)/=geom_len)) then
    call error("update_pes_forces: wrong argument size")
  end if

  ! update ------------------------------------------------
  pes_forces = arr

end subroutine update_pes_forces

!====================================================================

real(DBL) function get_scfconv()

  !--------------------------------------------------------
  ! Returns the scf convergence threshold.
  ! If there is one specified by the user, that's the default
  ! value returned by this function. Otherwise, it returns
  ! a default value based on the pes program in use.
  !--------------------------------------------------------

  real(DBL), parameter :: gaussian_conv = 1.0E-8_DBL
  real(DBL), parameter :: siesta_conv   = 1.0E-4_DBL

  if (flag_pesd_scfconv) then
    get_scfconv = pesd_scfconv
  else
    if (pes_program=="gaussian") then
      get_scfconv = gaussian_conv
    else if (pes_program=="siesta") then
      get_scfconv = siesta_conv
    else
      if (flag_new_pes_program) then
        call error("get_scfconv: by specifying #NEWPESPROGRAM,"//&
          &" the #SCFCONV keyword becomes mandatory")
      else
        call error("get_scfconv: unknown program """//&
          &trim(pes_program)//"""")
      end if
    end if
  end if

end function get_scfconv

!====================================================================

integer function get_scfcycle()
  
  !--------------------------------------------------------
  ! Returns the scf max cycles.
  ! If there is one specified by the user, that's the default
  ! value returned by this function. Otherwise, it returns
  ! a default value based on the pes program in use.
  !--------------------------------------------------------

  integer, parameter :: gaussian_cycle = 64
  integer, parameter :: siesta_cycle   = 50

  if (flag_pesd_scfcycle) then
    get_scfcycle = pesd_scfcycle
  else
    if (pes_program=="gaussian") then
      get_scfcycle = gaussian_cycle
    else if (pes_program=="siesta") then
      get_scfcycle = siesta_cycle
    else
      if (flag_new_pes_program) then
        call error("get_scfcycle: by specifying #NEWPESPROGRAM,"//&
          &" the #SCFCYCLE keyword becomes mandatory")
      else
        call error("get_scfcycle: unknown program """//&
          &trim(pes_program)//"""")
      end if
    end if
  end if

end function get_scfcycle

!====================================================================
! Private MPI
!====================================================================

#ifdef USE_MPI
subroutine mmpi_init_pes_module()

  !--------------------------------------------------------
  ! Init all slave processes
  !--------------------------------------------------------

  integer        :: i
  character(8)   :: istr
  integer        :: cmd
  character(30)  :: str30
  character(120) :: str120
  logical        :: flag
  integer        :: err_n
  character(120) :: err_msg

  ! master gets slaves into this subroutine ---------------
  if (proc_id==0) then
    cmd = MMPI_MSG_INIT_PES_MODULE
    call mpi_bcast(cmd,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)
  end if

  if (proc_id==0) then ! master stuffs --------------------
    ! mandatory variables ---------------------------------
    ! bcast pes_program
    str30 = pes_program
    call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)

    ! bcast pes_exec
    str120 = pes_exec
    call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)

    ! optional variables ----------------------------------
    ! bcast pes_proc
    flag = flag_pes_proc
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str30,*) pes_proc
      str30 = adjustl(str30)
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast pes_mem
    flag = flag_pes_mem
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str30,'(I29,A1)') pes_mem, pes_mem_scale
      str30 = adjustl(str30)
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! pes_data variables ----------------------------------
    ! bcast pesd_scfcycle
    flag = flag_pesd_scfcycle
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str30,*) pesd_scfcycle
      str30 = adjustl(str30)
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast pesd_scfconv
    flag = flag_pesd_scfconv
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str30,*) pesd_scfconv
      str30 = adjustl(str30)
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast pesd_scfvshift
    flag = flag_pesd_scfvshift
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str30,*) pesd_scfvshift
      str30 = adjustl(str30)
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast pesd_intgrid
    flag = flag_pesd_intgrid
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str120,'(A)') pesd_intgrid
      str120 = adjustl(str120)
      call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast pesd_additional_cmd
    flag = flag_pesd_additional_cmd
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str120,'(A)') pesd_additional_cmd
      str120 = adjustl(str120)
      call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast auxiliary_input_files
    flag = flag_pesd_auxiliary_input_files
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(istr,'(I8)') pesd_auxiliary_input_files_n
      istr = adjustl(istr)
      str120 = "# "//trim(istr)

      do i=1, pesd_auxiliary_input_files_n
        str120 = trim(str120)//" "//trim(pesd_auxiliary_input_files(i))
      end do

      call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast auxiliary_output_files
    flag = flag_pesd_auxiliary_output_files
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(istr,'(I8)') pesd_auxiliary_output_files_n
      istr = adjustl(istr)
      str120 = "# "//trim(istr)

      do i=1, pesd_auxiliary_output_files_n
        str120 = trim(str120)//" "//trim(pesd_auxiliary_output_files(i))
      end do

      call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! pes_input_template variables ------------------------
    call mmpi_sync_pes_it()

!DEBUG: used to check if message passing was successful
!    call end_main_exec() ! MPI Debug

  else                 ! slaves stuffs --------------------
    ! preliminary checks ----------------------------------
    if (flag_init_pes_module) then
      write(istr,'(I8)') proc_id
      istr = adjustl(istr)
      err_msg = "mmpi_init_pes_module: process "//trim(istr)//&
        &": module pes already initialized"
      call error(err_msg)
    end if

    if (flag_init_images.eqv..false.) then
      write(istr,'(I8)') proc_id
      istr = adjustl(istr)
      err_msg = "mmpi_init_pes_module: process "//trim(istr)//&
        &": module geometry not initialized"
      call error(err_msg)
    end if

    ! mandatory variables ---------------------------------
    ! bcast pes_program
    call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    call set_pes_program(str30)

    ! bcast pes_exec
    call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    call set_pes_exec(str120)

    ! optional variables ----------------------------------
    ! bcast pes_proc
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pes_proc(str30)
    end if

    ! bcast pes_mem
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pes_mem(str30)
    end if

    ! pes_data variables ----------------------------------
    ! bcast pesd_scfcycle
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pesd_scfcycle(str30)
    end if

    ! bcast pesd_scfconv
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pesd_scfconv(str30)
    end if

    ! bcast pesd_scfvshift
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pesd_scfvshift(str30)
    end if

    ! bcast pesd_intgrid
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pesd_intgrid(str120)
    end if

    ! bcast pesd_additional_cmd
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pesd_additional_cmd(str120)
    end if

    ! bcast auxiliary_input_files
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pesd_auxiliary_files(PESD_AUX_INPUT_FILES,str120)
    end if

    ! bcast auxiliary_output_files
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str120,len(str120),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_pesd_auxiliary_files(PESD_AUX_OUTPUT_FILES,str120)
    end if

    ! pes_input_template variables ------------------------
    call mmpi_sync_pes_it()

    ! allocation section ------------------------------------
    allocate(pes_energy(0:image_n+1),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr = adjustl(istr)
      call error("mmpi_init_pes_module: process "//&
        &trim(istr)//": "//trim(err_msg))
    end if

    allocate(pes_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr = adjustl(istr)
      call error("mmpi_init_pes_module: process "//&
        &trim(istr)//": "//trim(err_msg))
    end if

    ! finalize --------------------------------------------
    flag_init_pes_module = .true.

!DEBUG: used to check if message passing was successful
!    write(*,*) proc_id," prog  ",trim(pes_program)
!    write(*,*) proc_id," flag  ",flag_pes_program
!    write(*,*) proc_id," exec  ",trim(pes_exec)
!    write(*,*) proc_id," flag  ",flag_pes_exec
!    write(*,*) proc_id," pproc ",pes_proc
!    write(*,*) proc_id," flag  ",flag_pes_proc
!    write(*,*) proc_id," pmem  ",pes_mem
!    write(*,*) proc_id," pmems ",pes_mem_scale
!    write(*,*) proc_id," flag  ",flag_pes_mem
!    write(*,*) proc_id,"-------------------"
!    write(*,*) proc_id," pesd_scfcycle       ",pesd_scfcycle
!    write(*,*) proc_id," pesd_scfconv        ",pesd_scfconv
!    write(*,*) proc_id," pesd_scfvshift      ",pesd_scfvshift
!    write(*,*) proc_id," pesd_intgrid        ",pesd_intgrid
!    write(*,*) proc_id," pesd_additional_cmd ",pesd_additional_cmd

  end if

end subroutine mmpi_init_pes_module
#endif

!====================================================================

#ifdef USE_MPI
subroutine mmpi_compute_pes_forces()

  !--------------------------------------------------------
  ! MPI version of compute_pes_forces subroutine.
  ! This one is executed automatically if MPI is used
  ! and executing processes are more than one.
  !--------------------------------------------------------

  character(8)                           :: istr
  integer                                :: i
  integer                                :: total_computations
  integer                                :: cmd
  real(DBL)                              :: real_buff
  real(DBL), allocatable, dimension(:)   :: real_r1_buff
  real(DBL), allocatable, dimension(:)   :: energies_buff
  real(DBL), allocatable, dimension(:,:) :: forces_buff
  real(DBL), allocatable, dimension(:,:) :: image_geom_buff
  real(DBL)                              :: conv_threshold
  logical                                :: converged
  integer, dimension(MPI_STATUS_SIZE)    :: mstatus
  integer                                :: err_n
  character(120)                         :: err_msg

  ! master gets slaves into this subroutine ---------------
  if (proc_id==0) then
    cmd = MMPI_MSG_COMPUTE_PES_FORCES
    call mpi_bcast(cmd,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)
  end if

  !--------------------------------------------------------
  ! both master and slaves
  !--------------------------------------------------------

  ! preliminary checks ------------------------------------
  if (flag_init_pes_module.eqv..false.) then
    write(istr,'(I8)') proc_id
    istr = adjustl(istr)
    err_msg = "mmpi_compute_pes_forces: process "//trim(istr)//&
      &": module pes not initialized"
    call error(err_msg)
  end if

  ! allocation section ------------------------------------
  allocate (image_geom_buff(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    write(istr,'(I8)') proc_id
    istr = adjustl(istr)
    call error("mmpi_compute_pes_forces: process "//&
      &trim(istr)//": "//trim(err_msg))
  end if

  allocate (real_r1_buff(geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    write(istr,'(I8)') proc_id
    istr = adjustl(istr)
    call error("mmpi_compute_pes_forces: process "//&
      &trim(istr)//": "//trim(err_msg))
  end if

  if (proc_id==0) then
    allocate (energies_buff(image_n),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr = adjustl(istr)
      call error("mmpi_compute_pes_forces: process "//&
        &trim(istr)//": "//trim(err_msg))
    end if

    allocate (forces_buff(image_n,geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr = adjustl(istr)
      call error("mmpi_compute_pes_forces: process "//&
        &trim(istr)//": "//trim(err_msg))
    end if
  end if

  ! bcast image_geom --------------------------------------
  if (proc_id==0) then
    image_geom_buff = image_geom(1:image_n,:)
    call mpi_bcast(image_geom_buff,image_n*geom_len,MPI_REAL8,0,MPI_COMM_WORLD,err_n)
  else
    call mpi_bcast(image_geom_buff,image_n*geom_len,MPI_REAL8,0,MPI_COMM_WORLD,err_n)
    call update_images(image_geom_buff)
  end if

  ! working section ---------------------------------------
  total_computations = 0

  do i=1, image_n
    if (proc_id==mod(i,comm_sz)) then
      total_computations = total_computations+1
      conv_threshold = get_scfconv()

      do
        call get_pes_forces(i,proc_id,conv_threshold,converged)

        if (converged) then
          exit
        else
          conv_threshold = conv_threshold*10.0_DBL
          write(FILEOUT,'(5X,"compute_pes_forces: process ",I3,&
            &": convergence threshold on image ",I3,&
            &" reduced to ",ES8.1)') proc_id,i,conv_threshold
        end if
      end do
    end if
  end do

  ! collect results ---------------------------------------
  if (proc_id==0) then ! master stuffs --------------------

    ! write master results into buffers
    energies_buff = 0.0_DBL
    forces_buff   = 0.0_DBL

    do i=1, image_n 
      if (proc_id==mod(i,comm_sz)) then
        energies_buff(i) = pes_energy(i)
        forces_buff(i,:) = pes_forces(i,:)
      end if
    end do

    ! get energies from slaves
    do i=1, image_n-total_computations
      call mpi_recv(real_buff,1,MPI_REAL8,MPI_ANY_SOURCE,&
        &MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,err_n)
      energies_buff(mstatus(MPI_TAG)) = real_buff
    end do

    ! all energies received -------------------------------
    ! meet the slaves at the barrier
    call mpi_barrier(MPI_COMM_WORLD,err_n)

    ! get forces from slaves
    do i=1, image_n-total_computations
      call mpi_recv(real_r1_buff,size(real_r1_buff,1),MPI_REAL8,&
        &MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mstatus,err_n)
      forces_buff(mstatus(MPI_TAG),:) = real_r1_buff
    end do

    ! update results
    call update_pes_energy(energies_buff)
    call update_pes_forces(forces_buff)

  else                 ! slaves stuffs --------------------

    ! send energies to master
    do i=1, image_n
      if (proc_id==mod(i,comm_sz)) then
        real_buff = pes_energy(i)
        call mpi_send(real_buff,1,MPI_REAL8,0,i,MPI_COMM_WORLD,err_n)
      end if
    end do

    ! all energies sended ---------------------------------
    ! meet the master and other slaves at the barrier
    call mpi_barrier(MPI_COMM_WORLD,err_n)

    ! send forces to master
    do i=1, image_n
      if (proc_id==mod(i,comm_sz)) then
        real_r1_buff = pes_forces(i,:)
        call mpi_send(real_r1_buff,size(real_r1_buff,1),MPI_REAL8,&
          &0,i,MPI_COMM_WORLD,err_n)
      end if
    end do

  end if

  ! deallocation section ----------------------------------
  deallocate (image_geom_buff,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    write(istr,'(I8)') proc_id
    istr = adjustl(istr)
    call error("mmpi_compute_pes_forces: process "//&
      &trim(istr)//": "//trim(err_msg))
  end if

  deallocate (real_r1_buff,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    write(istr,'(I8)') proc_id
    istr = adjustl(istr)
    call error("mmpi_compute_pes_forces: process "//&
      &trim(istr)//": "//trim(err_msg))
  end if

  if (proc_id==0) then
    deallocate (energies_buff,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr = adjustl(istr)
      call error("mmpi_compute_pes_forces: process "//&
        &trim(istr)//": "//trim(err_msg))
    end if

    deallocate (forces_buff,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr = adjustl(istr)
      call error("mmpi_compute_pes_forces: process "//&
        &trim(istr)//": "//trim(err_msg))
    end if
  end if

end subroutine mmpi_compute_pes_forces
#endif

!====================================================================
! Private
!====================================================================

!====================================================================
! Generic subroutines for get_pes_forces
!====================================================================

subroutine set_dir(i,dirname,auxdirname)

  !--------------------------------------------------------
  ! Makes the working directory
  ! and copies in it the auxiliary files.
  !--------------------------------------------------------

  integer,      intent(IN)    :: i
  character(*), intent(INOUT) :: dirname
  character(*), intent(INOUT) :: auxdirname

  integer, save               :: calling_count = 0
  character(8)                :: i_str
  character(8)                :: exit_n_str
  character(300)              :: cmd
  integer                     :: j
  integer                     :: exit_n
  integer                     :: cmd_n

  ! preliminary settings ----------------------------------
  if (calling_count<=image_n) then
    calling_count = calling_count+1
  end if

  write(i_str,'(I8)') i
  i_str      = adjustl(i_str)
  dirname    = base_dirname//trim(i_str)
  auxdirname = base_auxdirname//trim(i_str)

  ! make directories for storing auxiliary output files ---
  if ((flag_pesd_auxiliary_output_files).and.(calling_count<=image_n)) then
    cmd = "mkdir "//trim(auxdirname)
    call execute_command_line(trim(cmd),&
      &wait=.true.,exitstat=exit_n,cmdstat=cmd_n)

    if (cmd_n/=0) then
      call error("set_dir: cannot execute command """//trim(cmd)//"""")
    end if

    if (exit_n/=0) then
      write(exit_n_str,'(I8)') exit_n
      exit_n_str = adjustl(exit_n_str)
      call error("set_dir: """//trim(cmd)//&
        &""" terminated with exit code: "//trim(exit_n_str))
    end if
  end if

  ! make the working directory ----------------------------
  cmd = "mkdir "//trim(dirname)
  call execute_command_line(trim(cmd),&
    &wait=.true.,exitstat=exit_n,cmdstat=cmd_n)

  if (cmd_n/=0) then
    call error("set_dir: cannot execute command """//trim(cmd)//"""")
  end if

  if (exit_n/=0) then
    write(exit_n_str,'(I8)') exit_n
    exit_n_str = adjustl(exit_n_str)
    call error("set_dir: """//trim(cmd)//&
      &""" terminated with exit code: "//trim(exit_n_str))
  end if

  ! copy the auxiliary input files ------------------------
  if (flag_pesd_auxiliary_input_files) then
    cmd = "cp"
    do j=1, pesd_auxiliary_input_files_n
      cmd = trim(cmd)//" "//trim(pesd_auxiliary_input_files(j))
    end do
    cmd = trim(cmd)//" "//trim(dirname)//"/."

    call execute_command_line(trim(cmd),&
      &wait=.true.,exitstat=exit_n,cmdstat=cmd_n)

    if (cmd_n/=0) then
      call error("set_dir: cannot execute command """//trim(cmd)//"""")
    end if

    if (exit_n/=0) then
      write(exit_n_str,'(I8)') exit_n
      exit_n_str = adjustl(exit_n_str)
      call error("set_dir: """//trim(cmd)//&
        &""" terminated with exit code: "//trim(exit_n_str))
    end if
  end if

  ! copy the auxiliary output files ------------------------
  if ((flag_pesd_auxiliary_output_files).and.(calling_count>image_n)) then
    cmd = "cp"
    do j=1, pesd_auxiliary_output_files_n
      cmd = trim(cmd)//" "//trim(auxdirname)//&
        &"/"//trim(pesd_auxiliary_output_files(j))
    end do
    cmd = trim(cmd)//" "//trim(dirname)//"/."

    call execute_command_line(trim(cmd),&
      &wait=.true.,exitstat=exit_n,cmdstat=cmd_n)

    if (cmd_n/=0) then
      call error("set_dir: cannot execute command """//trim(cmd)//"""")
    end if

    if (exit_n/=0) then
      write(exit_n_str,'(I8)') exit_n
      exit_n_str = adjustl(exit_n_str)
      call error("set_dir: """//trim(cmd)//&
        &""" terminated with exit code: "//trim(exit_n_str))
    end if
  end if

end subroutine set_dir

!====================================================================

subroutine get_auxiliary_files(dirname,auxdirname)

  character(*), intent(IN) :: dirname
  character(*), intent(IN) :: auxdirname

  integer                  :: j
  character(300)           :: cmd
  character(8)             :: exit_n_str
  integer                  :: exit_n
  integer                  :: cmd_n

  if (flag_pesd_auxiliary_output_files) then
    cmd = "cp"
    do j=1, pesd_auxiliary_output_files_n
      cmd = trim(cmd)//" "//trim(dirname)//&
        &"/"//trim(pesd_auxiliary_output_files(j))
    end do
    cmd = trim(cmd)//" "//trim(auxdirname)//"/."

    call execute_command_line(trim(cmd),&
      &wait=.true.,exitstat=exit_n,cmdstat=cmd_n)

    if (cmd_n/=0) then
      call error("get_auxiliary_files: cannot execute command """//trim(cmd)//"""")
    end if

    if (exit_n/=0) then
      write(exit_n_str,'(I8)') exit_n
      exit_n_str = adjustl(exit_n_str)
      call error("get_auxiliary_files: """//trim(cmd)//&
        &""" terminated with exit code: "//trim(exit_n_str))
    end if
  end if

end subroutine get_auxiliary_files

!====================================================================

subroutine remove_dir(dirname)

  character(*), intent(IN) :: dirname

  character(8)             :: exit_n_str
  character(140)           :: cmd
  integer                  :: exit_n
  integer                  :: cmd_n

  cmd = "rm -r "//trim(dirname)
  call execute_command_line(trim(cmd),&
    &wait=.true.,exitstat=exit_n,cmdstat=cmd_n)

  if (cmd_n/=0) then
    call error("remove_dir: cannot execute command """//trim(cmd)//"""")
  end if

  if (exit_n/=0) then
    write(exit_n_str,'(I8)') exit_n
    exit_n_str = adjustl(exit_n_str)
    call error("remove_dir: """//trim(cmd)//&
      &""" terminated with exit code: "//trim(exit_n_str))
  end if

end subroutine remove_dir

!====================================================================
! Gaussian Section
!====================================================================

subroutine write_gaussian_input(i,conv_threshold,fnumb_in,fname_in,fname_out,ig)

  integer,                           intent(IN)  :: i
  real(DBL),                         intent(IN)  :: conv_threshold
  integer,                           intent(IN)  :: fnumb_in
  character(*),                      intent(OUT) :: fname_in
  character(*),                      intent(OUT) :: fname_out
  real(DBL), dimension(:), optional, intent(IN)  :: ig

  integer                                        :: j
  integer                                        :: conv_val
  integer                                        :: indx
  character(30)                                  :: tmp_str
  character(8)                                   :: i_str
  character(120)                                 :: err_msg
  integer                                        :: err_n

  ! preliminary checks ------------------------------------
  if (present(ig)) then
    if (size(ig)/=geom_len) then
      call error("write_gaussian_input: wrong ig argument size")
    end if
  end if

  ! set fname_in and fname_out ----------------------------
  write(i_str,'(I8)') i
  i_str     = adjustl(i_str)
  fname_in  = base_name//trim(i_str)//".com"
  fname_out = base_name//trim(i_str)//".log"

  conv_val  = abs(nint(log10(conv_threshold)))

  ! open unit ---------------------------------------------
  open(unit=fnumb_in,file=fname_in,status='NEW',action='WRITE',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("write_gaussian_input: "//trim(err_msg))
  end if

  ! write stuffs ------------------------------------------
  ! #PESINPUTTEMPLATE 1
  call write_pes_it(fnumb_in,1)

  ! write SCF threshold
  write(tmp_str,'(I30)') abs(nint(log10(conv_threshold)))
  tmp_str = adjustl(tmp_str)
  write(fnumb_in,'(A,A,A)') "#scf(conver=",trim(tmp_str),")"

  ! write SCF cycles
  write(tmp_str,'(I30)') get_scfcycle()
  tmp_str = adjustl(tmp_str)
  write(fnumb_in,'(A,A,A)') "#scf(maxcycle=",trim(tmp_str),")"

  ! write forces
  write(fnumb_in,'(A)') "#force test"
  write(fnumb_in,*)

  ! #PESINPUTTEMPLATE 2
  call write_pes_it(fnumb_in,2)

  ! geometry
  do j=1, geom_len
    if (mod(j,3)==1) then
      write(fnumb_in,'(A3,1X)',advance='no') element(j/3+1)
    end if
    
    if (present(ig)) then
      write(fnumb_in,'(1X,F13.6)',advance='no') ig(j)
    else
      write(fnumb_in,'(1X,F13.6)',advance='no') image_geom(i,j)
    end if

    if (mod(j,3)==0) then
      write(fnumb_in,*)
    end if
  end do
  write(fnumb_in,*)

  ! #PESINPUTTEMPLATE 3 (optional)
  indx = get_pes_it_n(3)
  if (indx/=0) then
    call write_pes_it(fnumb_in,3)
    write(fnumb_in,*)
  end if

  ! close unit --------------------------------------------
  close(unit=fnumb_in,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("write_gaussian_input: "//trim(err_msg))
  end if

end subroutine write_gaussian_input

!====================================================================

subroutine exec_gaussian(fname_in,flag_conv)

  character(*), intent(IN)  :: fname_in
  logical,      intent(OUT) :: flag_conv

  integer, parameter        :: consecutive_failures_lim = 100
  integer, save             :: consecutive_failures     = 0
  character(140)            :: cmd
  integer                   :: exit_n
  integer                   :: cmd_n
  character(8)              :: i_str

  cmd = trim(pes_exec)//" "//fname_in(:len_trim(fname_in)-4)
  call execute_command_line(trim(cmd),&
    &wait=.true.,exitstat=exit_n,cmdstat=cmd_n)

  if (cmd_n/=0) then
    call error("exec_gaussian: cannot execute command """//trim(cmd)//"""")
  end if
  
  ! if gaussian does not converge, it returns an exit status /= 0
  if (exit_n/=0) then
    consecutive_failures = consecutive_failures+1
    if (consecutive_failures >= consecutive_failures_lim) then
      write(i_str,'(I8)') consecutive_failures
      i_str = adjustl(i_str)
      call error("exec_gaussian: Gaussian has failed consecutively "//trim(i_str)//" times")
    end if

    flag_conv = .false.
  else
    consecutive_failures = 0

    flag_conv = .true.
  end if

end subroutine exec_gaussian

!====================================================================

subroutine get_gaussian_output(i,fnumb_out,fname_out,pesf,pesg)

  integer,                           intent(IN)  :: i
  integer,                           intent(IN)  :: fnumb_out
  character(*),                      intent(IN)  :: fname_out
  real(DBL),               optional, intent(OUT) :: pesf
  real(DBL), dimension(:), optional, intent(OUT) :: pesg

  character(*), parameter                        :: norm_term = " Normal termination of Gaussian"
  integer,      parameter                        :: opt_arg   = 2
  integer                                        :: j, k
  character(200)                                 :: str, field
  logical, dimension(opt_arg)                    :: arg_presence
  integer                                        :: err_n
  character(120)                                 :: err_msg

  ! Checking optional argument ----------------------------
  arg_presence = .false.

  if (present(pesf)) then
    arg_presence(1) = .true.
  end if

  if (present(pesg)) then
    if (size(pesg)/=geom_len) then
      call error("get_gaussian_output: wrong pesg argument size")
    end if
    arg_presence(2) = .true.
  end if

  if (.not.(alltrue(arg_presence).or.&
            &allfalse(arg_presence))) then
    call error("get_gaussian_output: optional arguments were only partially supplied")
  end if

  ! Open output file, read from end -----------------------
  open(unit=fnumb_out,file=fname_out,status='OLD',action='READ',&
    &iostat=err_n,iomsg=err_msg,position='APPEND')
  if (err_n/=0) then
    call error("get_gaussian_output: "//trim(err_msg))
  end if
  backspace(fnumb_out,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("get_gaussian_output: "//trim(err_msg))
  end if
 
  ! Check for normal termination of Gaussian --------------
  read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
  if (err_n/=0) then
    call error("get_gaussian_output: "//trim(err_msg))
  else if (str(:31)/=norm_term) then
    call error("get_gaussian_output: bad termination of "//fname_out)
  end if

  ! Rewind file and check forces and energy ---------------
  rewind(unit=fnumb_out,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("get_gaussian_output: "//trim(err_msg))
  end if

  do
    read(fnumb_out,'(A200)',iostat=err_n) str

    ! Check end of file -----------------------------------
    if (err_n/=0) then
      exit
    end if

    ! Parse str -------------------------------------------
    if (str(:10)==" SCF Done:") then
      call get_field(str,field,5,err_n,err_msg)
      if (err_n/=0) then
        call error("get_gaussian_output: "//trim(err_msg))
      end if

      if (arg_presence(1).eqv..true.) then
        read(field,*,iostat=err_n,iomsg=err_msg) pesf
      else
        read(field,*,iostat=err_n,iomsg=err_msg) pes_energy(i)
      end if

      if (err_n/=0) then
        call error("get_gaussian_output: "//trim(err_msg))
      end if
    else if (str==" ***** Axes restored to original set *****") then
      do j=1, 4
        read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
        if (err_n/=0) then
          call error("get_gaussian_output: "//trim(err_msg))
        end if
      end do

      do j=1, geom_len/3
        read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
        if (err_n/=0) then
          call error("get_gaussian_output: "//trim(err_msg))
        end if
        
        do k=3,5
          call get_field(str,field,k,err_n,err_msg)
          if (err_n/=0) then
            call error("get_gaussian_output: "//trim(err_msg))
          end if

          if (arg_presence(1).eqv..true.) then
            read(field,*,iostat=err_n,iomsg=err_msg) pesg(3*(j-1)+(k-2))
          else
            read(field,*,iostat=err_n,iomsg=err_msg) pes_forces(i,3*(j-1)+(k-2))
          end if

          if (err_n/=0) then
            call error("get_gaussian_output: "//trim(err_msg))
          end if
        end do
      end do
    end if
  end do

  ! Convert forces from Hartree/Bohr to Hartree/Ang -------
  if (arg_presence(1).eqv..true.) then
    pesg = pesg*(1.0_DBL/BOHR_ON_ANG)
  else
    pes_forces(i,:) = pes_forces(i,:)*(1.0_DBL/BOHR_ON_ANG)
  end if
 
  ! Close output file -------------------------------------
  close(unit=fnumb_out,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("get_gaussian_output: "//trim(err_msg))
  end if

end subroutine get_gaussian_output

!====================================================================
! Siesta Section
!====================================================================

subroutine write_siesta_input(i,conv_threshold,fnumb_in,fname_in,fname_out,ig)

  integer,                           intent(IN)  :: i
  real(DBL),                         intent(IN)  :: conv_threshold
  integer,                           intent(IN)  :: fnumb_in
  character(*),                      intent(OUT) :: fname_in
  character(*),                      intent(OUT) :: fname_out
  real(DBL), dimension(:), optional, intent(IN)  :: ig

  character(8)                                   :: i_str
  integer                                        :: j
  integer                                        :: err_n
  character(120)                                 :: err_msg

  ! set fname_in and fname_out ----------------------------
  write(i_str,'(I8)') i
  i_str     = adjustl(i_str)
  fname_in  = base_name//trim(i_str)//".fdf"
  fname_out = base_name//trim(i_str)//".out"

  ! open unit ---------------------------------------------
  open(unit=fnumb_in,file=fname_in,status='NEW',action='WRITE',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("write_siesta_input: "//trim(err_msg))
  end if

  ! write stuffs ------------------------------------------
  ! #PESINPUTTEMPLATE 1
  call write_pes_it(fnumb_in,1)

  ! geometry
  do j=1, geom_len
    if (present(ig)) then
      write(fnumb_in,'(1X,F13.6)',advance='no') ig(j)
    else
      write(fnumb_in,'(1X,F13.6)',advance='no') image_geom(i,j)
    end if

    if (mod(j,3)==0) then
      write(fnumb_in,'(1X,A3)') elabel(j/3)
    end if
  end do

  ! #PESINPUTTEMPLATE 2
  call write_pes_it(fnumb_in,2)

  ! write SCF threshold
  write(fnumb_in,'("DM.Tolerance ",ES9.2)') conv_threshold
  write(fnumb_in,*)

  ! write SCF cycles
  write(fnumb_in,'("MaxSCFIterations ",I8)') get_scfcycle()
  write(fnumb_in,*)

  ! write forces
  write(fnumb_in,'("WriteForces true")')
  write(fnumb_in,*)

  ! close unit --------------------------------------------
  close(unit=fnumb_in,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("write_siesta_input: "//trim(err_msg))
  end if

end subroutine write_siesta_input

!====================================================================

subroutine exec_siesta(fnumb_in,fname_in,fname_out)

  integer,      intent(IN) :: fnumb_in
  character(*), intent(IN) :: fname_in
  character(*), intent(IN) :: fname_out

  character(*),  parameter :: my_name     = "exec_siesta"
  character(*),  parameter :: script_name = "siesta_omp.sh"
  character(140)           :: cmd
  integer                  :: exit_n
  integer                  :: cmd_n
  character(8)             :: exit_n_str
  character(8)             :: i_str
  integer                  :: err_n
  character(120)           :: err_msg

  ! determine how to execute siesta -----------------------
  if (pes_proc > 1) then ! -----------> ! run siesta in parallel
    if (flag_pes_program_with_mpi) then ! run siesta via MPI
      write(i_str,'(I8)') pes_proc
      i_str = adjustl(i_str)
      cmd = "mpirun -n "//trim(i_str)//" "//trim(pes_exec)//" < "//trim(fname_in)//" > "//trim(fname_out)
    else ! ---------------------------> ! run siesta via openMP
      ! open the script file
      open(unit=fnumb_in,file=script_name,status='NEW',action='WRITE',&
        &iostat=err_n,iomsg=err_msg,position='REWIND')
      if (err_n/=0) then
        call error(my_name//": "//trim(err_msg))
      end if

      ! write the script
      write(i_str,'(I8)') pes_proc
      i_str = adjustl(i_str)

      write(fnumb_in,'(A)') "#!/bin/bash"
      write(fnumb_in,'(A)') "export OMP_NUM_THREADS="//trim(i_str)
      write(fnumb_in,'(A)') trim(pes_exec)//" < "//trim(fname_in)//" > "//trim(fname_out)
      write(fnumb_in,'(A)') "exit_code=$?"
      write(fnumb_in,'(A)') "exit $exit_code"

      ! close the script file
      close(unit=fnumb_in,iostat=err_n,iomsg=err_msg)
      if (err_n/=0) then
        call error(my_name//": "//trim(err_msg))
      end if

      ! compose cmd
      cmd = "bash "//trim(script_name)
    end if
  else ! -----------------------------> ! run siesta in serial
    cmd = trim(pes_exec)//" < "//trim(fname_in)//" > "//trim(fname_out)
  end if

  ! execute siesta ----------------------------------------
  call execute_command_line(trim(cmd), wait=.true., exitstat=exit_n, cmdstat=cmd_n)

  if (cmd_n/=0) then
    call error(my_name//": cannot execute command """//trim(cmd)//"""")
  end if
  
  ! on error, siesta returns an exit status /= 0
  if (exit_n/=0) then
    write(exit_n_str,'(I8)') exit_n
    exit_n_str = adjustl(exit_n_str)
    call error(my_name//": "//trim(pes_exec)//" terminated with exit code: "//trim(exit_n_str))
  end if

end subroutine exec_siesta

!====================================================================

subroutine get_siesta_output(i,fnumb_out,fname_out,flag_conv)

  integer,      intent(IN)  :: i
  integer,      intent(IN)  :: fnumb_out
  character(*), intent(IN)  :: fname_out
  logical,      intent(OUT) :: flag_conv

  integer                   :: j
  integer                   :: k
  integer                   :: k_start
  integer                   :: default_scfcycle ! max scf cycles for the siesta computation
  integer                   :: read_scfcycle    ! last scf cycle read from siesta output
  logical                   :: correct_one
  character(200)            :: field
  character(200)            :: str
  integer                   :: err_n
  character(120)            :: err_msg

  ! variables initialization ------------------------------
  flag_conv        = .true.
  correct_one      = .false.
  read_scfcycle    = -1
  default_scfcycle = get_scfcycle()
  default_scfcycle = default_scfcycle - 1 ! If we tell siesta
    ! to do N scf cycles, we will find at most N - 1 scf cycles
    ! in the output file

  ! open unit ---------------------------------------------
  open(unit=fnumb_out,file=fname_out,status='OLD',action='READ',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("get_siesta_output: "//trim(err_msg))
  end if

  ! search for convergence --------------------------------
  ! read until the string "siesta: iscf" or "        iscf" is found
  do
    read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
    if (err_n/=0) then
      write(FILEOUT,*) "WAR get_siesta_output: ",&
        &"string ""siesta: iscf"" or ""        iscf"" not found"
      call error("get_siesta_output: "//trim(err_msg))
    end if

    if ((str(1:12)=="siesta: iscf").or.(str(1:12)=="        iscf")) then
      exit
    end if
  end do

  ! get the last scf iteration number
  do
    read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
    if (err_n/=0) then
      write(FILEOUT,*) "WAR get_siesta_output: ",&
        &"last scf iteration not found"
      call error("get_siesta_output: "//trim(err_msg))
    end if

    if (str=="") then
      if (read_scfcycle==-1) then
        call error("get_siesta_output: scf iterations not found")
      else if (read_scfcycle>=default_scfcycle) then
        flag_conv = .false.
      else
        flag_conv = .true.
      end if
      
      exit
    end if

    call get_field(str,field,2,err_n,err_msg)
    if (err_n/=0) then
      call error("get_siesta_output: "//trim(err_msg))
    end if

    if (isinteger(trim(adjustl(field)))) then
      read(field,*) read_scfcycle
    end if
  end do

  ! if converged, get energy and forces -------------------
  if (flag_conv) then
    ! read until the string "siesta: E_KS(eV)" is found,
    ! then get the total energy value
    do
      read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
      if (err_n/=0) then
        write(FILEOUT,*) "WAR get_siesta_output: ",&
          &"string ""siesta: E_KS(eV)"" not found"
        call error("get_siesta_output: "//trim(err_msg))
      end if

      if (str(1:16)=="siesta: E_KS(eV)") then
        call get_field(str,field,4,err_n,err_msg)
        if (err_n/=0) then
          call error("get_siesta_output: "//trim(err_msg))
        end if

        if (isreal(trim(adjustl(field)))) then
          read(field,*) pes_energy(i)
        else
          call error("get_siesta_output: wrong field """//&
            &trim(field)//""" while reading total energy")
        end if

        exit
      end if
    end do

    ! read until the string "siesta: Atomic forces" is found
    do
      read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
      if (err_n/=0) then
        write(FILEOUT,*) "WAR get_siesta_output: ",&
          &"string ""siesta: Atomic forces"" not found"
        call error("get_siesta_output: "//trim(err_msg))
      end if

      if (str(1:21)=="siesta: Atomic forces") then
        ! check if it's the correct one
        read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
        if (err_n/=0) then
          write(FILEOUT,*) "WAR get_siesta_output: ",&
            &"string ""siesta: Atomic forces"" not found"
          call error("get_siesta_output: "//trim(err_msg))
        end if

        call get_field(str,field,1,err_n,err_msg)
        if (err_n/=0) then
          call error("get_siesta_output: "//trim(err_msg))
        end if

        if (field=="siesta:") then
          k_start     = 3
          correct_one = .true.
        else if (isinteger(trim(adjustl(field)))) then
          k_start     = 2
          correct_one = .true.
        end if

        if (correct_one) then
          backspace(unit=fnumb_out,iostat=err_n,iomsg=err_msg)
          if (err_n/=0) then
            write(FILEOUT,*) "WAR get_siesta_output: ",&
              &"error while backspacing"
            call error("get_siesta_output: "//trim(err_msg))
          end if

          ! continue to get forces values
          exit
        end if
      end if
    end do

    ! get forces values
    do j=1, geom_len/3
      read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
      if (err_n/=0) then
        write(FILEOUT,*) "WAR get_siesta_output: ",&
          &"error while reading atomic forces"
        call error("get_siesta_output: "//trim(err_msg))
      end if

      do k=k_start, k_start+2
        call get_field(str,field,k,err_n,err_msg)
        if (err_n/=0) then
          call error("get_siesta_output: "//trim(err_msg))
        end if

        if (isreal(trim(adjustl(field)))) then
          read(field,*) pes_forces(i,3*(j-1)+(k-(k_start-1)))
        else
          call error("get_siesta_output: wrong field """//&
            &trim(field)//""" while reading atomic forces")
        end if
      end do
    end do
  end if

  ! Convert forces from eV/Ang to Hartree/Ang -------------
  pes_forces(i,:) = pes_forces(i,:)*(1.0_DBL/AU_ON_EV)

  ! close unit --------------------------------------------
  close(unit=fnumb_out,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("get_siesta_output: "//trim(err_msg))
  end if
 
end subroutine get_siesta_output

!====================================================================

end module pes

