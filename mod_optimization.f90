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

module optimization

  use utility
  use bfgs
  use bfgs_wrapper
  use idpp
  use geometry
  use pes
  use elastic
  use climbing
  use output

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: set_optmz_algo,   &
            set_optmz_nsteps, &
            set_optmz_tol,    &
            optmz_pes,        &
            optmz_idpp

  !--------------------------------------------------------
  ! ENUM
  integer, parameter :: ALGO_SD           = 0
  integer, parameter :: ALGO_BFGS         = 1
  integer, parameter :: ALGO_LBFGS        = 2
  integer, parameter :: ALGO_FIRE         = 3
  integer            :: optmz_algo        = ALGO_FIRE ! default algorithm

  logical            :: flag_optmz_nsteps = .false.
  logical            :: flag_optmz_tol    = .false.
  integer            :: optmz_nsteps
  real(DBL)          :: optmz_tol

contains

!====================================================================
! Public
!====================================================================

subroutine set_optmz_algo(str)

  ! Set the optimization algorithm.
  ! It must be called at most once.

  character(*), intent(INOUT) :: str

  logical, save               :: first_call = .true.

  if (first_call.eqv..false.) then
    call error("set_optmz_algo: subroutine called more than once")
  end if

  call tolower(str)

  select case (str)
  case ("sd","steepest_descent")
    optmz_algo = ALGO_SD
  case ("bfgs")
    optmz_algo = ALGO_BFGS
  case ("lbfgs")
    optmz_algo = ALGO_LBFGS
  case ("fire")
    optmz_algo = ALGO_FIRE
  case default
    call error("set_optmz_algo: unknown optimization algorithm """//&
      &trim(str)//"""")
  end select

  first_call = .false.

end subroutine set_optmz_algo

!====================================================================

subroutine set_optmz_nsteps(str)
  
  character(*), intent(IN) :: str

  if (flag_optmz_nsteps) then
    call error("set_optmz_nsteps: optimization steps already setted")
  end if

  if (isinteger(trim(adjustl(str)))) then
    read(str,*) optmz_nsteps
  else
    call error("set_optmz_nsteps: argument must be an integer")
  end if

  if (optmz_nsteps<0) then
    optmz_nsteps = -1
  end if

  flag_optmz_nsteps = .true.

end subroutine set_optmz_nsteps

!====================================================================

subroutine set_optmz_tol(str)
  
  character(*), intent(IN) :: str

  if (flag_optmz_tol) then
    call error("set_optmz_tol: optimization tolerance already setted")
  end if

  if (isreal(trim(adjustl(str)))) then
    read(str,*) optmz_tol
  else
    call error("set_optmz_tol: argument must be a real")
  end if

  flag_optmz_tol = .true.

end subroutine set_optmz_tol

!====================================================================

subroutine optmz_pes(flag_out)

  logical, intent(OUT) :: flag_out ! true if convergent, false otherwise

  character(8)         :: istr

  select case (optmz_algo)
  case (ALGO_SD)
    call optmz_steepest_descent(PES_MODE,flag_out)
  case (ALGO_BFGS)
    call optmz_bfgs(PES_MODE,flag_out)
  case (ALGO_LBFGS)
    call optmz_lbfgs(PES_MODE,flag_out)
  case (ALGO_FIRE)
    call optmz_fire(PES_MODE,flag_out)
  case default
    write(istr,'(I8)') optmz_algo
    istr = adjustl(istr)
    call error("optmz_pes: wrong optimization algorithm """//trim(istr)//"""")
  end select

end subroutine optmz_pes

!====================================================================

subroutine optmz_idpp(flag_out)

  logical, intent(OUT) :: flag_out ! true if convergent, false otherwise

  call optmz_fire(IDPP_MODE,flag_out,nsteps=1000,tol=idpp_tol)

end subroutine optmz_idpp

!====================================================================
! Private
!====================================================================

subroutine optmz_steepest_descent(mode,flag_out,nsteps,stepsize,tol,&
    &writegp,fixed,savelastgeom,verbose)

  integer,                   intent(IN)  :: mode
  logical,                   intent(OUT) :: flag_out ! true if convergent, false otherwise
  integer,   optional,       intent(IN)  :: nsteps
  real(DBL), optional,       intent(IN)  :: stepsize
  real(DBL), optional,       intent(IN)  :: tol
  integer,   optional,       intent(IN)  :: writegp
  logical,   optional,       intent(IN)  :: fixed
  logical,   optional,       intent(IN)  :: savelastgeom
  logical,   optional,       intent(IN)  :: verbose

  integer                                :: p_nsteps
  integer                                :: p_writegp
  real(DBL)                              :: p_stepsize
  real(DBL)                              :: p_tol
  logical                                :: p_fixed
  logical                                :: p_savelastgeom
  logical                                :: p_verbose

  integer                                :: i
  character(8)                           :: istr
  logical                                :: flag_converged
  logical,   allocatable, dimension(:)   :: total_conv
  real(DBL), allocatable, dimension(:,:) :: new_geom
  integer                                :: err_n
  character(120)                         :: err_msg

  ! checking arguments ------------------------------------
  if (present(nsteps)) then
    if (nsteps>=0) then
      p_nsteps = nsteps
    else
      p_nsteps = -1
    end if
  else if (flag_optmz_nsteps) then
    p_nsteps = optmz_nsteps
  else
    p_nsteps = 500
  end if

  if (present(stepsize).and.(stepsize>0)) then
    p_stepsize = stepsize
  else
    p_stepsize = 1.0E-2_DBL
  end if

  if (present(tol).and.((tol>0.0_DBL).and.(tol<=1.0_DBL))) then
    p_tol = tol
  else if (flag_optmz_tol) then
    p_tol = optmz_tol
  else
    p_tol = 3.5E-4_DBL
  end if

  if (present(writegp).and.(writegp>0)) then
    p_writegp = writegp
  else
    p_writegp = -1
  end if

  if (present(fixed)) then
    p_fixed = fixed
  else
    p_fixed = .true.
  end if

  if (present(savelastgeom)) then
    p_savelastgeom = savelastgeom
  else
    p_savelastgeom = .true.
  end if

  if (present(verbose)) then
    p_verbose = verbose
  else
    p_verbose = .false.
  end if

  ! preliminary checks ------------------------------------
  if (flag_init_images.eqv..false.) then
    call error("optmz_steepest_descent: images not initialized")
  end if

  select case(mode)
  case (PES_MODE)
    write(FILEOUT,*) "**  Optimization Steepest Descent -- PES Mode"
  case (IDPP_MODE)
    write(FILEOUT,*) "**  Optimization Steepest Descent -- IDPP Mode"
  case default
    write(istr,'(I8)') mode
    istr = adjustl(istr)
    call error("optmz_steepest_descent: mode """//trim(istr)//""" not valid")
  end select

  ! allocation section ------------------------------------
  allocate(new_geom(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_steepest_descent: "//trim(err_msg))
  end if

  allocate(total_conv(image_n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_steepest_descent: "//trim(err_msg))
  end if

  ! working section ---------------------------------------
  flag_converged = .false.
  new_geom       = image_geom(1:image_n,:)

  i = 1
  do
    ! check loop index ------------------------------------
    if ((p_nsteps/=-1).and.(i>p_nsteps)) then
      exit
    end if

    write(FILEOUT,'(1X,A,I8,A)') "**  Optimization Steepest Descent -- Iteration ",i,":"

    ! do the job ------------------------------------------
    call compute_total_forces(mode,p_fixed)
    call total_forces_modifiers(mode,p_nsteps,i,p_fixed)
    new_geom = new_geom+(p_stepsize*total_forces)
    call update_images(new_geom)

    ! write the results -----------------------------------
    call write_opt_results(mode,total_conv,p_tol)

    if ((p_writegp/=-1).and.(mod(i,p_writegp)==0)) then
      call write_gnuplot_pes_energy(i)
    end if

    if (p_savelastgeom) then
      call last_geom_bkp(.true.)
    end if

    if (p_verbose) then
      call write_all_images(FILEOUT,.true.)
      call write_parall_elastic_forces()
      call write_perpen_pes_forces()
      call write_total_forces()
    end if

    ! check exit condition --------------------------------
    if (alltrue(total_conv)) then
      flag_converged = .true.
    end if
    
    if (flag_converged) then
      exit
    end if

    ! update variables and loop index ---------------------
    i          = i+1
  end do

  ! deallocation section ----------------------------------
  deallocate(new_geom,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_steepest_descent: "//trim(err_msg))
  end if

  deallocate(total_conv,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_steepest_descent: "//trim(err_msg))
  end if

  ! write optimization result, set output flag ------------
  if (flag_converged) then
    write(FILEOUT,*) "**  Optimization Steepest Descent -> Convergence Achieved"
    flag_out = .true.
  else
    write(FILEOUT,*) "**  Optimization Steepest Descent -> Convergence NOT Achieved"
    flag_out = .false.
  end if

end subroutine optmz_steepest_descent

!====================================================================

subroutine optmz_bfgs(mode,flag_out,nsteps,tol,fixed,savelastgeom)

  ! This is the "global" version of bfgs algorithm applied to NEB.
  ! It means that the images' geometries are passed all at once (globally)
  ! using a single large array.

  integer,                     intent(IN)  :: mode
  logical,                     intent(OUT) :: flag_out  ! true if convergent, false otherwise
  integer, optional,           intent(IN)  :: nsteps
  real(DBL), optional,         intent(IN)  :: tol
  logical, optional,           intent(IN)  :: fixed
  logical, optional,           intent(IN)  :: savelastgeom

  integer                                  :: p_nsteps
  real(DBL)                                :: p_tol
  logical                                  :: p_fixed
  logical                                  :: p_savelastgeom

  character(*), parameter                  :: my_name = "optmz_bfgs"
  character(8)                             :: istr
  character(120)                           :: cmdstr
  integer                                  :: sz_imggeom ! image_n * geom_len
  real(DBL), dimension(:,:), allocatable   :: x1  ! sz_imggeom x 1
  real(DBL), dimension(:,:), allocatable   :: h0  ! sz_imggeom x sz_imggeom
  real(DBL), dimension(:,:), allocatable   :: h1  ! sz_imggeom x sz_imggeom
  real(DBL), dimension(:,:), allocatable   :: old_image_geom   ! image_n x geom_len
  real(DBL), dimension(:,:), allocatable   :: new_image_geom   ! image_n x geom_len
  real(DBL), dimension(:,:), allocatable   :: old_total_forces ! image_n x geom_len
  real(DBL), dimension(:,:), allocatable   :: new_total_forces ! image_n x geom_len
  logical, dimension(:), allocatable       :: total_conv
  logical                                  :: flag_converged
  integer                                  :: i
  integer                                  :: j
  integer                                  :: err_n
  character(120)                           :: err_msg

  ! checking arguments ------------------------------------
  if (present(nsteps)) then
    if (nsteps>=0) then
      p_nsteps = nsteps
    else
      p_nsteps = -1
    end if
  else if (flag_optmz_nsteps) then
    p_nsteps = optmz_nsteps
  else
    p_nsteps = 50
  end if

  if (present(tol).and.((tol>0.0_DBL).and.(tol<=1.0_DBL))) then
    p_tol = tol
  else if (flag_optmz_tol) then
    p_tol = optmz_tol
  else
    p_tol = 3.5E-4_DBL
  end if

  if (present(fixed)) then
    p_fixed = fixed
  else
    p_fixed = .true.
  end if

  if (present(savelastgeom)) then
    p_savelastgeom = savelastgeom
  else
    p_savelastgeom = .true.
  end if

  ! preliminary checks ------------------------------------
  if (flag_init_images.eqv..false.) then
    call error(my_name//": images not initialized")
  end if

  ! init sz_imggeom ---------------------------------------
  sz_imggeom = image_n * geom_len

  ! allocation section ------------------------------------
  allocate(x1(sz_imggeom,1),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(h0(sz_imggeom,sz_imggeom),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(h1(sz_imggeom,sz_imggeom),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(old_image_geom(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(new_image_geom(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(old_total_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(new_total_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(total_conv(image_n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! write optimization parameters -------------------------
  select case(mode)
  case (PES_MODE)
    write(FILEOUT,*) "**  Optimization BFGS -- PES Mode"
  case (IDPP_MODE)
    write(FILEOUT,*) "**  Optimization BFGS -- IDPP Mode"
  case default
    write(istr,'(I8)') mode
    istr = adjustl(istr)
    call error(my_name//": mode """//trim(istr)//""" not valid")
  end select

  if (p_nsteps==-1) then
    write(FILEOUT,'(5X,"Max optimization iterations: infinite")')
  else
    write(FILEOUT,'(5X,"Max optimization iterations: ",I10)') p_nsteps
  end if

  write(FILEOUT,'(5X,"Convergence threshold      : ",ES8.1)') p_tol

  ! working section ---------------------------------------
  flag_converged = .false.

  i = 1
  do
    ! check loop index ------------------------------------
    if ((p_nsteps/=-1).and.(i>p_nsteps)) then
      exit
    end if

    ! start -----------------------------------------------
    write(FILEOUT,'(1X,A,I8,A)') "**  Optimization BFGS -- Iteration ",i,":"

    ! main body -------------------------------------------
    if (i == 1) then
      ! init h0 to identity matrix
      h0 = 0.0_DBL
      do j = 1, size(h0,1)
        h0(j,j) = 1.0_DBL
      end do

      ! preliminary computation
      call compute_total_forces(mode,p_fixed)
      call total_forces_modifiers(mode,p_nsteps,i,p_fixed)
      old_image_geom   = image_geom(1:image_n,:)
      old_total_forces = total_forces

      cmdstr = "START"
      call bfgs_internal(                             &
        cmdstr,                                       &
        reshape(old_image_geom,(/sz_imggeom, 1/)),    &
        x1,                                           & ! used as output
        -reshape(old_total_forces,(/sz_imggeom, 1/)), &
        -reshape(new_total_forces,(/sz_imggeom, 1/)), & ! not used
        h0,                                           &
        h1,                                           &
        fixed_alpha = .true.,                         &
        reset_alpha = .true.                          &
      )
    else
      cmdstr = "START"
      call bfgs_internal(                             &
        cmdstr,                                       &
        reshape(old_image_geom,(/sz_imggeom, 1/)),    &
        x1,                                           & ! used as output
        -reshape(old_total_forces,(/sz_imggeom, 1/)), &
        -reshape(new_total_forces,(/sz_imggeom, 1/)), & ! not used
        h0,                                           &
        h1,                                           &
        fixed_alpha = .true.                          &
      )
    end if

    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": expected ""EVALUATE_DF1"", get """//trim(cmdstr)//"""")
    end if

    new_image_geom = reshape(x1,(/image_n, geom_len/))

    call update_images(new_image_geom)
    call compute_total_forces(mode,p_fixed)
    call total_forces_modifiers(mode,p_nsteps,i,p_fixed)
    new_total_forces = total_forces

    cmdstr = "EVALUATED"
    call bfgs_internal(                             &
      cmdstr,                                       &
      reshape(old_image_geom,(/sz_imggeom, 1/)),    &
      x1,                                           & ! this time used as input
      -reshape(old_total_forces,(/sz_imggeom, 1/)), &
      -reshape(new_total_forces,(/sz_imggeom, 1/)), & ! used as input
      h0,                                           &
      h1,                                           &
      fixed_alpha = .true.                          &
    )

    select case (cmdstr)
    case ("DONE")
      old_image_geom   = new_image_geom
      old_total_forces = new_total_forces
      h0 = h1

      ! write the results ---------------------------------
      call write_opt_results(mode,total_conv,p_tol)

      ! write geometry file -------------------------------
      if (p_savelastgeom) then
        call last_geom_bkp(.true.)
      end if

      ! check exit condition ------------------------------
      if (alltrue(total_conv)) then
        flag_converged = .true.
      end if

      if (flag_converged) then
        exit
      end if
    case ("SKIPPED")
      write(FILEOUT,'("         skipped")')
    case default
      call error(my_name//": expected ""DONE"" or ""SKIPPED"", get """//trim(cmdstr)//"""")
    end select

    ! update loop index -----------------------------------
    i = i+1
  end do

  ! deallocation section ----------------------------------
  deallocate(x1,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(h0,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(h1,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(old_image_geom,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(new_image_geom,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(old_total_forces,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(new_total_forces,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(total_conv,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! write optimization result, set output flag ------------
  if (flag_converged) then
    write(FILEOUT,*) "**  Optimization BFGS -> Convergence Achieved"
    flag_out = .true.
  else
    write(FILEOUT,*) "**  Optimization BFGS -> Convergence NOT Achieved"
    flag_out = .false.
  end if

end subroutine optmz_bfgs

!====================================================================

subroutine optmz_lbfgs(mode,flag_out,nsteps,memsize,prec,tol,verblvl,&
    &fixed,savelastgeom)

  integer,                   intent(IN)  :: mode
  logical,                   intent(OUT) :: flag_out         ! true if convergent, false otherwise
  integer,   optional,       intent(IN)  :: nsteps
  integer,   optional,       intent(IN)  :: memsize
  real(DBL), optional,       intent(IN)  :: prec
  real(DBL), optional,       intent(IN)  :: tol
  integer,   optional,       intent(IN)  :: verblvl
  logical,   optional,       intent(IN)  :: fixed
  logical,   optional,       intent(IN)  :: savelastgeom

  integer                                :: p_nsteps
  logical                                :: p_fixed
  logical                                :: p_savelastgeom

  integer                                :: n
  integer                                :: m
  integer                                :: iprint
  real(DBL)                              :: factr
  real(DBL)                              :: pgtol
  real(DBL), allocatable, dimension(:)   :: f                ! energy
  real(DBL), allocatable, dimension(:,:) :: x                ! coordinates
  real(DBL), allocatable, dimension(:,:) :: g                ! gradient
  character(60)                          :: task
  logical, allocatable, dimension(:)     :: image_conv
  logical                                :: flag_converged
  logical                                :: flag_first_bfgs
  integer                                :: i
  integer                                :: j
  character(8)                           :: j_str
  character(8)                           :: istr
  character(120)                         :: err_msg
  integer                                :: err_n

  ! checking arguments ------------------------------------
  if (present(nsteps)) then
    if (nsteps>=0) then
      p_nsteps = nsteps
    else
      p_nsteps = -1
    end if
  else if (flag_optmz_nsteps) then
    p_nsteps = optmz_nsteps
  else
    p_nsteps = 50
  end if

  if (present(memsize).and.((memsize>=3).and.(memsize<=20))) then
    m = memsize
  else
    m = 15
  end if

  if (present(prec).and.(prec>=1.0_DBL)) then
    factr = prec
  else
    factr = 1.0E+7_DBL
  end if

  if (present(tol).and.((tol>0.0_DBL).and.(tol<=1.0_DBL))) then
    pgtol = tol
  else if (flag_optmz_tol) then
    pgtol = optmz_tol
  else
    pgtol = 3.5E-4_DBL
  end if

  if (present(verblvl).and.(verblvl>=0)) then
    iprint = verblvl
  else
    iprint = -1
  end if

  if (present(fixed)) then
    p_fixed = fixed
  else
    p_fixed = .true.
  end if

  if (present(savelastgeom)) then
    p_savelastgeom = savelastgeom
  else
    p_savelastgeom = .true.
  end if

  ! preliminary checks ------------------------------------
  if (flag_init_images.eqv..false.) then
    call error("optmz_lbfgs: images not initialized")
  end if

  select case(mode)
  case (PES_MODE)
    write(FILEOUT,*) "**  Optimization l-BFGS -- PES Mode"
  case (IDPP_MODE)
    write(FILEOUT,*) "**  Optimization l-BFGS -- IDPP Mode"
  case default
    write(istr,'(I8)') mode
    istr = adjustl(istr)
    call error("optmz_lbfgs: mode """//trim(istr)//""" not valid")
  end select

  ! allocation section ------------------------------------
  allocate(f(image_n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_lbfgs: "//trim(err_msg))
  end if

  allocate(x(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_lbfgs: "//trim(err_msg))
  end if

  allocate(g(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_lbfgs: "//trim(err_msg))
  end if

  allocate(image_conv(image_n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_lbfgs: "//trim(err_msg))
  end if

  ! init section ------------------------------------------
  write(FILEOUT,*) "**  Optimization l-BFGS -- Initialization"
  call compute_total_forces(mode,p_fixed) ! needed because tangents estimation
    ! done by arbitrary_geom_total_forces requires energy values.

  n              = geom_len
  x              = image_geom(1:image_n,:)
  flag_converged = .false.

  ! optimization loop -------------------------------------
  i = 1
  external_lbfgs: do
    ! check loop index ------------------------------------
    if ((p_nsteps/=-1).and.(i>p_nsteps)) then
      exit
    end if

    write(FILEOUT,'(1X,A,I8,A)') "**  Optimization l-BFGS -- Iteration ",i,":"
    
    ! optimize every image with l-BFGS --------------------
    image_conv      = .false.
    flag_first_bfgs = .true.

    j = 1
    internal_lbfgs: do while (j<=image_n)
      if (flag_first_bfgs) then
        task            = "START"
        flag_first_bfgs = .false.
      end if

      ! call lbfgs subroutine -----------------------------
      call lbfgs(n,m,x(j,:),f(j),g(j,:),factr,pgtol,task,iprint)
      write(j_str,'(I8)') j
      write(FILEOUT,*) "lbfgs(image=",trim(adjustl(j_str)),") returned: ",trim(task)

      ! check lbfgs result --------------------------------
      if ((task(1:5)=="NEW_X").or.(task(1:4)=="CONV")) then
        if (task(1:4)=="CONV") then
          image_conv(j) = .true.
        end if
        j               = j+1
        flag_first_bfgs = .true.
      else if (task(1:2)=="FG") then
        call arbitrary_geom_total_forces(j,x(j,:),f(j),g(j,:))
        g(j,:) = -g(j,:) ! g is the gradient,
          ! but arbitrary_geom_total_forces returns the forces
      else if ((task(1:4)=="ABNO").or.(task(1:5)=="ERROR")) then
        call error("optmz_lbfgs: "//trim(task))
      else
        call error("optmz_lbfgs: unknown termination string: "//&
          &trim(task))
      end if
    end do internal_lbfgs

    ! update geometries and energies ----------------------
    call update_images(x)
    call update_pes_energy(f)

    ! write the results -----------------------------------
    write(FILEOUT,'(5X,"Image   Conv")')
    do j=1, image_n
      write(FILEOUT,'(7X,I3,4X,L1)') j, image_conv(j)
    end do

    if (p_savelastgeom) then
      call last_geom_bkp(.true.)
    end if

    ! check exit condition --------------------------------
    if (alltrue(image_conv)) then
      flag_converged=.true.
    end if

    if (flag_converged) then
      exit
    end if
    
    ! update loop index -----------------------------------
    i = i+1
  end do external_lbfgs

  ! deallocation section ----------------------------------
  deallocate(f,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_lbfgs: "//trim(err_msg))
  end if

  deallocate(x,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_lbfgs: "//trim(err_msg))
  end if

  deallocate(g,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_lbfgs: "//trim(err_msg))
  end if

  deallocate(image_conv,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_lbfgs: "//trim(err_msg))
  end if

  ! write optimization result, set output flag ------------
  if (flag_converged) then
    write(FILEOUT,*) "**  Optimization l-BFGS -> Convergence Achieved"
    flag_out = .true.
  else
    write(FILEOUT,*) "**  Optimization l-BFGS -> Convergence NOT Achieved"
    flag_out = .false.
  end if

end subroutine optmz_lbfgs

!====================================================================

subroutine optmz_fire(mode,flag_out,nsteps,maxstepsize,tol,&
    &fixed,savelastgeom)

  integer,                   intent(IN)  :: mode
  logical,                   intent(OUT) :: flag_out       ! true if convergent, false otherwise
  integer,   optional,       intent(IN)  :: nsteps
  real(DBL), optional,       intent(IN)  :: maxstepsize
  real(DBL), optional,       intent(IN)  :: tol
  logical,   optional,       intent(IN)  :: fixed
  logical,   optional,       intent(IN)  :: savelastgeom

  integer                                :: p_nsteps
  real(DBL)                              :: p_tol
  logical                                :: p_fixed
  logical                                :: p_savelastgeom

  real(DBL), parameter                   :: dtdec          = 0.5_DBL
  real(DBL), parameter                   :: dtinc          = 1.1_DBL
  real(DBL), parameter                   :: alpha_start    = 0.1_DBL
  real(DBL), parameter                   :: alpha_dec      = 0.99_DBL
  integer,   parameter                   :: stepi_min      = 5
  logical,   allocatable, dimension(:)   :: total_conv
  real(DBL), allocatable, dimension(:,:) :: x
  real(DBL), allocatable, dimension(:,:) :: dx
  real(DBL), allocatable, dimension(:,:) :: velocity
  real(DBL)                              :: dt
  real(DBL)                              :: dtmax
  real(DBL)                              :: power
  real(DBL)                              :: alpha
  logical                                :: flag_converged
  integer                                :: i
  integer                                :: stepi
  character(8)                           :: istr
  integer                                :: err_n
  character(120)                         :: err_msg

  ! checking arguments ------------------------------------
  if (present(nsteps)) then
    if (nsteps>=0) then
      p_nsteps = nsteps
    else
      p_nsteps = -1
    end if
  else if (flag_optmz_nsteps) then
    p_nsteps = optmz_nsteps
  else
    p_nsteps = 500
  end if

  if (present(maxstepsize).and.(maxstepsize>0)) then
    dtmax = min(maxstepsize,0.3_DBL)
  else
    dtmax = 0.3_DBL
  end if

  if (present(tol).and.((tol>0.0_DBL).and.(tol<=1.0_DBL))) then
    p_tol = tol
  else if (flag_optmz_tol) then
    p_tol = optmz_tol
  else
    p_tol = 3.5E-4_DBL
  end if

  if (present(fixed)) then
    p_fixed = fixed
  else
    p_fixed = .true.
  end if

  if (present(savelastgeom)) then
    p_savelastgeom = savelastgeom
  else
    p_savelastgeom = .true.
  end if

  ! preliminary checks ------------------------------------
  if (flag_init_images.eqv..false.) then
    call error("optmz_fire: images not initialized")
  end if

  ! allocation section ------------------------------------
  allocate(total_conv(image_n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_fire: "//trim(err_msg))
  end if

  allocate(x(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_fire: "//trim(err_msg))
  end if

  allocate(dx(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_fire: "//trim(err_msg))
  end if

  allocate(velocity(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_fire: "//trim(err_msg))
  end if

  ! write optimization parameters -------------------------
  select case(mode)
  case (PES_MODE)
    write(FILEOUT,*) "**  Optimization FIRE -- PES Mode"
  case (IDPP_MODE)
    write(FILEOUT,*) "**  Optimization FIRE -- IDPP Mode"
  case default
    write(istr,'(I8)') mode
    istr = adjustl(istr)
    call error("optmz_fire: mode """//trim(istr)//""" not valid")
  end select

  if (p_nsteps==-1) then
    write(FILEOUT,'(5X,"Max optimization iterations: infinite")')
  else
    write(FILEOUT,'(5X,"Max optimization iterations: ",I10)') p_nsteps
  end if

  write(FILEOUT,'(5X,"Convergence threshold      : ",ES8.1)') p_tol

  ! init section ------------------------------------------
  flag_converged = .false.
  dt             = dtmax/10.0
  alpha          = alpha_start
  stepi          = 0
  velocity       = 0.0_DBL
  x              = image_geom(1:image_n,:)

  ! working section ---------------------------------------
  i = 1
  do
    ! check loop index ------------------------------------
    if ((p_nsteps/=-1).and.(i>p_nsteps)) then
      exit
    end if

    ! get total forces ------------------------------------
    write(FILEOUT,'(1X,A,I8,A)') "**  Optimization FIRE -- Iteration ",i,":"
    call compute_total_forces(mode,p_fixed)
    call total_forces_modifiers(mode,p_nsteps,i,p_fixed)

    ! FIRE algorithm --------------------------------------
    if (sqrt(sum(velocity*velocity))/=0.0_DBL) then
      power = sum(total_forces*velocity)
      if (power>0.0_DBL) then
        velocity = (1.0_DBL-alpha)*velocity+&                                                          
          &alpha*total_forces/sqrt(sum(total_forces*total_forces))*&
          &sqrt(sum(velocity*velocity))

        if (stepi>=stepi_min) then
          dt    = min(dt*dtinc,dtmax)
          alpha = alpha*alpha_dec
        end if

        stepi = stepi+1
      else
        velocity = 0.0_DBL
        alpha    = alpha_start
        dt       = dt*dtdec
        stepi    = 0
      end if
    end if

    velocity = velocity+(dt*total_forces)
    dx       = dt*velocity
    x        = x+dx
    call update_images(x)
 
    ! write the results -----------------------------------
    call write_opt_results(mode,total_conv,p_tol)

    ! write geometry file ---------------------------------
    if (p_savelastgeom) then
      call last_geom_bkp(.true.)
    end if

    ! check exit condition --------------------------------
    if (alltrue(total_conv)) then
      flag_converged = .true.
    end if

    if (flag_converged) then
      exit
    end if
    
    ! update loop index -----------------------------------
    i = i+1
  end do

  ! deallocation section ----------------------------------
  deallocate(total_conv,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_fire: "//trim(err_msg))
  end if

  deallocate(x,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_fire: "//trim(err_msg))
  end if

  deallocate(dx,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_fire: "//trim(err_msg))
  end if

  deallocate(velocity,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("optmz_fire: "//trim(err_msg))
  end if

  ! write optimization result, set output flag ------------
  if (flag_converged) then
    write(FILEOUT,*) "**  Optimization FIRE -> Convergence Achieved"
    flag_out = .true.
  else
    write(FILEOUT,*) "**  Optimization FIRE -> Convergence NOT Achieved"
    flag_out = .false.
  end if

end subroutine optmz_fire

!====================================================================

subroutine total_forces_modifiers(mode,p_nsteps,i,p_fixed)

  !--------------------------------------------------------
  ! Called by optimization subroutines, it checks
  ! if proper conditions are satisfied
  ! to apply some forces modifier procedures,
  ! like climbing or descending image.
  !--------------------------------------------------------

  integer,  intent(IN) :: mode                       ! caller operating mode
  integer,  intent(IN) :: p_nsteps                   ! caller iteration total steps
  integer,  intent(IN) :: i                          ! caller current iteration
  logical,  intent(IN) :: p_fixed                    ! manipulates degrees of freedom

  real(DBL), parameter :: cap          = 0.1_DBL     ! climbing activation percentage
  integer,   parameter :: cam          = 10          ! climbing activation maximum steps
  logical              :: flag_climbing_iteration    ! can i do climbing?
  logical              :: flag_descending_iteration  ! can i do descending?
  logical              :: write_output = .true.      ! can called subroutines write output?

  ! works only in PES_MODE --------------------------------
  if (mode/=PES_MODE) then
    return
  end if

  ! checks if got climbing prerequisites ------------------
  if (flag_climbing_image) then
    if (flag_climbing_quick_start) then
      flag_climbing_iteration = .true.
    else
      if (p_nsteps==-1) then
        if (i>cam) then
          flag_climbing_iteration = .true.
        else
          flag_climbing_iteration = .false.
        end if
      else
        if (i>min(floor(p_nsteps*cap),cam)) then
          flag_climbing_iteration = .true.
        else
          flag_climbing_iteration = .false.
        end if
      end if
    end if

    if (flag_climbing_iteration) then
      write(FILEOUT,'(5X,"Climbing on current iteration: Active")')
      call exec_climbing(p_fixed,write_output)
    else
      write(FILEOUT,'(5X,"Climbing on current iteration: Not Active")')
    end if
  end if

  ! checks if got descending prerequisites ----------------
  if (flag_descending_image) then
    if (flag_descending_quick_start) then
      flag_descending_iteration = .true.
    else
      if (p_nsteps==-1) then
        if (i>cam) then
          flag_descending_iteration = .true.
        else
          flag_descending_iteration = .false.
        end if
      else
        if (i>min(floor(p_nsteps*cap),cam)) then
          flag_descending_iteration = .true.
        else
          flag_descending_iteration = .false.
        end if
      end if
    end if

    if (flag_descending_iteration) then
      write(FILEOUT,'(5X,"Descending on current iteration: Active")')
      call exec_descending(p_fixed,write_output)
    else
      write(FILEOUT,'(5X,"Descending on current iteration: Not Active")')
    end if
  end if

end subroutine total_forces_modifiers

!====================================================================

subroutine write_opt_results(mode,total_conv,p_tol)

  integer,               intent(IN)  :: mode
  logical, dimension(:), intent(OUT) :: total_conv
  real(DBL),             intent(IN)  :: p_tol

  character(*), parameter            :: my_name = "write_opt_results"
  integer                            :: j
  real(DBL)                          :: rtmp

  ! preliminary checks ------------------------------------
  if (size(total_conv,1)/=image_n) then
    call error(my_name//": wrong size of argument total_conv")
  end if

  ! write the results -------------------------------------
  total_conv = .false.
  write(FILEOUT,'(5X,"Image",12X,"Energy",9X,"Tot Force",4X,"Conv (",ES10.3,")")') p_tol
  do j=1, image_n
    rtmp = norm(total_forces(j,:))
    if (rtmp<p_tol) then
      total_conv(j) = .true.
    end if

    select case (mode)
    case (PES_MODE)
      write(FILEOUT,'(7X,I3,3X,F15.6,8X,ES10.3,7X,L1)') j, pes_energy(j), rtmp, total_conv(j)
    case (IDPP_MODE)
      write(FILEOUT,'(7X,I3,3X,F15.6,8X,ES10.3,7X,L1)') j, idpp_energy(j), rtmp, total_conv(j)
    end select
  end do

end subroutine write_opt_results

!====================================================================

end module optimization

