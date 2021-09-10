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

module bfgs

  use utility

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: bfgs_internal

contains

!====================================================================
! Public
!====================================================================

subroutine bfgs_internal(cmdstr,mode,x0,x1,df0,df1,h0,h1,reset_alpha)

  ! BFSG method, based on:
  !
  ! Nocedal - Numerical Optimization - 2nd edition
  ! Algorithm 6.1 (only the statements inside the loop)
  !
  ! Call to accelerated_backtracking_line_search() based on:
  !
  ! Herbol2017computational (DOI: 10.1021/acs.jctc.7b00360)
  ! Algorithm 9


  character(120),            intent(INOUT) :: cmdstr
  character(*),              intent(IN)    :: mode
  real(DBL), dimension(:,:), intent(IN)    :: x0
  real(DBL), dimension(:,:), intent(INOUT) :: x1
  real(DBL), dimension(:,:), intent(IN)    :: df0
  real(DBL), dimension(:,:), intent(IN)    :: df1
  real(DBL), dimension(:,:), intent(IN)    :: h0
  real(DBL), dimension(:,:), intent(OUT)   :: h1
  logical, optional,         intent(IN)    :: reset_alpha

  character(*), parameter                  :: my_name          = "bfgs_internal"
  real(DBL), parameter                     :: max_displacement = 0.2_DBL
  real(DBL), parameter                     :: alpha0           = 1.0_DBL
  real(DBL), save                          :: alpha            = alpha0
  real(DBL)                                :: h0_scaling
  real(DBL)                                :: p0_max
  integer                                  :: sz1_x0
  integer                                  :: sz2
  real(DBL), dimension(:,:), allocatable   :: p0
  real(DBL), dimension(:,:), allocatable   :: s0
  real(DBL), dimension(:,:), allocatable   :: y0
  real(DBL)                                :: r0
  real(DBL), dimension(:,:), allocatable   :: idnt
  real(DBL), dimension(:,:), allocatable   :: m1
  real(DBL), dimension(:,:), allocatable   :: m2
  logical                                  :: flag_fixed_alpha
  logical                                  :: flag_accelerated_bt
  logical                                  :: flag_max_displacement
  logical                                  :: flag_skip
  logical                                  :: flag_reset_n
  integer                                  :: i
  integer                                  :: err_n
  character(120)                           :: err_msg

  ! check optional argument -------------------------------
  flag_reset_n = .false.
  if (present(reset_alpha)) then
    if (reset_alpha.eqv..true.) then
      alpha = alpha0
      flag_reset_n = .true.
    end if
  end if

  ! mode check --------------------------------------------
  flag_fixed_alpha      = .false.
  flag_accelerated_bt   = .false.
  flag_max_displacement = .false.

  select case (mode)
  case ("fixed_alpha")
    flag_fixed_alpha      = .true.
  case ("accelerated_bt")
    flag_accelerated_bt   = .true.
  case ("max_displacement")
    flag_max_displacement = .true.
  case default
    call error(my_name//": unknown mode """//trim(mode)//"""")
  end select

  ! arguments' dimensions check ---------------------------
  sz1_x0 = size(x0,1)
  sz2    = size(x0,2)
  if (sz2 /= 1) then
    call error(my_name//": argument x0 must be a column vector (n x 1)")
  end if

  sz2 = size(x1,2)
  if (size(x1,1) /= sz1_x0) then
    call error(my_name//": size mismatch between arguments x0 and x1")
  end if
  if (sz2 /= 1) then
    call error(my_name//": argument x1 must be a column vector (n x 1)")
  end if

  sz2 = size(df0,2)
  if (size(df0,1) /= sz1_x0) then
    call error(my_name//": size mismatch between arguments x0 and df0")
  end if
  if (sz2 /= 1) then
    call error(my_name//": argument df0 must be a column vector (n x 1)")
  end if

  sz2 = size(df1,2)
  if (size(df1,1) /= sz1_x0) then
    call error(my_name//": size mismatch between arguments x0 and df1")
  end if
  if (sz2 /= 1) then
    call error(my_name//": argument df1 must be a column vector (n x 1)")
  end if

  sz2 = size(h0,2)
  if (size(h0,1) /= sz1_x0) then
    call error(my_name//": size mismatch between arguments x0 and h0")
  end if
  if (sz2 /= size(h0,1)) then
    call error(my_name//": argument h0 must be a square matrix (n x n)")
  end if

  sz2 = size(h1,2)
  if (size(h1,1) /= sz1_x0) then
    call error(my_name//": size mismatch between arguments x0 and h1")
  end if
  if (sz2 /= size(h1,1)) then
    call error(my_name//": argument h1 must be a square matrix (n x n)")
  end if

  ! allocation section ------------------------------------
  allocate(p0(sz1_x0,1),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(s0(sz1_x0,1),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(y0(sz1_x0,1),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(idnt(sz1_x0,sz1_x0),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(m1(sz1_x0,sz1_x0),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(m2(sz1_x0,sz1_x0),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! init identity matrix ----------------------------------
  idnt = 0.0_DBL
  do i = 1, sz1_x0
    idnt(i,i) = 1.0_DBL
  end do

  ! cmdstr parsing ----------------------------------------
  ! cmdstr = {START, EVALUATE_DF1, EVALUATED, SKIPPED, DONE, ERROR}
  select case (cmdstr)
  case ("START")
    ! compute search direction p0
    p0 = -matmul(h0,df0)

    if (flag_max_displacement) then
      alpha  = alpha0
      p0_max = maxval(p0)
      if (alpha*p0_max > max_displacement) then
        alpha = max_displacement / p0_max
      end if
    end if

    ! set x1
    x1 = x0 + alpha*p0

    ! need df1 to proceed
    cmdstr = "EVALUATE_DF1"
  case ("EVALUATED")
    if (flag_accelerated_bt) then
      call accelerated_backtracking_line_search(       &
        rms(reshape(df0,(/size(df0,1)*size(df0,2)/))), &
        rms(reshape(df1,(/size(df1,1)*size(df1,2)/))), &
        alpha,                                         &
        flag_skip,                                     &
        flag_reset_n                                   &
      )
    else
      flag_skip = .false.
    end if

    if (flag_skip) then
      cmdstr = "SKIPPED"
    else
      ! set s0, y0 and r0
      s0 = x1 - x0
      y0 = df1 - df0
      r0 = 1.0_DBL / sum(matmul(transpose(y0),s0))

      ! scale h0
      if (isidentity(h0)) then
        h0_scaling = sum(matmul(transpose(y0),s0)/matmul(transpose(y0),y0))
      else
        h0_scaling = 1.0_DBL
      end if

      ! compute h1
      m1 = idnt - r0*matmul(s0,transpose(y0))
      m2 = idnt - r0*matmul(y0,transpose(s0))
      h1 = matmul(matmul(m1,h0*h0_scaling),m2) + r0*matmul(s0,transpose(s0))

      ! done
      cmdstr = "DONE"
    end if
  case default
    cmdstr = "ERROR"
    call error(my_name//": unknown cmdstr "//trim(cmdstr))
  end select

  ! deallocation section ----------------------------------
  deallocate(p0,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(s0,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(y0,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(idnt,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(m1,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(m2,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

end subroutine bfgs_internal

!====================================================================
! Private
!====================================================================

subroutine accelerated_backtracking_line_search(df0_rms,df1_rms,alpha,flag_skip,flag_reset_n)

  ! Based on:
  !
  ! Herbol2017computational (DOI: 10.1021/acs.jctc.7b00360)
  ! Algorithm 4

  real(DBL),  intent(IN)    :: df0_rms
  real(DBL),  intent(IN)    :: df1_rms
  real(DBL),  intent(INOUT) :: alpha
  logical,    intent(OUT)   :: flag_skip
  logical,    intent(IN)    :: flag_reset_n

  character(*), parameter   :: my_name   = "accelerated_backtracking_line_search"
  real(DBL),    parameter   :: threshold = 1.0E-3_DBL
  real(DBL),    parameter   :: alpha0    = 1.0_DBL
  real(DBL),    parameter   :: gam       = 0.5_DBL
  integer,      parameter   :: n0        = 100  ! accelerate alpha after n0 successful iterations
  real(DBL)                 :: chk
  integer, save             :: n_back    = n0

  if (flag_reset_n.eqv..true.) then
    n_back = n0
  end if

  chk = (df1_rms - df0_rms) / abs(df1_rms + df0_rms)
  flag_skip = .false.

  if (chk > threshold) then
    alpha = alpha * gam  ! slow alpha
    flag_skip = .true.
    n_back = n0
  else
    n_back = n_back - 1
    if (n_back < 1) then
      n_back = n0
      if (alpha < alpha0) then
        alpha = alpha0  ! reset alpha
        flag_skip = .true.
      else
        alpha = alpha / gam  ! accelerate alpha
      end if
    end if
  end if

end subroutine accelerated_backtracking_line_search

!====================================================================

end module bfgs

