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

  ! drivers -----------------------------------------------
  public :: driver_bfgs_rosenbrock

contains

!====================================================================
! Public
!====================================================================

subroutine bfgs_internal(cmdstr,x0,x1,df0,df1,h0,h1,wolfe_cond1,wolfe_cond2)

  ! BFSG method, based on:
  ! Nocedal - Numerical Optimization - 2nd edition
  ! Algorithm 6.1 (only the statements inside the loop)

  character(120),            intent(INOUT) :: cmdstr
  real(DBL), dimension(:,:), intent(IN)    :: x0
  real(DBL), dimension(:,:), intent(INOUT) :: x1
  real(DBL), dimension(:,:), intent(IN)    :: df0
  real(DBL), dimension(:,:), intent(IN)    :: df1
  real(DBL), dimension(:,:), intent(IN)    :: h0
  real(DBL), dimension(:,:), intent(OUT)   :: h1
  logical,                   intent(IN)    :: wolfe_cond1
  logical,                   intent(IN)    :: wolfe_cond2

  character(*), parameter                  :: my_name = "bfgs_internal"
  integer                                  :: sz1_x0
  integer                                  :: sz2
  real(DBL), save                          :: alpha
  real(DBL), dimension(:,:), allocatable   :: p0
  real(DBL), dimension(:,:), allocatable   :: s0
  real(DBL), dimension(:,:), allocatable   :: y0
  real(DBL)                                :: r0
  real(DBL), dimension(:,:), allocatable   :: idnt
  real(DBL), dimension(:,:), allocatable   :: m1
  real(DBL), dimension(:,:), allocatable   :: m2
  integer                                  :: i
  integer                                  :: err_n
  character(120)                           :: err_msg

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
  ! cmdstr = {START, EVALUATE_DF1, EVALUATED, DONE, ERROR}
  select case (cmdstr)
  case ("START")
    ! compute search direction p0
    p0 = -matmul(h0,df0)

    ! get step length alpha satisfying strong wolfe conditions
    if ((wolfe_cond1.eqv..false.).and.(wolfe_cond2.eqv..false.)) then
      alpha = 1.0_DBL
    else
      call strong_wolfe_conditions(wolfe_cond1,wolfe_cond2,alpha)
    end if

    ! set x1
    x1 = x0 + alpha*p0

    ! need df1 to proceed
    cmdstr = "EVALUATE_DF1"
  case ("EVALUATED")
    ! set s0, y0 and r0
    s0 = x1 - x0
    y0 = df1 - df0
    r0 = sum((matmul(transpose(y0),s0))**(-1))

    ! compute h1
    m1 = idnt - r0*matmul(s0,transpose(y0))
    m2 = idnt - r0*matmul(y0,transpose(s0))
    h1 = matmul(matmul(m1,h0),m2) + r0*matmul(s0,transpose(s0))

    ! done
    cmdstr = "DONE"
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

subroutine strong_wolfe_conditions(cond1,cond2,alpha)

  ! Strong Wolfe conditions, based on:
  ! Nocedal - Numerical Optimization - 2nd edition
  ! Eq. (3.7a), (3.7b)

  logical,    intent(IN)  :: cond1
  logical,    intent(IN)  :: cond2
  real(DBL),  intent(OUT) :: alpha

  character(*), parameter :: my_name = "strong_wolfe_conditions"

  call error(my_name//": subroutine not yet implemented")

  ! dummy logic to suppress warnings
  if (cond1.eqv.cond2) alpha = 1.0_DBL

  !TODO implement logic

end subroutine strong_wolfe_conditions

!====================================================================
! Driver
!====================================================================

subroutine driver_bfgs_rosenbrock()

  character(*), parameter   :: my_name = "driver_bfgs_rosenbrock"
  character(*), parameter   :: header  = "(2X,""Iter"",3X,""x1"",14X,""x2"",14X,""f(x)"",12X, &
                                          &""d/dx1 f(x)"",6X,""d/dx2 f(x)"",6X,""Norm ("",ES8.1,"")"")"
  character(*), parameter   :: sep     = "(103(""-""))"
  character(*), parameter   :: format1 = "(2X,I4,3X,5(ES13.6,3X))"
  character(*), parameter   :: format2 = "(2X,I4,3X,6(ES13.6,3X))"
  real(DBL), parameter      :: tol     = 1.0E-5_DBL
  character(120)            :: cmdstr
  real(DBL), dimension(2,1) :: old_x
  real(DBL)                 :: old_f
  real(DBL), dimension(2,1) :: old_d
  real(DBL), dimension(2,2) :: old_h
  real(DBL), dimension(2,1) :: new_x
  real(DBL)                 :: new_f
  real(DBL), dimension(2,1) :: new_d
  real(DBL), dimension(2,2) :: new_h
  logical                   :: wolfe_cond1
  logical                   :: wolfe_cond2
  real(DBL), dimension(2)   :: buff
  real(DBL)                 :: nrm
  integer                   :: i

  write(*,sep)
  write(*,*) my_name//": BFGS without Wolfe conditions"
  write(*,header) tol

  wolfe_cond1 = .false.
  wolfe_cond2 = .false.
  old_x(1,1) = -1.2_DBL
  old_x(2,1) =  1.0_DBL

  i = 1
  do
    if (i == 1) then
      old_h      = 0.0_DBL
      old_h(1,1) = 1.0_DBL
      old_h(2,2) = 1.0_DBL

      call rosenbrock_f(reshape(old_x,(/2/)),old_f)
      call rosenbrock_d(reshape(old_x,(/2/)),buff)
      old_d = reshape(buff,(/2,1/))
      write(*,format1) i, old_x, old_f, old_d
    end if

    cmdstr = "START"
    call bfgs_internal( &
      cmdstr,           &
      old_x,            &
      new_x,            &
      old_d,            &
      new_d,            &
      old_h,            &
      new_h,            &
      wolfe_cond1,      &
      wolfe_cond2       &
    )

    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": expected ""EVALUATE_DF1"", get """//trim(cmdstr)//"""")
    end if

    call rosenbrock_f(reshape(new_x,(/2/)),new_f)
    call rosenbrock_d(reshape(new_x,(/2/)),buff)
    new_d = reshape(buff,(/2,1/))

    cmdstr = "EVALUATED"
    call bfgs_internal( &
      cmdstr,           &
      old_x,            &
      new_x,            &
      old_d,            &
      new_d,            &
      old_h,            &
      new_h,            &
      wolfe_cond1,      &
      wolfe_cond2       &
    )

    if (cmdstr /= "DONE") then
      call error(my_name//": expected ""DONE"", get """//trim(cmdstr)//"""")
    end if

    nrm = sqrt(sum(matmul(transpose(new_d),new_d)))
    write(*,format2) i, new_x, new_f, new_d, nrm

    if (nrm < tol) then
      exit
    end if

    old_x = new_x
    old_f = new_f
    old_d = new_d
    old_h = new_h

    i = i+1
  end do

  write(*,sep)

end subroutine driver_bfgs_rosenbrock

!====================================================================

subroutine rosenbrock_f(x,f)

  real(DBL), dimension(2), intent(IN)  :: x
  real(DBL),               intent(OUT) :: f

  f = 100.0_DBL*(x(2) - x(1)**2)**2 + (1.0_DBL - x(1))**2

end subroutine rosenbrock_f

!====================================================================

subroutine rosenbrock_d(x,d)

  real(DBL), dimension(2), intent(IN)  :: x
  real(DBL), dimension(2), intent(OUT) :: d

  d(1) = 400.0_DBL*x(1)**3 - 400.0_DBL*x(1)*x(2) - 2.0_DBL + 2.0_DBL*x(1)
  d(2) = 200.0_DBL*x(2) - 200.0_DBL*x(1)**2

end subroutine rosenbrock_d

!====================================================================

end module bfgs