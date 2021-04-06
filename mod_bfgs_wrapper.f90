!  Copyright (C) 2020-2021     Marco Bertini
!  
!  This file is part of neb.x.
!  
!  neb.x is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  neb.x is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with neb.x.  If not, see <https://www.gnu.org/licenses/>.

module bfgs_wrapper

  !==================================================================
  !   BFGS Wrapper Module
  !==================================================================
  ! Wrapper module for liblbfgsb.a
  !==================================================================

  use utility

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: lbfgs

contains

!====================================================================
! Public
!====================================================================

subroutine lbfgs(n,m,x,f,g,factr,pgtol,task,iprint)

  !--------------------------------------------------------
  ! This subroutine is a wrapper for the library subroutine setulb.
  ! It performs the limited-memory BFGS algorithm
  ! without the boundary conditions.
  !--------------------------------------------------------

  integer,                     intent(IN)    :: n
  integer,                     intent(IN)    :: m
  real(DBL), dimension(n),     intent(INOUT) :: x
  real(DBL),                   intent(INOUT) :: f
  real(DBL), dimension(n),     intent(INOUT) :: g
  real(DBL),                   intent(IN)    :: factr
  real(DBL),                   intent(IN)    :: pgtol
  character(60),               intent(INOUT) :: task
  integer,                     intent(IN)    :: iprint

  integer,   allocatable, dimension(:), save :: nbd
  integer,   allocatable, dimension(:), save :: iwa
  real(DBL), allocatable, dimension(:), save :: l
  real(DBL), allocatable, dimension(:), save :: u
  real(DBL), allocatable, dimension(:), save :: wa
  character(60),                        save :: csave
  logical,   dimension(4),              save :: lsave
  integer,   dimension(44),             save :: isave
  real(DBL), dimension(29),             save :: dsave
  integer                                    :: err_n
  character(120)                             :: err_msg

  ! allocation and initialization -------------------------
  if (.not.allocated(nbd)) then
    allocate(nbd(n), stat=err_n, errmsg=err_msg)
    if (err_n/=0) then
      call error("lbfgs: "//trim(err_msg))
    end if
    nbd=0
  end if

  if (.not.allocated(iwa)) then
    allocate(iwa(3*n), stat=err_n, errmsg=err_msg)
    if (err_n/=0) then
      call error("lbfgs: "//trim(err_msg))
    end if
  end if

  if (.not.allocated(l)) then
    allocate(l(n), stat=err_n, errmsg=err_msg)
    if (err_n/=0) then
      call error("lbfgs: "//trim(err_msg))
    end if
    l=0.0_DBL
  end if

  if (.not.allocated(u)) then
    allocate(u(n), stat=err_n, errmsg=err_msg)
    if (err_n/=0) then
      call error("lbfgs: "//trim(err_msg))
    end if
    u=0.0_DBL
  end if

  if (.not.allocated(wa)) then
    allocate(wa(2*n*m+5*n+11*m*m+8*m), stat=err_n, errmsg=err_msg)
    if (err_n/=0) then
      call error("lbfgs: "//trim(err_msg))
    end if
  end if

  ! call the setulb library subroutine --------------------
  call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,&
    &task,iprint,csave,lsave,isave,dsave)

end subroutine lbfgs

!====================================================================

end module bfgs_wrapper

