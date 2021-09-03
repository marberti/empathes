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

module lbfgs

  use utility

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: init_lbfgs,    &
            lbfgs_internal

  !--------------------------------------------------------
  logical                                :: flag_init_lbfgs = .false.
  integer                                :: lbfgs_memory
  integer                                :: lbfgs_vectors_size
  integer                                :: store_vectors_counter
  real(DBL), dimension(:,:), allocatable :: s_vectors
  real(DBL), dimension(:,:), allocatable :: y_vectors

contains

!====================================================================
! Public
!====================================================================

subroutine init_lbfgs(mem,sz)

  integer,     intent(IN) :: mem
  integer,     intent(IN) :: sz

  character(*), parameter :: my_name = "init_lbfgs"
  integer                 :: err_n
  character(120)          :: err_msg

  ! preliminary check -------------------------------------
  if (flag_init_lbfgs) then
    call error(my_name//": module already initialized")
  end if

  ! arguments checks --------------------------------------
  if (mem <= 0) then
    call error(my_name//": argument mem must be a non-zero positive integer")
  end if

  if (sz <= 0) then
    call error(my_name//": argument sz must be a non-zero positive integer")
  end if

  ! init global variables ---------------------------------
  lbfgs_memory          = mem
  lbfgs_vectors_size    = sz
  store_vectors_counter = 0

  ! allocation section ------------------------------------
  allocate(s_vectors(lbfgs_memory,lbfgs_vectors_size),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(y_vectors(lbfgs_memory,lbfgs_vectors_size),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  flag_init_lbfgs = .true.

end subroutine init_lbfgs

!====================================================================

subroutine lbfgs_internal()

  character(*), parameter :: my_name = "lbfgs_internal"

  ! preliminary check -------------------------------------
  if (flag_init_lbfgs.eqv..false.) then
    call error(my_name//": module not initialized")
  end if

end subroutine lbfgs_internal

!====================================================================
! Private
!====================================================================

subroutine lbfgs_get_direction()

end subroutine lbfgs_get_direction

!====================================================================

subroutine store_vectors(s_vec,y_vec)

  real(DBL), dimension(:), intent(IN) :: s_vec
  real(DBL), dimension(:), intent(IN) :: y_vec

  integer                             :: i

  store_vectors_counter = store_vectors_counter + 1
  i = mod(store_vectors_counter,lbfgs_memory)

  s_vectors(i,:) = s_vec
  y_vectors(i,:) = y_vec

end subroutine store_vectors

!====================================================================

end module lbfgs

