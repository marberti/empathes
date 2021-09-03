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
  real(DBL), dimension(:),   allocatable :: rho
  integer,   dimension(:),   allocatable :: sorted_indexes

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

  allocate(rho(lbfgs_memory),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(sorted_indexes(lbfgs_memory),stat=err_n,errmsg=err_msg)
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

subroutine lbfgs_get_direction(df,h0,r)

  ! Based on:
  !
  ! Nocedal - Numerical Optimization - 2nd edition
  ! Algorithm 7.4

  real(DBL), dimension(:),   intent(IN)    :: df
  real(DBL), dimension(:,:), intent(IN)    :: h0
  real(DBL), dimension(:),   intent(INOUT) :: r

  character(*), parameter                  :: my_name = "lbfgs_get_direction"
  real(DBL), dimension(:), allocatable     :: q
  real(DBL), dimension(:), allocatable     :: a
  real(DBL)                                :: b
  integer                                  :: i
  integer                                  :: j
  integer                                  :: err_n
  character(120)                           :: err_msg

  ! allocation section ------------------------------------
  allocate(q(lbfgs_vectors_size),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(a(lbfgs_memory),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! working section ---------------------------------------
  call compute_rho()

  q = df

  do j = 1, min(store_vectors_counter,lbfgs_memory)
    i = sorted_indexes(j)

    a(i) = rho(i) * dot_product(s_vectors(i,:),q)
    q    = q - a(i)*y_vectors(i,:)
  end do

  r = reshape(matmul(h0,reshape(q,(/size(q), 1/))), (/size(r)/))

  do j = min(store_vectors_counter,lbfgs_memory), 1, -1
    i = sorted_indexes(j)

    b = rho(i) * dot_product(y_vectors(i,:), r)
    r = r + (a(i) - b)*s_vectors(i,:)
  end do

  ! deallocation section ----------------------------------
  deallocate(q,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(a,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

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

  call sort_indexes()

end subroutine store_vectors

!====================================================================

subroutine sort_indexes()

  ! This subroutine updates the sorted_indexes array,
  ! that contains the indexes of s_vectors and y_vectors
  ! from the newest to the oldest.

  integer :: indx
  integer :: i

  indx = mod(store_vectors_counter,lbfgs_memory)

  do i = 1, lbfgs_memory
    sorted_indexes(i) = indx

    indx = indx - 1
    if (indx == 0) then
      indx = lbfgs_memory
    end if
  end do

end subroutine sort_indexes

!====================================================================

subroutine compute_rho()

end subroutine compute_rho

!====================================================================

end module lbfgs

