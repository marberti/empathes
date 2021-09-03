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
  logical :: flag_init_lbfgs = .false.

contains

!====================================================================
! Public
!====================================================================

subroutine init_lbfgs()

  character(*), parameter :: my_name = "init_lbfgs"

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

end module lbfgs

