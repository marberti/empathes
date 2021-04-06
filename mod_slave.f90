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

module slave
! start: module wide conditional compilation --------------
#ifdef USE_MPI

  use mpi
  use utility
  use geometry
  use pes

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: slave_on_idle

contains

!====================================================================

subroutine slave_on_idle()

  integer :: cmd
  integer :: err_n

  do
    call mpi_bcast(cmd,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)

    select case(cmd)
    case (MMPI_MSG_END_MAIN_EXEC)
      call end_main_exec() ! slave
    case (MMPI_MSG_INIT_IMAGES)
      call init_images()
    case (MMPI_MSG_INIT_PES_MODULE)
      call init_pes_module()
    case (MMPI_MSG_COMPUTE_PES_FORCES)
      call compute_pes_forces()
    case default
      call error("slave_on_idle: received unknown command")
    end select
  end do

end subroutine slave_on_idle

!====================================================================

! end: module wide conditional compilation ----------------
#endif
end module slave

