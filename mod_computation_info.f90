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

module computation_info

  use utility
  use geometry
  use rotation
  use idpp
  use pes
  use pes_data
  use elastic
  use climbing
  use optimization

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: write_computation_info

contains

!====================================================================

subroutine write_computation_info()

  ! write all the computation details
  ! related to the current execution.

  write(FILEOUT,'(1X,"INF ",17("info"))')

  ! geometry info -----------------------------------------

  ! rotation info -----------------------------------------
  
  ! idpp info ---------------------------------------------

  ! pes info ----------------------------------------------

  ! pes_data info -----------------------------------------

  ! elastic info ------------------------------------------

  ! climbing info -----------------------------------------

  ! optimization info -------------------------------------

  write(FILEOUT,'(1X,"INF ",17("info"))')

end subroutine write_computation_info

!====================================================================

end module computation_info

