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

