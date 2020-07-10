module pes_data

  use utility

  implicit none
  save
  private

  ! protected flags ---------------------------------------
  public    :: flag_pesd_scfcycle, flag_pesd_scfconv,&
    &flag_pesd_scfvshift, flag_pesd_intgrid, flag_pesd_additional_cmd 
  protected :: flag_pesd_scfcycle,&
               flag_pesd_scfconv,&
               flag_pesd_scfvshift,&
               flag_pesd_intgrid,&
               flag_pesd_additional_cmd 
  ! protected variables -----------------------------------
  public    :: pesd_scf_flag_count, pesd_scfcycle, pesd_scfconv,&
    &pesd_scfvshift, pesd_intgrid, pesd_additional_cmd
  protected :: pesd_scf_flag_count,&
               pesd_scfcycle,&
               pesd_scfconv,&
               pesd_scfvshift,&
               pesd_intgrid,&
               pesd_additional_cmd
  ! public procedures -------------------------------------
  public :: set_pesd_scfcycle, set_pesd_scfconv, set_pesd_scfvshift,&
    &set_pesd_intgrid, set_pesd_additional_cmd

  !--------------------------------------------------------
  logical :: flag_pesd_scfcycle       = .false.
  logical :: flag_pesd_scfconv        = .false.
  logical :: flag_pesd_scfvshift      = .false.
  logical :: flag_pesd_intgrid        = .false.
  logical :: flag_pesd_additional_cmd = .false.
  integer :: pesd_scf_flag_count=0 ! automatically incremented by set_* subroutines

  integer :: pesd_scfcycle
  integer :: pesd_scfvshift
  real(DBL) :: pesd_scfconv
  character(120) :: pesd_intgrid
  character(120) :: pesd_additional_cmd

contains

!====================================================================
! Public
!====================================================================

subroutine set_pesd_scfcycle(str)

  character(*) :: str

  if (len_trim(str)==0) then
    call error("set_pesd_scfcycle: empty argument")
  end if

  if (.not.isinteger(trim(adjustl(str)))) then
    call error("set_pesd_scfcycle: wrong argument type")
  end if

  read(str,*) pesd_scfcycle

  if (pesd_scfcycle<0) then
    call error("set_pesd_scfcycle: argument must be a non-zero positive integer")
  end if

  pesd_scf_flag_count=pesd_scf_flag_count+1
  flag_pesd_scfcycle=.true.

end subroutine set_pesd_scfcycle

!====================================================================

subroutine set_pesd_scfconv(str)

  character(*) :: str

  if (len_trim(str)==0) then
    call error("set_pesd_scfconv: empty argument")
  end if

  if (.not.isreal(trim(adjustl(str)))) then
    call error("set_pesd_scfconv: wrong argument type")
  end if

  read(str,*) pesd_scfconv

  pesd_scf_flag_count=pesd_scf_flag_count+1
  flag_pesd_scfconv=.true.

end subroutine set_pesd_scfconv

!====================================================================

subroutine set_pesd_scfvshift(str)

  character(*) :: str

  if (len_trim(str)==0) then
    call error("set_pesd_scfvshift: empty argument")
  end if

  if (.not.isinteger(trim(adjustl(str)))) then
    call error("set_pesd_scfvshift: wrong argument type")
  end if

  read(str,*) pesd_scfvshift

  if (pesd_scfvshift<0) then
    call error("set_pesd_scfvshift: argument must be a non-zero positive integer")
  end if

  pesd_scf_flag_count=pesd_scf_flag_count+1
  flag_pesd_scfvshift=.true.

end subroutine set_pesd_scfvshift

!====================================================================

subroutine set_pesd_intgrid(str)

  character(*) :: str

  if (len_trim(str)==0) then
    call error("set_pesd_intgrid: empty argument")
  end if

  pesd_intgrid=trim(adjustl(str))

  flag_pesd_intgrid=.true.

end subroutine set_pesd_intgrid

!====================================================================

subroutine set_pesd_additional_cmd(str)

  character(*) :: str

  if (len_trim(str)==0) then
    call error("set_pesd_additional_cmd: empty argument")
  end if

  pesd_additional_cmd=trim(adjustl(str))

  flag_pesd_additional_cmd=.true.

end subroutine set_pesd_additional_cmd

!====================================================================

end module pes_data

