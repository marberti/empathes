module pes_data

  use utility

  implicit none
  save
  private

  ! protected flags ---------------------------------------
  public    :: flag_pesd_scfcycle, flag_pesd_scfconv,&
    &flag_pesd_scfvshift, flag_pesd_intgrid, flag_pesd_additional_cmd,&
    &flag_pesd_auxiliary_files
  protected :: flag_pesd_scfcycle,&
               flag_pesd_scfconv,&
               flag_pesd_scfvshift,&
               flag_pesd_intgrid,&
               flag_pesd_additional_cmd,&
               flag_pesd_auxiliary_files
  ! protected variables -----------------------------------
  public    :: pesd_scf_flag_count, pesd_scfcycle, pesd_scfconv,&
    &pesd_scfvshift, pesd_intgrid, pesd_additional_cmd,&
    &pesd_auxiliary_files, pesd_auxiliary_files_n
  protected :: pesd_scf_flag_count,&
               pesd_scfcycle,&
               pesd_scfconv,&
               pesd_scfvshift,&
               pesd_intgrid,&
               pesd_additional_cmd,&
               pesd_auxiliary_files,&
               pesd_auxiliary_files_n
  ! public procedures -------------------------------------
  public :: set_pesd_scfcycle, set_pesd_scfconv, set_pesd_scfvshift,&
    &set_pesd_intgrid, set_pesd_additional_cmd, set_pesd_auxiliary_files

  !--------------------------------------------------------
  logical :: flag_pesd_scfcycle        = .false.
  logical :: flag_pesd_scfconv         = .false.
  logical :: flag_pesd_scfvshift       = .false.
  logical :: flag_pesd_intgrid         = .false.
  logical :: flag_pesd_additional_cmd  = .false.
  logical :: flag_pesd_auxiliary_files = .false.
  integer :: pesd_scf_flag_count=0 ! automatically incremented
                                   ! by set_pesd_scf* subroutines

  integer :: pesd_scfcycle
  integer :: pesd_scfvshift
  real(DBL) :: pesd_scfconv
  character(120) :: pesd_intgrid
  character(120) :: pesd_additional_cmd
  character(120), allocatable, dimension(:) :: pesd_auxiliary_files
  integer :: pesd_auxiliary_files_n

contains

!====================================================================
! Public
!====================================================================

subroutine set_pesd_scfcycle(str)

  character(*), intent(IN) :: str

  if (flag_pesd_scfcycle) then
    call error("set_pesd_scfcycle: scf cycles already setted")
  end if

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

  character(*), intent(IN) :: str

  if (flag_pesd_scfconv) then
    call error("set_pesd_scfconv: scf convergence already setted")
  end if

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

  character(*), intent(IN) :: str

  if (flag_pesd_scfvshift) then
    call error("set_pesd_scfvshift: scf vshift already setted")
  end if

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

  character(*), intent(IN) :: str

  if (flag_pesd_intgrid) then
    call error("set_pesd_intgrid: integration grid already setted")
  end if

  if (len_trim(str)==0) then
    call error("set_pesd_intgrid: empty argument")
  end if

  pesd_intgrid=trim(adjustl(str))

  flag_pesd_intgrid=.true.

end subroutine set_pesd_intgrid

!====================================================================

subroutine set_pesd_additional_cmd(str)

  character(*), intent(IN) :: str

  if (flag_pesd_additional_cmd) then
    call error("set_pesd_additional_cmd: additional command already setted")
  end if

  if (len_trim(str)==0) then
    call error("set_pesd_additional_cmd: empty argument")
  end if

  pesd_additional_cmd=trim(adjustl(str))

  flag_pesd_additional_cmd=.true.

end subroutine set_pesd_additional_cmd

!====================================================================

subroutine set_pesd_auxiliary_files(str)

  ! read the entire input file string of the kind:
  ! #KEYWORD int_N str_1 str_2 ... str_N

  character(*), intent(IN) :: str
  character(120) :: field
  character(8)   :: i_str
  integer :: i
  integer :: err_n
  character(120) :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_pesd_auxiliary_files) then
    call error("set_pesd_auxiliary_files: auxiliary files already setted")
  end if

  ! get number of aux files -------------------------------
  call get_field(str,field,2,err_n,err_msg)
  if (err_n/=0) then
    call error("set_pesd_auxiliary_files: "//trim(err_msg))
  end if

  if (isinteger(trim(adjustl(field)))) then
    read(field,*) pesd_auxiliary_files_n
  else
    call error("set_pesd_auxiliary_files: wrong argument """&
      &//trim(field)//""", expected integer")
  end if

  if (pesd_auxiliary_files_n<=0) then
    call error("set_pesd_auxiliary_files: argument must be a non-zero positive integer")
  end if

  ! allocation section ------------------------------------
  allocate(pesd_auxiliary_files(pesd_auxiliary_files_n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("set_pesd_auxiliary_files: "//trim(err_msg))
  end if

  ! get aux files -----------------------------------------
  do i=1, pesd_auxiliary_files_n
    call get_field(str,pesd_auxiliary_files(i),i+2,err_n,err_msg)
    if (err_n/=0) then
      call error("set_pesd_auxiliary_files: "//trim(err_msg))
    end if
  end do

  ! check that no more files are specified ----------------
  call get_field(str,field,pesd_auxiliary_files_n+3,err_n,err_msg)
  if (err_n==0) then ! yes, it's correct ==
    write(i_str,'(I8)') pesd_auxiliary_files_n
    i_str=adjustl(i_str)
    call error("set_pesd_auxiliary_files: specified more than "&
      &//trim(i_str)//" auxiliary file(s)")
  end if

  ! finalize ----------------------------------------------
  flag_pesd_auxiliary_files=.true.

end subroutine set_pesd_auxiliary_files

!====================================================================

end module pes_data

