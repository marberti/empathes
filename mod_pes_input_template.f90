module pes_input_template

  use utility

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: read_pes_it,  &
            write_pes_it, &
            get_pes_it_n

  !--------------------------------------------------------
  integer, parameter                                 :: IT_STR_LEN = 200

  type                                               :: pes_it_t
    integer                                          :: n
    integer                                          :: lines
    character(IT_STR_LEN), allocatable, dimension(:) :: s
  end type pes_it_t

  type (pes_it_t), allocatable, dimension(:)         :: pes_it
  integer                                            :: pes_it_len = 0

contains

!====================================================================
! Public
!====================================================================

subroutine read_pes_it(n_str,fnumb,ending)

  ! Reads the #PESINPUTTEMPLATE block

  character(*), intent(IN) :: n_str
  integer,      intent(IN) :: fnumb
  character(*), intent(IN) :: ending

  character(*), parameter  :: my_name = "read_pes_it: "
  integer                  :: n
  integer                  :: i
  integer                  :: lines
  character(IT_STR_LEN)    :: str
  integer                  :: err_n
  character(120)           :: err_msg

  ! preliminary checks ------------------------------------
  if (isinteger(trim(adjustl(n_str)))) then
    read(n_str,*) n
  else
    call error(my_name//"expected integer after #PESINPUTTEMPLATE, got """//trim(n_str)//"""")
  end if

  if (n<=0) then
    call error(my_name//"argument must be a non-zero positive integer")
  end if

  ! get number of lines in the block ----------------------
  lines = get_lines(fnumb,ending)
  if (lines==0) then
    call error(my_name//"block #PESINPUTTEMPLATE "//trim(adjustl(n_str))//" is empty")
  end if

  ! add another element in pes_it array -------------------
  call add_pes_it(n,lines)
  
  ! rewind the file ---------------------------------------
  do i=1, lines
    backspace(unit=fnumb,iostat=err_n,iomsg=err_msg)
    if (err_n/=0) then
      call error(my_name//trim(err_msg))
    end if
  end do

  ! copy the block ----------------------------------------
  do i=1, lines
    read(fnumb,'(A200)',iostat=err_n,iomsg=err_msg) str
    if (err_n/=0) then
      call error(my_name//trim(err_msg))
    end if
    pes_it(pes_it_len)%s(i) = str
  end do

end subroutine read_pes_it

!====================================================================

subroutine write_pes_it(fnumb,n)

  integer,     intent(IN) :: fnumb
  integer,     intent(IN) :: n

  character(*), parameter :: my_name = "write_pes_it: "
  character(8)            :: i_str
  integer                 :: i
  integer                 :: indx
  logical                 :: is_open

  ! preliminary checks ------------------------------------
  indx = get_pes_it_n(n)
  if (indx==0) then
    write(i_str,'(I8)') n
    i_str = adjustl(i_str)
    call error(my_name//"input template number "//&
      &trim(i_str)//" not specified")
  end if

  inquire(unit=fnumb,opened=is_open)

  if (is_open.eqv..false.) then
    call error(my_name//"output file was not opened")
  end if

  ! write -------------------------------------------------
  do i=1, pes_it(indx)%lines
    write(fnumb,'(A)') trim(pes_it(indx)%s(i))
  end do

end subroutine write_pes_it

!====================================================================

integer function get_pes_it_n(n)

  ! return 0 if n is not found, a positive integer otherwise

  integer, intent(IN) :: n

  integer             :: i

  get_pes_it_n = 0

  do i=1, pes_it_len
    if (pes_it(i)%n==n) then
      get_pes_it_n = i
      return
    end if
  end do

end function get_pes_it_n

!====================================================================
! Private
!====================================================================

subroutine add_pes_it(n,lines)

  integer,                       intent(IN) :: n
  integer,                       intent(IN) :: lines

  character(*), parameter                   :: my_name = "add_pes_it: "
  type(pes_it_t), allocatable, dimension(:) :: tmp_it
  character(8)                              :: i_str
  integer                                   :: indx
  integer                                   :: err_n
  character(120)                            :: err_msg

  ! preliminary checks ------------------------------------
  indx = get_pes_it_n(n)
  if (indx/=0) then
    write(i_str,'(I8)') n
    i_str = adjustl(i_str)
    call error(my_name//"input template number "//&
      &trim(i_str)//" already specified")
  end if

  ! logic -------------------------------------------------
  pes_it_len = pes_it_len+1

  if (pes_it_len>1) then
    call move_alloc(pes_it,tmp_it)
  end if

  allocate(pes_it(pes_it_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//trim(err_msg))
  end if

  if (pes_it_len>1) then
    pes_it(:pes_it_len-1) = tmp_it
  end if

  pes_it(pes_it_len)%n     = n
  pes_it(pes_it_len)%lines = lines

  allocate(pes_it(pes_it_len)%s(lines),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//trim(err_msg))
  end if

end subroutine add_pes_it

!====================================================================

end module pes_input_template

