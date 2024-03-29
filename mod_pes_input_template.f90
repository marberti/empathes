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

module pes_input_template

#ifdef USE_MPI
  use mpi
#endif
  use utility

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: read_pes_it,  &
            write_pes_it, &
            get_pes_it_n
#ifdef USE_MPI
  public :: mmpi_sync_pes_it
#endif

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

#ifdef USE_MPI
subroutine mmpi_sync_pes_it()

  integer        :: i
  character(8)   :: i_str
  integer        :: lines
  integer        :: err_n
  character(120) :: err_msg

  ! both master and slaves are here -----------------------

  ! pes_it_len bcast --------------------------------------
  call mpi_bcast(pes_it_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)

  if (pes_it_len==0) then
    return
  end if

  ! slaves allocate pes_it --------------------------------
  if (proc_id/=0) then
    allocate (pes_it(pes_it_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(i_str,'(I8)') proc_id
      i_str = adjustl(i_str)
      call error("mmpi_init_pes_module: process "//&
        &trim(i_str)//": "//trim(err_msg))
    end if
  end if

  ! master sends pes_it data ------------------------------
  do i=1, pes_it_len
    ! bcast n and lines
    call mpi_bcast(pes_it(i)%n,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)
    call mpi_bcast(pes_it(i)%lines,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)

    lines = pes_it(i)%lines

    ! slaves allocate pes_it buffer strings
    if (proc_id/=0) then
      allocate (pes_it(i)%s(lines),stat=err_n,errmsg=err_msg)
      if (err_n/=0) then
        write(i_str,'(I8)') proc_id
        i_str = adjustl(i_str)
        call error("mmpi_init_pes_module: process "//&
          &trim(i_str)//": "//trim(err_msg))
      end if
    end if

    ! bcast buffer strings
    call mpi_bcast(pes_it(i)%s,lines*IT_STR_LEN,&
      &MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
  end do

end subroutine mmpi_sync_pes_it
#endif

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

