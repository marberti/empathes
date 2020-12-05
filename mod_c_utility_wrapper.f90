module c_utility_wrapper

  use, intrinsic :: iso_c_binding
  use utility

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: f_chdir, &
            f_error

contains

!====================================================================

subroutine f_chdir(dir)

  character(*),                    intent(IN) :: dir

  integer, parameter                          :: cdir_len = 120
  character(kind=C_CHAR), dimension(cdir_len) :: cdir
  integer                                     :: i
  integer                                     :: loop_end
  character(200)                              :: err_msg

  interface
    subroutine c_chdir(cdir) bind(c)
      use, intrinsic :: iso_c_binding
      character(kind=C_CHAR), dimension(120) :: cdir
    end subroutine c_chdir
  end interface

  ! check string lenght -----------------------------------
  if (len_trim(dir)>=cdir_len) then
    err_msg = "f_chdir: too long argument """//trim(dir)//""""
    call error(trim(err_msg))
  end if

  ! copy Fortran string into C characters array -----------
  loop_end = len_trim(dir)+1
  do i=1, loop_end
    if (i==loop_end) then
      cdir(i) = C_NULL_CHAR
    else
      cdir(i) = dir(i:i)
    end if
  end do

  ! call c function ---------------------------------------
  write(*,*) "f_chdir: changing directory to """,trim(dir),""""
  call c_chdir(cdir)

end subroutine f_chdir

!====================================================================

subroutine f_error(msg) bind(c)

  character(kind=C_CHAR), dimension(300), intent(IN) :: msg

  integer, parameter                                 :: err_msg_len = 300
  integer                                            :: i
  integer                                            :: j
  character(err_msg_len)                             :: err_msg

  ! copy C characters array into Fortran string -----------
  do i=1, err_msg_len
    if (msg(i)==C_NULL_CHAR) then
      j = i
      exit
    end if

    err_msg(i:i) = msg(i)
  end do

  do i=j, err_msg_len
    err_msg(i:i) = " "
  end do

  ! call error --------------------------------------------
  call error(trim(err_msg))

end subroutine f_error

!====================================================================

end module c_utility_wrapper

