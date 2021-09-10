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
  public :: init_lbfgs,     &
            finalize_lbfgs, &
            lbfgs_internal

  ! drivers -----------------------------------------------
  public :: driver_lbfgs

  !--------------------------------------------------------
  logical                                :: flag_init_lbfgs = .false.
  integer                                :: lbfgs_memory
  integer                                :: lbfgs_vectors_size
  integer                                :: stored_vectors_counter
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
  lbfgs_memory           = mem
  lbfgs_vectors_size     = sz
  stored_vectors_counter = 0

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

  ! init global variables ---------------------------------
  s_vectors       = 0.0_DBL
  y_vectors       = 0.0_DBL
  rho             = 0.0_DBL
  sorted_indexes  = 0

  flag_init_lbfgs = .true.

end subroutine init_lbfgs

!====================================================================

subroutine finalize_lbfgs()

  character(*), parameter :: my_name = "finalize_lbfgs"
  integer                 :: err_n
  character(120)          :: err_msg

  ! preliminary check -------------------------------------
  if (flag_init_lbfgs.eqv..false.) then
    call error(my_name//": module not initialized")
  end if

  ! deallocation section ----------------------------------
  deallocate(s_vectors,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(y_vectors,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(rho,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(sorted_indexes,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! init global variables ---------------------------------
  flag_init_lbfgs = .false.

end subroutine finalize_lbfgs

!====================================================================

subroutine lbfgs_internal(cmdstr,x0,x1,df0,df1,reset_alpha)

  ! BFSG method, based on:
  !
  ! Nocedal - Numerical Optimization - 2nd edition
  ! Algorithm 7.5 (only the statements inside the loop)

  character(120),          intent(INOUT) :: cmdstr
  real(DBL), dimension(:), intent(IN)    :: x0
  real(DBL), dimension(:), intent(INOUT) :: x1
  real(DBL), dimension(:), intent(IN)    :: df0
  real(DBL), dimension(:), intent(IN)    :: df1
  logical, optional,       intent(IN)    :: reset_alpha

  character(*), parameter                :: my_name = "lbfgs_internal"
  logical, parameter                     :: flag_max_displacement = .true.
  real(DBL), parameter                   :: max_displacement = 0.2_DBL
  real(DBL), parameter                   :: alpha0           = 1.0_DBL
  real(DBL), save                        :: alpha            = alpha0
  real(DBL)                              :: p_max
  real(DBL), dimension(:),   allocatable :: p
  real(DBL)                              :: h0
  integer                                :: i
  integer                                :: err_n
  character(120)                         :: err_msg

  ! initialization check ----------------------------------
  if (flag_init_lbfgs.eqv..false.) then
    call error(my_name//": module not initialized")
  end if

  ! check optional argument -------------------------------
  if (present(reset_alpha)) then
    if (reset_alpha.eqv..true.) then
      alpha = alpha0
    end if
  end if

  ! arguments' dimensions check ---------------------------
  if (size(x0) /= lbfgs_vectors_size) then
    call error(my_name//": size of argument x0 differs from initialization size")
  end if

  if (size(x1) /= lbfgs_vectors_size) then
    call error(my_name//": size of argument x1 differs from initialization size")
  end if

  if (size(df0) /= lbfgs_vectors_size) then
    call error(my_name//": size of argument df0 differs from initialization size")
  end if

  if (size(df1) /= lbfgs_vectors_size) then
    call error(my_name//": size of argument df1 differs from initialization size")
  end if

  ! allocation section ------------------------------------
  allocate(p(lbfgs_vectors_size),stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! cmdstr parsing ----------------------------------------
  ! cmdstr = {START, EVALUATE_DF1, EVALUATED, DONE, ERROR}
  select case (cmdstr)
  case ("START")
    ! set h0
    if (stored_vectors_counter > 0) then
      i  = sorted_indexes(1)
      h0 = dot_product(s_vectors(i,:),y_vectors(i,:)) / dot_product(y_vectors(i,:),y_vectors(i,:))
    else
      h0 = 1.0_DBL
    end if

    ! get direction
    call lbfgs_get_direction(df0,h0,p)
    p = -p

    if (flag_max_displacement) then
      alpha = alpha0
      p_max = maxval(abs(p))
      if (alpha*p_max > max_displacement) then
        alpha = max_displacement / p_max
      end if
    end if

    ! new coordinates
    x1 = x0 + alpha*p

    ! need df1 to proceed
    cmdstr = "EVALUATE_DF1"
  case ("EVALUATED")
    ! save the new s and y vectors
    call store_vectors(x1 - x0, df1 - df0)
    cmdstr = "DONE"
  case default
    cmdstr = "ERROR"
    call error(my_name//": unknown cmdstr "//trim(cmdstr))
  end select

  ! deallocation section ----------------------------------
  deallocate(p,stat=err_n,errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
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
  real(DBL),                 intent(IN)    :: h0
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
  q = df

  do j = 1, min(stored_vectors_counter,lbfgs_memory)
    i = sorted_indexes(j)

    a(i) = rho(i) * dot_product(s_vectors(i,:),q)
    q    = q - a(i)*y_vectors(i,:)
  end do

  r = h0 * q

  do j = min(stored_vectors_counter,lbfgs_memory), 1, -1
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

  ! in addition to storing the two input argument vectors,
  ! this subroutine computes the corresponding rho value and sorts the indexes

  real(DBL), dimension(:), intent(IN) :: s_vec
  real(DBL), dimension(:), intent(IN) :: y_vec

  integer                             :: i

  stored_vectors_counter = stored_vectors_counter + 1
  i = mod(stored_vectors_counter,lbfgs_memory)
  if (i == 0) then
    i = lbfgs_memory
  end if

  s_vectors(i,:) = s_vec
  y_vectors(i,:) = y_vec
  rho(i)         = 1.0_DBL / dot_product(s_vec,y_vec)

  call sort_indexes()

end subroutine store_vectors

!====================================================================

subroutine sort_indexes()

  ! This subroutine updates the sorted_indexes array,
  ! that contains the indexes of s_vectors and y_vectors
  ! from the newest to the oldest.

  integer :: indx
  integer :: i

  indx = mod(stored_vectors_counter,lbfgs_memory)
  if (indx == 0) then
    indx = lbfgs_memory
  end if

  do i = 1, lbfgs_memory
    sorted_indexes(i) = indx

    indx = indx - 1
    if (indx == 0) then
      indx = lbfgs_memory
    end if
  end do

end subroutine sort_indexes

!====================================================================
! Drivers
!====================================================================

subroutine driver_lbfgs()

  character(*), parameter          :: my_name    = "driver_lbfgs"
  character(*), parameter          :: sep        = '(80("-"))'
  character(*), parameter          :: vecform    = '(1X,I3,2(3X,F5.2))'  ! length -> vector_len
  character(*), parameter          :: vecform2   = '(1X,A3,2(3X,F5.2))'  ! length -> vector_len
  character(*), parameter          :: indxform   = '(1X,"ind",3(3X,I5))' ! length -> memory
  integer, parameter               :: memory     = 3
  integer, parameter               :: vector_len = 2
  integer, parameter               :: max_cyc_n  = 5
  character(120)                   :: cmdstr
  integer                          :: i
  integer                          :: cyc_n
  character(8)                     :: istr
  real(DBL), dimension(vector_len) :: x0
  real(DBL), dimension(vector_len) :: x1
  real(DBL), dimension(vector_len) :: df0
  real(DBL), dimension(vector_len) :: df1

  write(*,sep)
  write(*,*) my_name//": call inti_lbfgs()"
  call init_lbfgs(memory,vector_len)

  write(istr,'(I8)') memory
  istr = adjustl(istr)
  write(*,*) my_name//": lbfgs_memory       -> ", lbfgs_memory,       "(", trim(istr), ")"
  write(istr,'(I8)') vector_len
  istr = adjustl(istr)
  write(*,*) my_name//": lbfgs_vectors_size -> ", lbfgs_vectors_size, "(", trim(istr), ")"

  x0  = 2.0_DBL
  df0 = 3.0_DBL

  do cyc_n = 1, max_cyc_n
    write(*,sep)
    write(*,*) my_name//": cycle                  -> ", cyc_n

    write(*,*) my_name//": stored_vectors_counter -> ", stored_vectors_counter

    write(*,*) my_name//": s_vectors"
    do i = 1, lbfgs_memory
      write(*,vecform) i, s_vectors(i,:)
    end do

    write(*,*) my_name//": y_vectors"
    do i = 1, lbfgs_memory
      write(*,vecform) i, y_vectors(i,:)
    end do

    write(*,*) my_name//": sorted_indexes"
    write(*,indxform) sorted_indexes

    write(*,*) my_name//":"
    write(*,vecform2) "x0",  x0
    write(*,vecform2) "df0", df0

    write(*,*) my_name//": call lbfgs_internal(), START"
    cmdstr = "START"
    if (cyc_n == 1) then
      call lbfgs_internal(cmdstr,x0,x1,df0,df1,reset_alpha=.true.)
    else
      call lbfgs_internal(cmdstr,x0,x1,df0,df1)
    end if
    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": received command string """//trim(cmdstr)//""", expected ""EVALUATE_DF1""")
    end if
    write(*,*) my_name//": stored_vectors_counter -> ", stored_vectors_counter

    df1 = 0.95_DBL * df0
    write(*,*) my_name//":"
    write(*,vecform2) "x1",  x1
    write(*,vecform2) "df1", df1

    write(*,*) my_name//": call lbfgs_internal(), EVALUATED"
    cmdstr = "EVALUATED"
    call lbfgs_internal(cmdstr,x0,x1,df0,df1)
    if (cmdstr /= "DONE") then
      call error(my_name//": received command string """//trim(cmdstr)//""", expected ""DONE""")
    end if
    write(*,*) my_name//": stored_vectors_counter -> ", stored_vectors_counter

    write(*,*) my_name//": s_vectors"
    do i = 1, lbfgs_memory
      write(*,vecform) i, s_vectors(i,:)
    end do

    write(*,*) my_name//": y_vectors"
    do i = 1, lbfgs_memory
      write(*,vecform) i, y_vectors(i,:)
    end do

    write(*,*) my_name//": sorted_indexes"
    write(*,indxform) sorted_indexes

    x0  = x1
    df0 = df1
  end do

  write(*,sep)
  write(*,*) my_name//": call finalize_lbfgs()"
  call finalize_lbfgs()

  write(*,sep)

end subroutine driver_lbfgs

!====================================================================

end module lbfgs

