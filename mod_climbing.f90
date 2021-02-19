module climbing

  !==================================================================
  !   Climbing Module
  !==================================================================
  !   It contains useful procedures
  ! for the climbing image variant of the NEB (CI-NEB).
  !   If the climbing image operating mode is used,
  ! this method can be applied to an arbitrary number of images
  ! chosen by the user.
  !   It also contains an experimental descending image procedure.
  ! This deals with the minimization of all the minima
  ! encountered along the MEP.
  !   If the descending image operating mode is used,
  ! all the images located in an energy minimum will be descended.
  !==================================================================

  use utility
  use geometry
  use pes
  use elastic

  implicit none
  save
  private

  ! protected variables -----------------------------------
  public    :: flag_climbing_image,         &
               flag_climbing_quick_start,   &
               flag_descending_image,       &
               flag_descending_quick_start, &
               climbing_image_n
  protected :: flag_climbing_image,         &
               flag_climbing_quick_start,   &
               flag_descending_image,       &
               flag_descending_quick_start, &
               climbing_image_n
  ! public procedures -------------------------------------
  public    :: set_climbing_image,          &
               set_climbing_quick_start,    &
               set_descending_image,        &
               set_descending_quick_start,  &
               exec_climbing,               &
               exec_descending

  !--------------------------------------------------------
  logical :: flag_climbing_image         = .false.
  logical :: flag_climbing_quick_start   = .false.
  logical :: flag_descending_image       = .false.
  logical :: flag_descending_quick_start = .false.
  logical :: flag_all_maxima_climbing    = .false.
  integer :: climbing_image_n            = 0

contains

!====================================================================
! Public
!====================================================================

subroutine set_climbing_image(str)

  character(*), intent(INOUT) :: str

  character(*), parameter     :: my_name    = "set_climbing_image"
  logical, save               :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  call tolower(str)

  if (isinteger(trim(adjustl(str)))) then
    read(str,*) climbing_image_n

    if (climbing_image_n<0) then
      call error("set_climbing_image: argument must be a positive integer")
    end if

    if (climbing_image_n==0) then
      flag_climbing_image = .false.
    else
      flag_climbing_image = .true.
    end if
  else
    select case (str)
    case ("y","yes","t","true","on")
      flag_climbing_image = .true.
      climbing_image_n    = 1    
    case ("n","no","f","false","off")
      flag_climbing_image = .false.
      climbing_image_n    = 0    
    case ("all")
      flag_climbing_image      = .true.
      flag_all_maxima_climbing = .true.
    case default
      call error("set_climbing_image: argument """//trim(str)//&
        &""" is not valid")
    end select
  end if

  first_call = .false.

end subroutine set_climbing_image

!====================================================================

subroutine set_climbing_quick_start(flag)

  logical, intent(IN) :: flag

  character(*), parameter     :: my_name    = "set_climbing_quick_start"
  logical, save               :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  flag_climbing_quick_start = flag

  first_call = .false.

end subroutine set_climbing_quick_start

!====================================================================

subroutine set_descending_image(flag)

  logical, intent(IN) :: flag

  character(*), parameter     :: my_name    = "set_descending_image"
  logical, save               :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  flag_descending_image = flag

  first_call = .false.

end subroutine set_descending_image

!====================================================================

subroutine set_descending_quick_start(flag)

  logical, intent(IN) :: flag

  character(*), parameter     :: my_name    = "set_descending_quick_start"
  logical, save               :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  flag_descending_quick_start = flag

  first_call = .false.

end subroutine set_descending_quick_start

!====================================================================

subroutine exec_climbing(fixed,write_output)

  logical,                  intent(IN) :: fixed
  logical,                  intent(IN) :: write_output

  logical,   allocatable, dimension(:) :: mask
  real(DBL), allocatable, dimension(:) :: forces
  character(8)                         :: istr
  integer                              :: i
  integer                              :: indx
  integer                              :: atoms
  integer                              :: tot_maxima
  integer                              :: up_lim
  real(DBL)                            :: dp
  integer                              :: err_n
  character(120)                       :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_climbing_image.eqv..false.) then
    call error("exec_climbing: subroutine called with false climbing image flag")
  end if

  ! allocation section ------------------------------------
  allocate(mask(0:image_n+1),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("exec_climbing: "//trim(err_msg))
  end if

  allocate(forces(geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("exec_climbing: "//trim(err_msg))
  end if

  ! working section ---------------------------------------
  call get_maxima(pes_energy,mask)
  tot_maxima = 0
  do i=0, image_n+1
    if (mask(i)) then
      tot_maxima = tot_maxima+1
    end if
  end do

  if (flag_all_maxima_climbing) then
    up_lim = tot_maxima
  else
    up_lim = min(climbing_image_n,tot_maxima)
  end if

  do i=1, up_lim
    indx   = maxloc(pes_energy,1,mask)-1
    dp     = dot_product(pes_forces(indx,:),norm_tangent(indx,:))
    forces = dp*norm_tangent(indx,:)
    forces = pes_forces(indx,:)-2*forces

    if (fixed) then
      atoms = geom_len/3
      if (atoms>=1) forces(1:3) = 0.0_DBL
      if (atoms>=2) forces(5:6) = 0.0_DBL
      if (atoms>=3) forces(9:9) = 0.0_DBL
    end if
    
    if (write_output) then
      write(istr,'(I8)') indx
      istr = adjustl(istr)
      write(FILEOUT,'(5X,"Climbing on image: ",A)') trim(istr)
    end if

    call set_total_forces_on_i(indx,forces)
    mask(indx) = .false.
  end do

  ! deallocation section ------------------------------------
  deallocate(mask,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("exec_climbing: "//trim(err_msg))
  end if

  deallocate(forces,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("exec_climbing: "//trim(err_msg))
  end if

end subroutine exec_climbing

!====================================================================

subroutine exec_descending(fixed,write_output)

  logical,                  intent(IN) :: fixed
  logical,                  intent(IN) :: write_output

  logical,   allocatable, dimension(:) :: mask
  real(DBL), allocatable, dimension(:) :: forces
  character(8)                         :: istr
  integer                              :: i
  integer                              :: indx
  integer                              :: atoms
  integer                              :: tot_minima
  integer                              :: err_n
  character(120)                       :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_descending_image.eqv..false.) then
    call error("exec_descending: subroutine called with false descending image flag")
  end if

  ! allocation section ------------------------------------
  allocate(mask(0:image_n+1),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("exec_descending: "//trim(err_msg))
  end if

  allocate(forces(geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("exec_descending: "//trim(err_msg))
  end if

  ! working section ---------------------------------------
  call get_minima(pes_energy,mask)
  tot_minima = 0
  do i=0, image_n+1
    if (mask(i)) then
      tot_minima = tot_minima+1
    end if
  end do

  do i=1, tot_minima
    indx   = minloc(pes_energy,1,mask)-1
    forces = pes_forces(indx,:)

    if (fixed) then
      atoms = geom_len/3
      if (atoms>=1) forces(1:3) = 0.0_DBL
      if (atoms>=2) forces(5:6) = 0.0_DBL
      if (atoms>=3) forces(9:9) = 0.0_DBL
    end if
    
    if (write_output) then
      write(istr,'(I8)') indx
      istr = adjustl(istr)
      write(FILEOUT,'(5X,"Descending on image: ",A)') trim(istr)
    end if

    call set_total_forces_on_i(indx,forces)
    mask(indx) = .false.
  end do

  ! deallocation section ----------------------------------
  deallocate(mask,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("exec_descending: "//trim(err_msg))
  end if

  deallocate(forces,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("exec_descending: "//trim(err_msg))
  end if

end subroutine exec_descending

!====================================================================

end module climbing

