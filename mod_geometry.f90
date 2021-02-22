module geometry

#ifdef USE_MPI
  use mpi
#endif
  use utility
  use rotation

  implicit none
  save
  private

  ! public variables --------------------------------------
  public    :: start_geom,              &
               end_geom
  ! protected variables -----------------------------------
  public    :: element,                 &
               elabel,                  &
               geom_len,                &
               geom_charge,             &
               geom_multip,             &
               image_n,                 &
               image_geom,              &
               flag_init_images,        &
               flag_only_interpolation, &
               flag_geom_charge,        &
               flag_geom_multip
  protected :: element,                 &
               elabel,                  &
               geom_len,                &
               geom_charge,             &
               geom_multip,             &
               image_n,                 &
               image_geom,              &
               flag_init_images,        &
               flag_only_interpolation, &
               flag_geom_charge,        &
               flag_geom_multip
  ! public procedures -------------------------------------
  public    :: allocate_geom,           &
               allocate_image_geom,     &
               set_geom_len,            &
               set_image_n,             &
               init_images,             &
               update_images,           &
               update_image_geom,       &
               init_element,            &
               update_element,          &
               init_elabel,             &
               update_elabel,           &
               set_geom_charge,         &
               set_geom_multip,         &
               set_only_interpolation,  &
               set_geometries_file

  !--------------------------------------------------------
  logical                                   :: flag_init_images        = .false.
  logical                                   :: flag_only_interpolation = .false.
  logical                                   :: flag_geom_charge        = .false.
  logical                                   :: flag_geom_multip        = .false.
  logical                                   :: flag_init_element       = .false.
  logical                                   :: flag_init_elabel        = .false.
  logical                                   :: flag_set_image_n        = .false.
  logical                                   :: flag_set_geom_len       = .false.
  logical                                   :: flag_geometries_file    = .false.
  integer                                   :: geom_len   ! length of start and end geometries = 3*atoms
  integer                                   :: image_n    ! number of images
  integer                                   :: geom_charge
  integer                                   :: geom_multip
  character(3), allocatable, dimension(:)   :: element
  character(3), allocatable, dimension(:)   :: elabel
  real(DBL),    allocatable, dimension(:)   :: start_geom
  real(DBL),    allocatable, dimension(:)   :: end_geom
  real(DBL),    allocatable, dimension(:,:) :: image_geom ! geometry for each image (image,coordinate)

contains

!====================================================================
! Public
!====================================================================

subroutine allocate_geom(geom,len)

  character(*), intent(IN) :: geom
  integer,      intent(IN) :: len

  integer                  :: err_n
  character(120)           :: err_msg

  ! preliminary checks ------------------------------------
  if (len<=0) then
    call error("allocate_geom: argument ""len"" must be a non-zero positive integer")
  end if

  ! allocation section ------------------------------------
  if (geom=="#START") then
    allocate(start_geom(len),stat=err_n,errmsg=err_msg)
  else if (geom=="#END") then
    allocate(end_geom(len),stat=err_n,errmsg=err_msg)
  else
    call error("allocate_geom: unknown argument "//trim(geom))
  end if

  if (err_n/=0) then
    call error("allocate_geom: "//trim(err_msg))
  end if

end subroutine allocate_geom

!====================================================================

subroutine allocate_image_geom()

  integer        :: err_n
  character(120) :: err_msg

  ! preliminary checks ------------------------------------
  if (allocated(image_geom)) then
    call error("allocate_image_geom: image_geom already allocated")
  end if

  if (flag_set_image_n.eqv..false.) then
    call error("allocate_image_geom: images number not setted")
  end if

  if (flag_set_geom_len.eqv..false.) then
    call error("allocate_image_geom: geometry length not setted")
  end if

  ! allocation section ------------------------------------
  allocate(image_geom(0:image_n+1,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("allocate_image_geom: "//trim(err_msg))
  end if

end subroutine allocate_image_geom

!====================================================================

subroutine set_geom_len(atoms)

  integer, intent(IN) :: atoms

  ! preliminary checks ------------------------------------
  if (flag_set_geom_len) then
    call error("set_geom_len: geometry length already setted")
  end if

  if (atoms<=0) then
    call error("set_geom_len: argument must be a non-zero positive integer")
  end if

  ! set value ---------------------------------------------
  geom_len = 3*atoms

  flag_set_geom_len = .true.

end subroutine set_geom_len

!====================================================================

subroutine set_image_n(str)

  character(*), intent(IN) :: str

  ! preliminary checks ------------------------------------
  if (flag_set_image_n) then
    call error("set_image_n: images number already setted")
  end if

  ! read value --------------------------------------------
  if (isinteger(trim(adjustl(str)))) then
    read(str,*) image_n
  else
    call error("set_image_n: argument """//str//""" is not valid")
  end if

  if (image_n<=0) then
    call error("set_image_n: argument must be a non-zero positive integer")
  end if

  flag_set_image_n = .true.

end subroutine set_image_n

!====================================================================

subroutine init_images()

  !--------------------------------------------------------
  ! Initializes the images geometries.
  ! If a geometries file has been supplied,
  ! this subroutine does some checks and sets flag_init_images,
  ! otherwise it does some kind of initialization.
  ! Currently only linear interpolation is supported.
  !--------------------------------------------------------

  integer                              :: i
  real(DBL), allocatable, dimension(:) :: tmp_start
  real(DBL), allocatable, dimension(:) :: tmp_end
  real(DBL), allocatable, dimension(:) :: delta
  integer                              :: err_n
  character(120)                       :: err_msg

  ! if I'm a slave, i go to... ----------------------------
#ifdef USE_MPI
  if (proc_id/=0) then
    call mmpi_init_images()
    return
  end if
#endif

  ! else, if I'm the master... ----------------------------

  ! conditions checking -----------------------------------
  if (flag_init_images) then
    call error("init_images: images already initialized")
  end if

  if (flag_geometries_file) then
    !TODO do some checks
  else
    ! allocation section ----------------------------------
    call allocate_image_geom()

    if (flag_rotation) then
      allocate(tmp_start(geom_len),stat=err_n,errmsg=err_msg)
      if (err_n/=0) then
        call error("init_images: "//trim(err_msg))
      end if
    end if

    if (flag_rotation) then
      allocate(tmp_end(geom_len),stat=err_n,errmsg=err_msg)
      if (err_n/=0) then
        call error("init_images: "//trim(err_msg))
      end if
    end if

    allocate(delta(geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("init_images: "//trim(err_msg))
    end if

    ! initial geometries rotation -------------------------
    if (flag_rotation) then
      tmp_start = start_geom
      tmp_end   = end_geom
      call rotate_geometry(tmp_start,start_geom)
      call rotate_geometry(tmp_end,end_geom)
    end if

    ! initialization --------------------------------------
    image_geom(0,:)         = start_geom
    image_geom(image_n+1,:) = end_geom

    ! images generated by linear interpolation ------------
    delta = end_geom-start_geom
    do i=1, image_n
      image_geom(i,:) = start_geom+((i/real(image_n+1,DBL))*delta)
    end do

    ! deallocation section --------------------------------
    if (allocated(tmp_start)) then
      deallocate(tmp_start,stat=err_n,errmsg=err_msg)
      if (err_n/=0) then
        call error("init_images: "//trim(err_msg))
      end if
    end if

    if (allocated(tmp_end)) then
      deallocate(tmp_end,stat=err_n,errmsg=err_msg)
      if (err_n/=0) then
        call error("init_images: "//trim(err_msg))
      end if
    end if

    deallocate(delta,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("init_images: "//trim(err_msg))
    end if
  end if

  flag_init_images = .true.

  ! Master finished its initialization. -------------------
  ! If slave processes exist,
  ! master initializes them.
#ifdef USE_MPI
  if (comm_sz>1) then
    call mmpi_init_images()
  end if
#endif

end subroutine init_images

!====================================================================

subroutine update_images(arr)

  !--------------------------------------------------------
  ! Takes an array arr as input and save its values in
  ! image_geom array.
  ! arr must be shaped as (image_n,geom_len)
  !--------------------------------------------------------

  real(DBL), dimension(:,:), intent(IN) :: arr

  ! preliminary checks ------------------------------------
  if (flag_init_images.eqv..false.) then
    call error("update_images: images not initialized")
  end if

  if ((size(arr,1)/=image_n).or.(size(arr,2)/=geom_len)) then
    call error("update_images: wrong argument size")
  end if

  ! update ------------------------------------------------
  image_geom(1:image_n,:) = arr

end subroutine update_images

!====================================================================

subroutine update_image_geom(i,arr)

  integer,                 intent(IN) :: i
  real(DBL), dimension(:), intent(IN) :: arr

  ! preliminary checks ------------------------------------
  if (.not.allocated(image_geom)) then
    call error("update_image_geom: image_geom not allocated")
  end if

  if ((i<0).or.(i>(image_n+1))) then
    call error("update_image_geom: argument ""i"" out of bounds")
  end if

  if (size(arr,1)/=geom_len) then
    call error("update_image_geom: wrong argument ""arr"" size")
  end if

  ! update ------------------------------------------------
  image_geom(i,:) = arr

end subroutine update_image_geom

!====================================================================

subroutine init_element(n)

  integer, intent(IN) :: n

  integer             :: err_n
  character(120)      :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_init_element) then
    call error("init_element: element already initialized")
  end if

  if (n<=0) then
    call error("init_element: argument must be a non-zero positive integer")
  end if

  ! allocation section ------------------------------------
  allocate(element(n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_element: "//trim(err_msg))
  end if

  flag_init_element = .true.

end subroutine init_element

!====================================================================

subroutine update_element(arr)

  character(3), dimension(:), intent(IN) :: arr

  ! preliminary checks ------------------------------------
  if (.not.allocated(element)) then
    call error("update_element: array element not allocated")
  end if

  if (size(arr,1)/=size(element,1)) then
    call error("update_element: wrong argument size")
  end if

  ! update ------------------------------------------------
  element = arr

end subroutine update_element

!====================================================================

subroutine init_elabel(n)

  integer, intent(IN) :: n

  integer             :: err_n
  character(120)      :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_init_elabel) then
    call error("init_elabel: element label already initialized")
  end if

  if (n<=0) then
    call error("init_elabel: argument must be a non-zero positive integer")
  end if

  ! allocation section ------------------------------------
  allocate(elabel(n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_elabel: "//trim(err_msg))
  end if

  flag_init_elabel = .true.

end subroutine init_elabel

!====================================================================

subroutine update_elabel(arr)

  character(3), dimension(:), intent(IN) :: arr

  ! preliminary checks ------------------------------------
  if (.not.allocated(elabel)) then
    call error("update_elabel: array elabel not allocated")
  end if

  if (size(arr,1)/=size(elabel,1)) then
    call error("update_elabel: wrong argument size")
  end if

  ! update ------------------------------------------------
  elabel = arr

end subroutine update_elabel

!====================================================================

subroutine set_geom_charge(str)
  
  character(*), intent(IN) :: str

  logical, save            :: first_call = .true.

  ! preliminary checks ------------------------------------
  if (first_call.eqv..false.) then
    call error("set_geom_charge: subroutine called more than once")
  end if

  ! reading value -----------------------------------------
  if (isinteger(trim(adjustl(str)))) then
    read(str,*) geom_charge
  else
    call error("set_geom_charge: argument """//str//""" is not valid")
  end if

  flag_geom_charge = .true.

  first_call = .false.

end subroutine set_geom_charge

!====================================================================

subroutine set_geom_multip(str)

  character(*), intent(IN) :: str

  logical, save            :: first_call = .true.

  ! preliminary checks ------------------------------------
  if (first_call.eqv..false.) then
    call error("set_geom_multip: subroutine called more than once")
  end if

  ! reading value -----------------------------------------
  if (isinteger(trim(adjustl(str)))) then
    read(str,*) geom_multip
  else
    call error("set_geom_multip: argument """//str//""" is not valid")
  end if

  if (geom_multip<=0) then
    call error("set_geom_multip: argument must be a non-zero positive integer")
  end if

  flag_geom_multip = .true.

  first_call = .false.

end subroutine set_geom_multip

!====================================================================

subroutine set_only_interpolation(flag)

  logical,     intent(IN) :: flag

  character(*), parameter :: my_name    = "set_only_interpolation"
  logical, save           :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  flag_only_interpolation = flag

  first_call = .false.

end subroutine set_only_interpolation

!====================================================================

subroutine set_geometries_file(flag)

  logical,     intent(IN) :: flag

  character(*), parameter :: my_name    = "set_geometries_file"
  logical, save           :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  flag_geometries_file = flag

  first_call = .false.

end subroutine set_geometries_file

!====================================================================
! Private MPI
!====================================================================

#ifdef USE_MPI
subroutine mmpi_init_images()

  !--------------------------------------------------------
  ! Init all slave processes
  !--------------------------------------------------------

  character(8)                            :: istr
  integer                                 :: cmd
  integer                                 :: sz
  integer                                 :: i
  character(30)                           :: str30
  character(3), allocatable, dimension(:) :: e_buff
  logical                                 :: flag
  integer                                 :: err_n
  character(120)                          :: err_msg

  ! master gets slaves into this subroutine ---------------
  if (proc_id==0) then
    cmd = MMPI_MSG_INIT_IMAGES
    call mpi_bcast(cmd,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)
  end if

  if (proc_id==0) then ! master stuffs --------------------
    ! mandatory variables ---------------------------------
    ! bcast image_n
    write(str30,*) image_n
    str30 = adjustl(str30)
    call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)

    ! bcast geom_len
    i = geom_len/3
    call mpi_bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)

    ! bcast element
    sz = geom_len/3
    allocate(e_buff(sz),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr    = adjustl(istr)
      err_msg = "mmpi_init_images: process "//trim(istr)//&
        &": "//err_msg
      call error(err_msg)
    end if

    e_buff = element
    call mpi_bcast(e_buff,sz*len(e_buff),&
      &MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)

    ! optional variables ----------------------------------
    ! bcast geom_charge
    flag = flag_geom_charge
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str30,*) geom_charge
      str30 = adjustl(str30)
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast geom_multip
    flag = flag_geom_multip
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      write(str30,*) geom_multip
      str30 = adjustl(str30)
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! bcast elabel
    flag = flag_init_elabel
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      e_buff = elabel
      call mpi_bcast(e_buff,sz*len(e_buff),&
        &MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    end if

    ! finalize --------------------------------------------
    deallocate(e_buff,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr    = adjustl(istr)
      err_msg = "mmpi_init_images: process "//trim(istr)//&
        &": "//err_msg
      call error(err_msg)
    end if

!DEBUG: used to check if message passing was successful
!    call end_main_exec() ! MPI Debug

  else                 ! slaves stuffs --------------------
    ! preliminary checks ----------------------------------
    if (flag_init_images) then 
      write(istr,'(I8)') proc_id
      istr    = adjustl(istr)
      err_msg = "mmpi_init_images: process "//trim(istr)//&
        &": module geometry already initialized"
      call error(err_msg)
    end if

    ! mandatory variables ---------------------------------
    ! bcast image_n
    call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    call set_image_n(str30)

    ! bcast geom_len
    call mpi_bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)
    call set_geom_len(i)

    ! bcast element
    sz = geom_len/3
    allocate(e_buff(sz),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr    = adjustl(istr)
      err_msg = "mmpi_init_images: process "//trim(istr)//&
        &": "//err_msg
      call error(err_msg)
    end if

    call mpi_bcast(e_buff,sz*len(e_buff),&
      &MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
    call init_element(sz)
    call update_element(e_buff)

    ! optional variables ----------------------------------
    ! bcast geom_charge
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_geom_charge(str30)
    end if

    ! bcast geom_multip
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(str30,len(str30),MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call set_geom_multip(str30)
    end if

    ! bcast elabel
    call mpi_bcast(flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,err_n)
    if (flag) then
      call mpi_bcast(e_buff,sz*len(e_buff),&
        &MPI_CHARACTER,0,MPI_COMM_WORLD,err_n)
      call init_elabel(sz)
      call update_elabel(e_buff)
    end if

    ! allocation section ----------------------------------
    call allocate_image_geom()

    ! finalize --------------------------------------------
    deallocate(e_buff,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      write(istr,'(I8)') proc_id
      istr    = adjustl(istr)
      err_msg = "mmpi_init_images: process "//trim(istr)//&
        &": "//err_msg
      call error(err_msg)
    end if

    flag_init_images = .true.

!DEBUG: used to check if message passing was successful
!    write(*,*) proc_id,"geom_len    ",geom_len
!    write(*,*) proc_id,"geom_charge ",geom_charge
!    write(*,*) proc_id,"geom_multip ",geom_multip
!    write(*,*) proc_id,"element"
!    write(*,*) element
!    if (flag_init_elabel) then
!      write(*,*) proc_id,"elabel"
!      write(*,*) elabel
!    end if

  end if

end subroutine mmpi_init_images
#endif

!====================================================================

end module geometry

