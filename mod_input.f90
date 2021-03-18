module input

  use utility
  use rotation
  use geometry
  use pes_data
  use pes_input_template
  use pes
  use elastic
  use climbing
  use idpp
  use optimization

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: read_input

contains

!====================================================================
! Public
!====================================================================

subroutine read_input(fname_in)

  !--------------------------------------------------------
  ! Open fname_in input file, read all the informations
  ! and when it's done close the stream and set variables
  !--------------------------------------------------------

  character(*),                intent(IN) :: fname_in

  integer, parameter                      :: fnumb_in = 100
  integer                                 :: ri_start_atoms
  integer                                 :: ri_end_atoms
  character(200)                          :: cmd_str
  character(50)                           :: keyword
  character(50)                           :: arg
  character(50)                           :: ri_pes_program
  character(3), allocatable, dimension(:) :: ri_start_elem
  character(3), allocatable, dimension(:) :: ri_end_elem
  character(3), allocatable, dimension(:) :: ri_start_elab
  character(3), allocatable, dimension(:) :: ri_end_elab
  logical                                 :: got_error
  integer                                 :: err_n
  character(120)                          :: err_msg
  ! got flags check that there are no duplicate keywords in the input file
  ! some of them are also used to perform input consistency checks
  logical                                 :: got_start
  logical                                 :: got_end
  logical                                 :: got_charge
  logical                                 :: got_multip
  logical                                 :: got_images
  logical                                 :: got_start_energy
  logical                                 :: got_end_energy
  logical                                 :: got_new_pes_program
  logical                                 :: got_pes_program
  logical                                 :: got_pes_exec
  logical                                 :: got_geometries_file
  logical                                 :: got_start_elabel
  logical                                 :: got_end_elabel
  logical                                 :: got_idpp
  logical                                 :: got_idpp_conv
  logical                                 :: got_scfcycle
  logical                                 :: got_scfconv
  logical                                 :: got_pes_proc
  logical                                 :: got_pes_mem
  logical                                 :: got_autorotation
  logical                                 :: got_aux_input_files
  logical                                 :: got_aux_output_files
  logical                                 :: got_climbing
  logical                                 :: got_climbing_quick_start
  logical                                 :: got_descending
  logical                                 :: got_descending_quick_start
  logical                                 :: got_only_interpolation
  logical                                 :: got_opt_algo
  logical                                 :: got_opt_cycle
  logical                                 :: got_opt_conv
  logical                                 :: got_springmode
  logical                                 :: got_scfvshift
  logical                                 :: got_intgrid
  logical                                 :: got_additional_cmd

  ! open file ---------------------------------------------
  open(unit=fnumb_in,file=fname_in,status='OLD',action='READ',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("read_input: "//trim(err_msg))
  end if

  ! set variables -----------------------------------------
  got_error                  = .false.
  got_start                  = .false.
  got_end                    = .false.
  got_charge                 = .false.
  got_multip                 = .false.
  got_images                 = .false.
  got_start_energy           = .false.
  got_end_energy             = .false.
  got_new_pes_program        = .false.
  got_pes_program            = .false.
  got_pes_exec               = .false.
  got_geometries_file        = .false.
  got_start_elabel           = .false.
  got_end_elabel             = .false.
  got_idpp                   = .false.
  got_idpp_conv              = .false.
  got_scfcycle               = .false.
  got_scfconv                = .false.
  got_pes_proc               = .false.
  got_pes_mem                = .false.
  got_autorotation           = .false.
  got_aux_input_files        = .false.
  got_aux_output_files       = .false.
  got_climbing               = .false.
  got_climbing_quick_start   = .false.
  got_descending             = .false.
  got_descending_quick_start = .false.
  got_only_interpolation     = .false.
  got_opt_algo               = .false.
  got_opt_cycle              = .false.
  got_opt_conv               = .false.
  got_springmode             = .false.
  got_scfvshift              = .false.
  got_intgrid                = .false.
  got_additional_cmd         = .false.

  ri_start_atoms             = -1
  ri_end_atoms               = -1

  ! read the input file -----------------------------------
  do
    ! Get a string ----------------------------------------
    read(fnumb_in,'(A200)',iostat=err_n) cmd_str
    
    ! Check end of file -----------------------------------
    if (err_n/=0) then
      exit
    end if

    ! Parse command ---------------------------------------
    ! Mutual exclusivity is checked line-by-line

    ! Skip empty strings
    if (len_trim(cmd_str)==0) then
      cycle
    end if

    ! Get the keyword
    call get_field(cmd_str,keyword,1,err_n,err_msg)
    if (err_n/=0) then
      call error("read_input: "//trim(err_msg))
    end if
    
    ! Skip comment strings
    if (keyword(1:1)=="!") then
      cycle
    else if (keyword(1:1)/="#") then
      call error("read_input: in input file, string """&
        &//trim(cmd_str)//""" is not a valid keyword")
    end if

    ! Parse the keyword
    select case (keyword)
    case ("#START")
      if (got_start) then
        call error("read_input: #START specified more than once")
      end if

      if (got_geometries_file) then
        call error("read_input: #GEOMETRIESFILE and #START are mutually exclusive")
      end if

      call read_geometries_from_input(fnumb_in,keyword,ri_start_atoms,ri_start_elem,ri_start_elab,got_start_elabel)
      got_start = .true.

    case ("#END")
      if (got_end) then
        call error("read_input: #END specified more than once")
      end if

      if (got_geometries_file) then
        call error("read_input: #GEOMETRIESFILE and #END are mutually exclusive")
      end if

      call read_geometries_from_input(fnumb_in,keyword,ri_end_atoms,ri_end_elem,ri_end_elab,got_end_elabel)
      got_end = .true.

    case ("#GEOMETRIESFILE")
      if (got_geometries_file) then
        call error("read_input: #GEOMETRIESFILE specified more than once")
      end if

      if (got_start) then
        call error("read_input: #START and #GEOMETRIESFILE are mutually exclusive")
      end if

      if (got_end) then
        call error("read_input: #END and #GEOMETRIESFILE are mutually exclusive")
      end if

      if (got_images) then
        call error("read_input: #IMAGES and #GEOMETRIESFILE are mutually exclusive")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call read_geometries_from_file(trim(arg))
      got_geometries_file = .true.

    case ("#CHARGE")
      if (got_charge) then
        call error("read_input: #CHARGE specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_geom_charge(trim(arg))
      got_charge = .true.

    case ("#MULTIP")
      if (got_multip) then
        call error("read_input: #MULTIP specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_geom_multip(trim(arg))
      got_multip = .true.

    case ("#IMAGES")
      if (got_images) then
        call error("read_input: #IMAGES specified more than once")
      end if

      if (got_geometries_file) then
        call error("read_input: #GEOMETRIESFILE and #IMAGES are mutually exclusive")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_image_n(trim(arg))
      got_images = .true.

    case ("#STARTENERGY")
      if (got_start_energy) then
        call error("read_input: #STARTENERGY specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_start_energy(trim(arg))
      got_start_energy = .true.

    case ("#ENDENERGY")
      if (got_end_energy) then
        call error("read_input: #ENDENERGY specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_end_energy(trim(arg))
      got_end_energy = .true.

    case ("#IDPP")
      if (got_idpp) then
        call error("read_input: #IDPP specified more than once")
      end if

      call set_idpp(.true.)
      got_idpp = .true.

    case ("#NEWPESPROGRAM")
      if (got_new_pes_program) then
        call error("read_input: #NEWPESPROGRAM specified more than once")
      end if

      call set_new_pes_program(.true.)
      got_new_pes_program = .true.

    case ("#PESPROGRAM")
      if (got_pes_program) then
        call error("read_input: #PESPROGRAM specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      ri_pes_program  = arg
      got_pes_program = .true.

    case ("#PESEXEC")
      if (got_pes_exec) then
        call error("read_input: #PESEXEC specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_pes_exec(arg)
      got_pes_exec = .true.

    case ("#PESINPUTTEMPLATE")
      ! no duplicate checking because this keyword can be specified more than once

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if

      call read_pes_it(arg,fnumb_in,"#ENDPESINPUTTEMPLATE")

    case ("#PESPROC")
      if (got_pes_proc) then
        call error("read_input: #PESPROC specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_pes_proc(arg)
      got_pes_proc = .true.

    case ("#PESMEM")
      if (got_pes_mem) then
        call error("read_input: #PESMEM specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_pes_mem(arg)
      got_pes_mem = .true.

    case ("#AUTOROTATION")
      if (got_autorotation) then
        call error("read_input: #AUTOROTATION specified more than once")
      end if

      call set_rotation(.true.)
      got_autorotation = .true.

    case ("#CLIMBING")
      if (got_climbing) then
        call error("read_input: #CLIMBING specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_climbing_image(arg)
      got_climbing = .true.

    case ("#CLIMBINGQUICKSTART")
      if (got_climbing_quick_start) then
        call error("read_input: #CLIMBINGQUICKSTART specified more than once")
      end if

      call set_climbing_quick_start(.true.)
      got_climbing_quick_start = .true.

    case ("#DESCENDING")
      if (got_descending) then
        call error("read_input: #DESCENDING specified more than once")
      end if

      call set_descending_image(.true.)
      got_descending = .true.

    case ("#DESCENDINGQUICKSTART")
      if (got_descending_quick_start) then
        call error("read_input: #DESCENDINGQUICKSTART specified more than once")
      end if

      call set_descending_quick_start(.true.)
      got_descending_quick_start = .true.

    case ("#ONLYINTERPOLATION")
      if (got_only_interpolation) then
        call error("read_input: #ONLYINTERPOLATION specified more than once")
      end if

      call set_only_interpolation(.true.)
      got_only_interpolation = .true.

    case ("#SPRINGMODE")
      if (got_springmode) then
        call error("read_input: #SPRINGMODE specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if

      call set_spring_mode(arg)
      got_springmode=.true.

    case ("#SCFCYCLE")
      if (got_scfcycle) then
        call error("read_input: #SCFCYCLE specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_pesd_scfcycle(arg)
      got_scfcycle = .true.

    case ("#SCFCONV")
      if (got_scfconv) then
        call error("read_input: #SCFCONV specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_pesd_scfconv(arg)
      got_scfconv = .true.

    case ("#SCFVSHIFT")
      if (got_scfvshift) then
        call error("read_input: #SCFVSHIFT specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_pesd_scfvshift(arg)
      got_scfvshift = .true.

    case ("#INTGRID")
      if (got_intgrid) then
        call error("read_input: #INTGRID specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_pesd_intgrid(arg)
      got_intgrid = .true.

    case ("#ADDITIONALCMD")
      if (got_additional_cmd) then
        call error("read_input: #ADDITIONALCMD specified more than once")
      end if

      read(fnumb_in,'(A200)',iostat=err_n) cmd_str
      if (err_n/=0) then
        call error("read_input: cannot read additional command")
      end if

      call set_pesd_additional_cmd(trim(cmd_str))
      got_additional_cmd = .true.

    case ("#OPTALGORITHM")
      if (got_opt_algo) then
        call error("read_input: #OPTALGORITHM specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if
 
      call set_optmz_algo(arg)
      got_opt_algo = .true.

    case ("#OPTCYCLE")
      if (got_opt_cycle) then
        call error("read_input: #OPTCYCLE specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if

      call set_optmz_nsteps(arg)
      got_opt_cycle = .true.
 
    case ("#OPTCONV")
      if (got_opt_conv) then
        call error("read_input: #OPTCONV specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if

      call set_optmz_tol(arg)
      got_opt_conv = .true.
 
    case ("#IDPPCONV")
      if (got_idpp_conv) then
        call error("read_input: #IDPPCONV specified more than once")
      end if

      call get_field(cmd_str,arg,2,err_n,err_msg)
      if (err_n/=0) then
        call error("read_input: "//trim(err_msg))
      end if

      call set_idpp_tol(arg)
      got_idpp_conv = .true.

    case ("#AUXINPUTFILES")
      if (got_aux_input_files) then
        call error("read_input: #AUXINPUTFILES specified more than once")
      end if

      call set_pesd_auxiliary_files(PESD_AUX_INPUT_FILES,cmd_str)
      got_aux_input_files = .true.

    case ("#AUXOUTPUTFILES")
      if (got_aux_output_files) then
        call error("read_input: #AUXOUTPUTFILES specified more than once")
      end if

      call set_pesd_auxiliary_files(PESD_AUX_OUTPUT_FILES,cmd_str)
      got_aux_output_files = .true.

    ! keywords to be ignored
    case ("#ENDPESINPUTTEMPLATE")

    case default
      call error("read_input: unknown keyword """//trim(keyword)//"""")

    end select
  end do

  close(unit=fnumb_in,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("read_input: "//trim(err_msg))
  end if

  ! set pes program ---------------------------------------
  call set_pes_program(ri_pes_program)

  ! input checks ------------------------------------------
  call check_neb_kw(got_new_pes_program,got_pes_program,&
    &got_pes_exec,got_geometries_file,got_start,got_start_energy,got_end,&
    &got_end_energy,got_images,got_scfcycle,got_scfconv,got_idpp)

  call consistency_check(got_geometries_file,ri_start_atoms,&
    &ri_end_atoms,ri_start_elem,ri_end_elem,got_start_elabel,got_end_elabel,&
    &ri_start_elab,ri_end_elab)

  ! set elements and labels -------------------------------
  if (.not.got_geometries_file) then
    call set_geom_len(ri_start_atoms)
    if (allocated(ri_start_elem)) then
      call allocate_element()
      call update_element(ri_start_elem)
    end if

    if (got_start_elabel) then
      call allocate_elabel()
      call update_elabel(ri_start_elab)
    end if
  end if

  ! check the mandatory keywords for the external program -
  call check_gaussian_mandatory_kw()
  call check_siesta_mandatory_kw()

  ! deallocation section ----------------------------------
  if (allocated(ri_start_elem)) then
    deallocate(ri_start_elem,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("read_input: "//trim(err_msg))
    end if
  end if

  if (allocated(ri_end_elem)) then
    deallocate(ri_end_elem,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("read_input: "//trim(err_msg))
    end if
  end if

  if (allocated(ri_start_elab)) then
    deallocate(ri_start_elab,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("read_input: "//trim(err_msg))
    end if
  end if

  if (allocated(ri_end_elab)) then
    deallocate(ri_end_elab,stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("read_input: "//trim(err_msg))
    end if
  end if

end subroutine read_input

!====================================================================
! Private
!====================================================================

subroutine read_geometries_from_input(fnumb_in,point,atoms,elem,elab,found_elab)

  integer,                                 intent(IN)    :: fnumb_in
  character(*),                            intent(IN)    :: point
  integer,                                 intent(OUT)   :: atoms
  character(*), allocatable, dimension(:), intent(INOUT) :: elem
  character(*), allocatable, dimension(:), intent(INOUT) :: elab
  logical,                                 intent(OUT)   :: found_elab

  character(*), parameter                                :: my_name = "read_geometries_from_input"
  integer                                                :: g_len
  character(120)                                         :: str
  integer                                                :: err_n
  character(120)                                         :: err_msg

  ! read the number of atoms ------------------------------
  read(fnumb_in,*,iostat=err_n) str
  if (err_n/=0) then
    call error(my_name//": atom number not specified")
  end if

  if (isinteger(trim(adjustl(str)))) then
    read(str,*) atoms
  else
    call error(my_name//": got """//trim(str)//""", expected integer")
  end if

  if (atoms<=0) then
    call error(my_name//": non-positive atom number in input file")
  end if

  ! allocate start_geom or end_geom -----------------------
  g_len = 3*atoms

  call allocate_geom(point,g_len)

  ! allocation section ------------------------------------
  if (allocated(elem)) then
    call error(my_name//": design error, passed an already allocated elem")
  else
    allocate(elem(atoms),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error(my_name//": "//trim(err_msg))
    end if
  end if

  if (allocated(elab)) then
    call error(my_name//": design error, passed an already allocated elab")
  else
    allocate(elab(atoms),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error(my_name//": "//trim(err_msg))
    end if
  end if

  ! read geometry -----------------------------------------
  backspace(unit=fnumb_in,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//err_msg)
  end if

  select case (point)
  case ("#START")
    call read_xyz(fnumb_in,elem,start_geom,elab,found_elab,.false.)
  case ("#END")
    call read_xyz(fnumb_in,elem,end_geom,elab,found_elab,.false.)
  case default
    call error(my_name//": unknown argument "//trim(point))
  end select

end subroutine read_geometries_from_input

!====================================================================
! Subroutines that check
!====================================================================

subroutine check_neb_kw(got_new_pes_program,got_pes_program,&
    &got_pes_exec,got_geometries_file,got_start,got_start_energy,got_end,&
    &got_end_energy,got_images,got_scfcycle,got_scfconv,got_idpp)

  logical,     intent(IN) :: got_new_pes_program
  logical,     intent(IN) :: got_pes_program
  logical,     intent(IN) :: got_pes_exec
  logical,     intent(IN) :: got_geometries_file
  logical,     intent(IN) :: got_start
  logical,     intent(IN) :: got_start_energy
  logical,     intent(IN) :: got_end
  logical,     intent(IN) :: got_end_energy
  logical,     intent(IN) :: got_images
  logical,     intent(IN) :: got_scfcycle
  logical,     intent(IN) :: got_scfconv
  logical,     intent(IN) :: got_idpp

  character(*), parameter :: my_name   = "check_neb_kw"
  logical                 :: got_error

  got_error = .false.

  if (.not.got_pes_program) then
    write(FILEOUT,* )"WAR "//my_name//": #PESPROGRAM was not specified"
    got_error = .true.
  end if

  if (.not.got_pes_exec) then
    write(FILEOUT,*) "WAR "//my_name//": #PESEXEC was not specified"
    got_error = .true.
  end if

  if (.not.got_geometries_file) then
    if (.not.got_start) then
      write(FILEOUT,*) "WAR "//my_name//": #START was not specified"
      got_error = .true.
    end if

    if (.not.got_end) then
      write(FILEOUT,*) "WAR "//my_name//": #END was not specified"
      got_error = .true.
    end if

    if (.not.got_images) then
      write(FILEOUT,*) "WAR "//my_name//": #IMAGES was not specified"
      got_error = .true.
    end if
  end if

  if (got_geometries_file) then
    if (got_idpp) then
      write(FILEOUT,*) "WAR "//my_name//": #IDPP can not be used with #GEOMETRIESFILE"
      got_error = .true.
    end if
  end if

  if (.not.got_start_energy) then
    write(FILEOUT,*) "WAR "//my_name//": #STARTENERGY was not specified"
    got_error = .true.
  end if

  if (.not.got_end_energy) then
    write(FILEOUT,*) "WAR "//my_name//": #ENDENERGY was not specified"
    got_error = .true.
  end if

  if (got_new_pes_program) then
    if (.not.got_scfcycle) then
      write(FILEOUT,*) "WAR "//my_name//": #NEWPESPROGRAM used but #SCFCYCLE was not specified"
      got_error = .true.
    end if

    if (.not.got_scfconv) then
      write(FILEOUT,*) "WAR "//my_name//": #NEWPESPROGRAM used but #SCFCONV was not specified"
      got_error = .true.
    end if
  end if

  if (got_error) then
    call error(my_name//": found errors in the input file, please correct and try again")
  end if

end subroutine check_neb_kw

!====================================================================

subroutine check_gaussian_mandatory_kw()

  character(*), parameter :: my_name   = "check_gaussian_mandatory_kw"
  logical                 :: got_error
  integer                 :: indx
  integer                 :: i
  character(8)            :: i_str

  got_error = .false.

  if (pes_program/="gaussian") then
    return
  end if

  ! check that both #PESINPUTTEMPLATE 1 and 2 blocks were specified
  do i=1, 2
    indx = get_pes_it_n(i)

    if (indx==0) then
      write(i_str,'(I8)') i
      i_str=adjustl(i_str)
      write(FILEOUT,*) "WAR "//my_name//": ""#PESINPUTTEMPLATE "//trim(i_str)//""" block not specified"
      got_error = .true.
    end if
  end do

  if (got_error) then
    call error(my_name//": some necessary keywords were not specified")
  end if

end subroutine check_gaussian_mandatory_kw

!====================================================================

subroutine check_siesta_mandatory_kw()

  character(*), parameter :: my_name   = "check_siesta_mandatory_kw"
  logical                 :: got_error
  integer                 :: indx
  integer                 :: i
  character(8)            :: i_str

  got_error = .false.

  if (pes_program/="siesta") then
    return
  end if

  ! check that both #PESINPUTTEMPLATE 1 and 2 blocks were specified
  do i=1, 2
    indx = get_pes_it_n(i)

    if (indx==0) then
      write(i_str,'(I8)') i
      i_str=adjustl(i_str)
      write(FILEOUT,*) "WAR "//my_name//": ""#PESINPUTTEMPLATE "//trim(i_str)//""" block not specified"
      got_error = .true.
    end if
  end do

  ! check that element labels were specified
  if (flag_elabel.eqv..false.) then
    write(FILEOUT,*) "WAR "//my_name//": element labels not specified in the geometry blocks"
    got_error = .true.
  end if

  if (got_error) then
    call error(my_name//": some necessary keywords were not specified")
  end if

end subroutine check_siesta_mandatory_kw

!====================================================================

subroutine consistency_check(got_geometries_file,ri_start_atoms,&
    &ri_end_atoms,ri_start_elem,ri_end_elem,got_start_elabel,got_end_elabel,&
    &ri_start_elab,ri_end_elab)

  logical,                                 intent(IN) :: got_geometries_file
  integer,                                 intent(IN) :: ri_start_atoms
  integer,                                 intent(IN) :: ri_end_atoms
  character(*), allocatable, dimension(:), intent(IN) :: ri_start_elem
  character(*), allocatable, dimension(:), intent(IN) :: ri_end_elem
  logical,                                 intent(IN) :: got_start_elabel
  logical,                                 intent(IN) :: got_end_elabel
  character(*), allocatable, dimension(:), intent(IN) :: ri_start_elab
  character(*), allocatable, dimension(:), intent(IN) :: ri_end_elab

  integer                                             :: i
  character(8)                                        :: i_str
  logical                                             :: flag_alloc_s
  logical                                             :: flag_alloc_e
  logical                                             :: flag_not_consistent

  flag_not_consistent = .false.

  ! those checks must be done only if ---------------------
  ! geometries file was not specified
  if (got_geometries_file.eqv..false.) then
    ! check number of atoms -------------------------------
    if (ri_start_atoms<=0) then
      write(FILEOUT,*) "WAR consistency_check: non-positive start number of atoms"
      flag_not_consistent = .true.
    end if

    if (ri_end_atoms<=0) then
      write(FILEOUT,*) "WAR consistency_check: non-positive end number of atoms"
      flag_not_consistent = .true.
    end if

    if (ri_start_atoms/=ri_end_atoms) then
      write(FILEOUT,*) "WAR consistency_check: "//&
        &"start and end geometries contain a different number of atoms"
      flag_not_consistent = .true.
    end if

    ! check elements array and order ------------------------
    flag_alloc_s = allocated(ri_start_elem)
    flag_alloc_e = allocated(ri_end_elem)

    if ((flag_alloc_s).and.(.not.flag_alloc_e)) then
      write(FILEOUT,*) "WAR consistency_check: "//&
        &"only start elements array is present"
      flag_not_consistent = .true.
    else if ((.not.flag_alloc_s).and.(flag_alloc_e)) then
      write(FILEOUT,*) "WAR consistency_check: "//&
        &"only end elements array is present"
      flag_not_consistent = .true.
    end if

    if (flag_alloc_s.and.flag_alloc_e) then
      if (size(ri_start_elem,1)/=ri_start_atoms) then
        write(FILEOUT,*) "WAR consistency_check: "//&
          &"wrong start elements array size"
        flag_not_consistent = .true.
      end if
      if (size(ri_end_elem,1)/=ri_end_atoms) then
        write(FILEOUT,*) "WAR consistency_check: "//&
          &"wrong end elements array size"
        flag_not_consistent = .true.
      end if

      if (flag_not_consistent.eqv..false.) then
        do i=1, size(ri_start_elem,1)
          if (ri_start_elem(i)/=ri_end_elem(i)) then
            write(i_str,'(I8)') i
            i_str = adjustl(i_str)
            write(FILEOUT,*) "WAR consistency_check: "//&
              &"elements in line "//trim(i_str)//" are different"
            flag_not_consistent = .true.
          end if
        end do
      end if
    end if

    ! check element label -----------------------------------
    if (got_start_elabel.neqv.got_end_elabel) then
      write(FILEOUT,*) "WAR consistency_check: element labels must be specified "//&
        &"in all the geometry blocks or none"
      flag_not_consistent = .true.
    else
      if (got_start_elabel.eqv..true.) then
        if (size(ri_start_elab,1)/=ri_start_atoms) then
          write(FILEOUT,*) "WAR consistency_check: wrong start elabel array size"
          flag_not_consistent = .true.
        end if

        if (size(ri_end_elab,1)/=ri_end_atoms) then
          write(FILEOUT,*) "WAR consistency_check: wrong end elabel array size"
          flag_not_consistent = .true.
        end if

        if (flag_not_consistent.eqv..false.) then
          do i=1, size(ri_start_elab,1)
            if (ri_start_elab(i)/=ri_end_elab(i)) then
              write(i_str,'(I8)') i
              i_str = adjustl(i_str)
              write(FILEOUT,*) "WAR consistency_check: "//&
                &"labels in line "//trim(i_str)//" are different"
              flag_not_consistent = .true.
            end if
          end do
        end if
      end if
    end if
  end if

  ! check climbing image number ---------------------------
  if (climbing_image_n>image_n) then
    write(FILEOUT,*) "WAR consistency_check: "//&
      &"specified more climbing images than normal ones"
    flag_not_consistent = .true.
  end if

  ! check climbing image quick start ----------------------
  if (flag_climbing_quick_start) then
    if (flag_climbing_image.eqv..false.) then
      write(FILEOUT,*) "WAR consistency_check: "//&
        &"#CLIMBINGQUICKSTART requires #CLIMBING"
      flag_not_consistent = .true.
    end if
  end if

  ! result ------------------------------------------------
  if (flag_not_consistent) then
    call error("consistency_check: inconsistent input")
  end if

end subroutine consistency_check

!====================================================================
! Subroutines that read geometries file
!====================================================================

subroutine read_geometries_from_file(gf_fname)

  !--------------------------------------------------------
  ! Reads geometries from gf_fname and
  ! executes consistency checks on the input file.
  ! Finally, sets flag_geometries_file on true.
  !--------------------------------------------------------

  character(*),                intent(IN) :: gf_fname

  character(*), parameter                 :: my_name  = "read_geometries_from_file"
  integer, parameter                      :: gf_fnumb = 101
  integer                                 :: i
  integer                                 :: j
  integer                                 :: ngeom
  integer                                 :: natom
  character(120)                          :: str
  character(3), allocatable, dimension(:) :: elem_arr
  character(3), allocatable, dimension(:) :: elem_dfl
  character(3), allocatable, dimension(:) :: elab_arr
  character(3), allocatable, dimension(:) :: elab_dfl
  real(DBL),    allocatable, dimension(:) :: geom_arr
  logical                                 :: found_elab
  logical                                 :: found_elab_dfl
  integer                                 :: err_n
  character(120)                          :: err_msg

  ! open unit ---------------------------------------------
  open(unit=gf_fnumb,file=gf_fname,status='OLD',action='READ',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! get the info ------------------------------------------
  call get_geometries_info(gf_fnumb,ngeom,natom)

  if (ngeom<=2) then
    call error(my_name//": "//trim(gf_fname)//" must contain 3 or more geometries")
  end if

  if (natom<=0) then
    call error(my_name//": bad format in geometries file")
  end if

  write(str,'(I8)') ngeom-2
  str = adjustl(str)
  call set_geom_len(natom)
  call set_image_n(trim(str))
  call allocate_image_geom()
  call allocate_element()

  ! allocation section ------------------------------------
  allocate(elem_arr(geom_len/3),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(elem_dfl(geom_len/3),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(elab_arr(geom_len/3),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(elab_dfl(geom_len/3),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(geom_arr(geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! read all the geometries -------------------------------
  rewind(unit=gf_fnumb,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  do i=0, image_n+1
    call read_xyz(gf_fnumb,elem_arr,geom_arr,elab_arr,found_elab,.true.)
    
    if (i==0) then
      elem_dfl = elem_arr
      call update_element(elem_dfl)

      found_elab_dfl = found_elab
      if (found_elab_dfl) then
        call allocate_elabel()
        elab_dfl = elab_arr
        call update_elabel(elab_dfl)
      end if
    else
      if (found_elab.neqv.found_elab_dfl) then
        call error(my_name//": labels must be specified for all geometries or none")
      end if

      do j=1, size(elem_arr,1)
        if (elem_dfl(j)/=elem_arr(j)) then
          call error(my_name//": elements in input geometries are inconsistent")
        end if
      end do

      if (found_elab_dfl) then
        do j=1, size(elab_arr,1)
          if (elab_dfl(j)/=elab_arr(j)) then
            call error(my_name//": labels in input geometries are inconsistent")
          end if
        end do
      end if
    end if

    call update_image_geom(i,geom_arr)
  end do

  ! deallocate stuffs -------------------------------------
  deallocate(elem_arr,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(elem_dfl,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(elab_arr,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(elab_dfl,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(geom_arr,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! close unit --------------------------------------------
  close(unit=gf_fnumb,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  call set_geometries_file(.true.)

end subroutine read_geometries_from_file

!====================================================================

subroutine get_geometries_info(fnumb,ngeom,natom)

  !--------------------------------------------------------
  ! Returns some info on file unit fnumb, like the number
  ! of geometries (ngeom) and the number of atoms
  ! for each geometry (natom).
  ! It also executes some consistency checks on fnumb format.
  ! Being a private function, it assumes that fnumb is opened.
  !--------------------------------------------------------

  integer, intent(IN)  :: fnumb
  integer, intent(OUT) :: ngeom
  integer, intent(OUT) :: natom

  character(120)       :: str
  integer              :: i
  integer              :: tmp_atom
  integer              :: err_n

  ngeom = 0
  do
    read(fnumb,'(A120)',iostat=err_n) str
    if (err_n/=0) then
      exit
    end if

    str = adjustl(str)
    if (.not.isinteger(trim(adjustl(str)))) then
      call error("get_geometries_info: bad format in geometries file")
    end if

    if (ngeom==0) then
      read(str,*) natom
    else
      read(str,*) tmp_atom
      if (tmp_atom/=natom) then
        call error("get_geometries_info: geometries must have same number of atoms")
      end if
    end if

    do i=1,natom+1
      read(fnumb,'(A120)',iostat=err_n) str
      if (err_n/=0) then
        call error("get_geometries_info: bad format in geometries file")
      end if
    end do

    ngeom = ngeom+1
  end do

end subroutine get_geometries_info

!====================================================================

subroutine read_xyz(fnumb,elem_arr,geom_arr,elab_arr,found_elab,comment_line)

  !--------------------------------------------------------
  ! Reads an xyz format geometry from fnumb.
  ! Writes elements encountered in elem_arr,
  ! coordinates in geom_arr, and if they are present,
  ! element lables in elab_arr.
  ! Can be used with or without a comment line after
  ! the first one that contains the number of atoms.
  !--------------------------------------------------------

  integer,                    intent(IN)  :: fnumb
  character(*), dimension(:), intent(OUT) :: elem_arr
  real(DBL),    dimension(:), intent(OUT) :: geom_arr
  character(*), dimension(:), intent(OUT) :: elab_arr
  logical,                    intent(OUT) :: found_elab
  logical,                    intent(IN)  :: comment_line

  character(*), parameter                 :: my_name = "read_xyz"
  logical                                 :: is_open
  logical                                 :: elabels_are_specified
  integer                                 :: i
  integer                                 :: j
  integer                                 :: i_start
  integer                                 :: atoms
  character(8)                            :: i_str
  character(120)                          :: str
  character(40)                           :: field
  integer                                 :: err_n
  character(120)                          :: err_msg

  ! preliminary checks ------------------------------------
  inquire(unit=fnumb,opened=is_open)
  if (.not.is_open) then
    call error(my_name//": input stream not opened")
  end if

  ! preliminary sets --------------------------------------
  elabels_are_specified = .true.

  ! read the number of atoms ------------------------------
  read(fnumb,'(A120)',iostat=err_n) str
  if (err_n/=0) then
    call error(my_name//": line with number of atoms not found")
  end if
  str = adjustl(str)

  if (.not.isinteger(trim(adjustl(str)))) then
    call error(my_name//": line with number of atoms is not an integer")
  end if
  read(str,*) atoms

  ! check size of argument arrays -------------------------
  if (size(elem_arr,1)/=atoms) then
    call error(my_name//": wrong size of elem_arr argument")
  end if

  if (size(geom_arr,1)/=(atoms*3)) then
    call error(my_name//": wrong size of geom_arr argument")
  end if

  if (size(elab_arr,1)/=atoms) then
    call error(my_name//": wrong size of elab_arr argument")
  end if

  ! set loop starting point -------------------------------
  if (comment_line.eqv..true.) then
    i_start = 0
  else
    i_start = 1
  end if

  ! get the geometry --------------------------------------
  do i=i_start, atoms
    read(fnumb,'(A120)',iostat=err_n) str
    if (err_n/=0) then
      write(i_str,'(I8)') i
      i_str = adjustl(i_str)
      call error(my_name//": missing line "//trim(i_str)//" while reading geometry")
    end if
    str = adjustl(str)

    if (i==0) then
      cycle
    end if

    call get_field(str,field,1,err_n,err_msg)
    if (err_n/=0) then
      call error(my_name//": "//trim(err_msg))
    end if

    elem_arr(i) = field

    do j=2, 4
      call get_field(str,field,j,err_n,err_msg)
      if (err_n/=0) then
        call error(my_name//": "//trim(err_msg))
      end if

      if (.not.isreal(trim(adjustl(field)))) then
        call error(my_name//": wrong formatting for coordinates, real numbers were expected")
      end if
      read(field,*) geom_arr(3*(i-1)+j-1)
    end do

    call get_field(str,field,5,err_n,err_msg)
    if (err_n/=0) then ! field was not found
      if (i==1) then
        elabels_are_specified = .false.
      else
        if (elabels_are_specified) then
          call error(my_name//": labels must be specified for all the elements or none; "//trim(err_msg))
        end if
      end if
    else               ! field was found
      if (elabels_are_specified) then
        elab_arr(i) = field
      else
        call error(my_name//": labels must be specified for all the elements or none")
      end if
    end if
  end do

  found_elab = elabels_are_specified

end subroutine read_xyz

!====================================================================

end module input

