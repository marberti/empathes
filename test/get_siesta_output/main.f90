program main

  use iso_fortran_env
  use utility

  implicit none

  integer, parameter :: file_count = 3
  character(*), dimension(file_count), parameter :: force_file  = &
    &["test1.forces","test2.forces","test3.forces"]
  character(*), dimension(file_count), parameter :: output_file = &
    &["test1.output","test2.output","test3.output"]

  integer :: i
  logical :: flag_conv
  logical :: same_forces
  logical :: verbose
  integer :: geom_len
  real(DBL), dimension(:),   allocatable :: pes_energy
  real(DBL), dimension(:,:), allocatable :: pes_forces
  integer :: err_n
  character(120) :: err_msg

!  verbose = .true.
  verbose = .false.

  ! open FILEOUT ------------------------------------------
  open(unit=FILEOUT,file="output.ignore",action='WRITE',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("main: "//trim(err_msg))
  end if

  ! tests -------------------------------------------------
  do i=1, file_count
    write(*,'("Test ",I2,"...")') i

    call set_geom_len(force_file(i))
    call my_alloc()
    call read_forces(force_file(i))
    call get_siesta_output(2,".",100,output_file(i),flag_conv)
    call check_forces(same_forces,verbose)

    if (same_forces) then
      write(*,'("Test ",I2,": passed")') i
    else
      write(*,'("Test ",I2,": not passed")') i
    end if

    call my_dealloc()
  end do

  ! close FILEOUT -----------------------------------------
  close(unit=FILEOUT,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("main: "//trim(err_msg))
  end if

contains

!====================================================================

subroutine set_geom_len(fname)

  character(*), intent(IN) :: fname
  integer, parameter :: unit_n = 50
  integer :: i
  integer :: atoms
  integer :: err_n
  character(120) :: err_msg

  ! open unit ---------------------------------------------
  open(unit=unit_n,file=fname,status='OLD',action='READ',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("read_forces: "//trim(err_msg))
  end if

  ! read number of atoms ----------------------------------
  read(unit_n,*) atoms
  geom_len = atoms*3

  ! close unit --------------------------------------------
  close(unit=unit_n,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("read_forces: "//trim(err_msg))
  end if

end subroutine set_geom_len

!====================================================================

subroutine read_forces(fname)

  character(*), intent(IN) :: fname
  integer, parameter :: unit_n = 50
  integer :: i
  integer :: atoms
  integer :: err_n
  character(120) :: err_msg

  ! open unit ---------------------------------------------
  open(unit=unit_n,file=fname,status='OLD',action='READ',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("read_forces: "//trim(err_msg))
  end if

  ! read number of atoms ----------------------------------
  read(unit_n,*) atoms

  ! read forces -------------------------------------------
  do i = 1, atoms
    read(unit_n,*) pes_forces(1,3*(i-1)+1),&
      &pes_forces(1,3*(i-1)+2), pes_forces(1,3*(i-1)+3)
  end do

  ! close unit --------------------------------------------
  close(unit=unit_n,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("read_forces: "//trim(err_msg))
  end if

end subroutine read_forces

!====================================================================

integer function get_scfcycle()

  get_scfcycle = 10000

end function get_scfcycle

!====================================================================

subroutine my_alloc()

  allocate(pes_energy(2))
  allocate(pes_forces(2,geom_len))

  pes_energy = 0.0_DBL
  pes_forces = 0.0_DBL

end subroutine my_alloc

!====================================================================

subroutine my_dealloc()

  deallocate(pes_energy)
  deallocate(pes_forces)

end subroutine my_dealloc

!====================================================================
    
subroutine check_forces(same_forces,verbose)

  logical, intent(OUT) :: same_forces
  logical, intent(IN)  :: verbose
  integer :: i

  same_forces = .true.
  do i=1, geom_len
    if (verbose) write(*,'(I4,":",2(" ",F9.6))',advance="no")&
      &i, pes_forces(1,i), pes_forces(2,i)

    if (pes_forces(1,i) /= pes_forces(2,i)) then
      if (verbose) write(*,*) " different"
      same_forces = .false.
      if (.not.verbose) exit
    else
      if (verbose) write(*,*) " same"
    end if
  end do

end subroutine check_forces

!====================================================================

subroutine get_siesta_output(i,dirname,fnumb_out,fname_out,flag_conv)

  integer,      intent(IN)  :: i
  character(*), intent(IN)  :: dirname
  integer,      intent(IN)  :: fnumb_out
  character(*), intent(IN)  :: fname_out
  logical,      intent(OUT) :: flag_conv
  integer :: j
  integer :: k
  integer :: k_start
  integer :: default_scfcycle ! max scf cycles for the siesta computation
  integer :: read_scfcycle    ! last scf cycle read from siesta output
  logical :: correct_one
  character(200) :: field
  character(200) :: str
  character(120) :: path
  integer :: err_n
  character(120) :: err_msg

  ! variables initialization ------------------------------
  flag_conv        = .true.
  correct_one      = .false.
  path             = trim(dirname)//"/"//trim(fname_out)
  read_scfcycle    = -1
  default_scfcycle = get_scfcycle()
  default_scfcycle = default_scfcycle - 1 ! If we tell siesta
    ! to do N scf cycles, we will find at most N - 1 scf cycles
    ! in the output file

  ! open unit ---------------------------------------------
  open(unit=fnumb_out,file=path,status='OLD',action='READ',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("get_siesta_output: "//trim(err_msg))
  end if

  ! search for convergence --------------------------------
  ! read until the string "siesta: iscf" or "        iscf" is found
  do
    read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
    if (err_n/=0) then
      write(FILEOUT,*) "WAR get_siesta_output: ",&
        &"string ""siesta: iscf"" or ""        iscf"" not found"
      call error("get_siesta_output: "//trim(err_msg))
    end if

    if ((str(1:12)=="siesta: iscf").or.(str(1:12)=="        iscf")) then
      exit
    end if
  end do

  ! get the last scf iteration number
  do
    read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
    if (err_n/=0) then
      write(FILEOUT,*) "WAR get_siesta_output: ",&
        &"last scf iteration not found"
      call error("get_siesta_output: "//trim(err_msg))
    end if

    if (str=="") then
      if (read_scfcycle==-1) then
        call error("get_siesta_output: scf iterations not found")
      else if (read_scfcycle>=default_scfcycle) then
        flag_conv=.false.
      else
        flag_conv=.true.
      end if
      
      exit
    end if

    call get_field(str,field,2,err_n,err_msg)
    if (err_n/=0) then
      call error("get_siesta_output: "//trim(err_msg))
    end if

    if (isinteger(trim(adjustl(field)))) then
      read(field,*) read_scfcycle
    end if
  end do

  ! if converged, get energy and forces -------------------
  if (flag_conv) then
    ! read until the string "siesta: E_KS(eV)" is found,
    ! then get the total energy value
    do
      read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
      if (err_n/=0) then
        write(FILEOUT,*) "WAR get_siesta_output: ",&
          &"string ""siesta: E_KS(eV)"" not found"
        call error("get_siesta_output: "//trim(err_msg))
      end if

      if (str(1:16)=="siesta: E_KS(eV)") then
        call get_field(str,field,4,err_n,err_msg)
        if (err_n/=0) then
          call error("get_siesta_output: "//trim(err_msg))
        end if

        if (isreal(trim(adjustl(field)))) then
          read(field,*) pes_energy(i)
        else
          call error("get_siesta_output: wrong field """//&
            &trim(field)//""" while reading total energy")
        end if

        exit
      end if
    end do

    ! read until the string "siesta: Atomic forces" is found
    do
      read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
      if (err_n/=0) then
        write(FILEOUT,*) "WAR get_siesta_output: ",&
          &"string ""siesta: Atomic forces"" not found"
        call error("get_siesta_output: "//trim(err_msg))
      end if

      if (str(1:21)=="siesta: Atomic forces") then
        ! check if it's the correct one
        read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
        if (err_n/=0) then
          write(FILEOUT,*) "WAR get_siesta_output: ",&
            &"string ""siesta: Atomic forces"" not found"
          call error("get_siesta_output: "//trim(err_msg))
        end if

        call get_field(str,field,1,err_n,err_msg)
        if (err_n/=0) then
          call error("get_siesta_output: "//trim(err_msg))
        end if

        if (field=="siesta:") then
          k_start = 3
          correct_one = .true.
        else if (isinteger(trim(adjustl(field)))) then
          k_start = 2
          correct_one = .true.
        end if

        if (correct_one) then
          backspace(unit=fnumb_out,iostat=err_n,iomsg=err_msg)
          if (err_n/=0) then
            write(FILEOUT,*) "WAR get_siesta_output: ",&
              &"error while backspacing"
            call error("get_siesta_output: "//trim(err_msg))
          end if

          ! continue to get forces values
          exit
        end if
      end if
    end do

    ! get forces values
    do j=1, geom_len/3
      read(fnumb_out,'(A200)',iostat=err_n,iomsg=err_msg) str
      if (err_n/=0) then
        write(FILEOUT,*) "WAR get_siesta_output: ",&
          &"error while reading atomic forces"
        call error("get_siesta_output: "//trim(err_msg))
      end if

      do k=k_start, k_start+2
        call get_field(str,field,k,err_n,err_msg)
        if (err_n/=0) then
          call error("get_siesta_output: "//trim(err_msg))
        end if

        if (isreal(trim(adjustl(field)))) then
          read(field,*) pes_forces(i,3*(j-1)+(k-(k_start-1)))
        else
          call error("get_siesta_output: wrong field """//&
            &trim(field)//""" while reading atomic forces")
        end if
      end do
    end do
  end if

  ! Convert forces from eV/Ang to Hartree/Ang -------------
!  pes_forces(i,:)=pes_forces(i,:)*(1.0_DBL/AU_ON_EV)

  ! close unit --------------------------------------------
  close(unit=fnumb_out,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("get_siesta_output: "//trim(err_msg))
  end if
 
end subroutine get_siesta_output

end program main
