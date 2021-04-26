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

module output

  use utility
  use geometry
  use pes
  use elastic

  implicit none
  save
  private

  ! public procedures -------------------------------------
  public :: input_dumping,               &
            write_all_images,            &
            write_pes_energy,            &
            write_compare_pes_energy,    &
            write_pes_forces,            &
            write_parallelization_info,  &
            write_pes_info,              &
            write_total_forces,          &
            write_parall_elastic_forces, &
            write_perpen_pes_forces,     &
            write_gnuplot_pes_energy,    &
            last_geom_bkp,               &
            write_transition_state,      &
            write_neb_input_template

contains

!====================================================================
! Public
!====================================================================

subroutine input_dumping(fname_in)

  character(*), intent(IN) :: fname_in

  integer, parameter       :: fnumb_in = 100
  character(200)           :: str
  integer                  :: err_n
  character(120)           :: err_msg

  ! open input file ---------------------------------------
  open(unit=fnumb_in,file=fname_in,status='OLD',action='READ',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("input_dumping: "//trim(err_msg))
  end if

  ! dump the imput ----------------------------------------
  write(FILEOUT,'(1X,"INP ",14("input"))')
  do
    read(fnumb_in,'(A200)',iostat=err_n) str

    if (err_n/=0) then
      exit
    end if

    write(FILEOUT,*) trim(str)
  end do
  write(FILEOUT,'(1X,"INP ",14("input"))')

  ! close input file --------------------------------------
  close(unit=fnumb_in,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("input_dumping: "//trim(err_msg))
  end if

end subroutine input_dumping

!====================================================================

subroutine write_all_images(fnumb,write_elabel)

  integer,     intent(IN) :: fnumb
  logical,     intent(IN) :: write_elabel

  character(*), parameter :: my_name = "write_all_images"
  integer                 :: i
  logical                 :: is_open

  inquire(unit=fnumb,opened=is_open)
  if (is_open.eqv..false.) then
    call error(my_name//": output stream not opened")
  end if

  if (fnumb==FILEOUT) then
    write(fnumb,*) "**  Images Geometry:"
  end if

  do i=0, image_n+1
    call write_image(i,fnumb,write_elabel)
  end do

end subroutine write_all_images

!====================================================================

subroutine write_pes_energy()

  integer :: i

  if (.not.allocated(pes_energy)) then
    call error("write_pes_energy: energy array not allocated")
  end if

  write(FILEOUT,*) "**  Pes Energy:"
  do i=0, image_n+1
    if (i==0) then
      write(FILEOUT,'(5X,A)',advance='no') "Start     :"
    else if (i==image_n+1) then
      write(FILEOUT,'(5X,A)',advance='no') "End       :"
    else
      write(FILEOUT,'(5X,A,I3,A)',advance='no') "Image ",i," :"
    end if
    write(FILEOUT,*) pes_energy(i)
  end do

end subroutine write_pes_energy

!====================================================================

subroutine write_compare_pes_energy(arr)

  real(DBL), dimension(:), intent(IN) :: arr

  integer                             :: i

  if (.not.allocated(pes_energy)) then
    call error("write_compare_pes_energy: energy array not allocated")
  end if

  if (size(arr,1)/=image_n) then
    call error("write_compare_pes_energy: wrong argument size")
  end if

  write(FILEOUT,*) "**  Compare Energies"
  write(FILEOUT,'(5X,"Image",5X,"Real",12X,"Approx")')
  do i=1, image_n
    write(FILEOUT,'(6X,I3,3X,F15.6,1X,F15.6)') i,pes_energy(i),arr(i)
  end do

end subroutine write_compare_pes_energy

!====================================================================

subroutine write_pes_forces()

  integer :: i
  integer :: j

  if (.not.(allocated(pes_forces))) then
    call error("write_pes_forces: forces array not allocated")
  end if

  write(FILEOUT,*) "**  Pes Forces:"
  do i=1, image_n
    write(FILEOUT,'(5X,"Image",1X,I3)') i
    do j=1, geom_len
      write(FILEOUT,'(7X,F16.9)',advance='no') pes_forces(i,j)
      if (mod(j,3)==0) then
        write(FILEOUT,*)
      end if
    end do
    write(FILEOUT,*)
  end do

end subroutine write_pes_forces

!====================================================================

subroutine write_parallelization_info()

  call write_procs_info()

end subroutine write_parallelization_info

!====================================================================

subroutine write_pes_info()

  write(FILEOUT,*) "**  Pes Info:"
  write(FILEOUT,'(5X,"pes_program : ",A)') trim(pes_program)

end subroutine write_pes_info

!====================================================================

subroutine write_total_forces()

  integer :: i
  integer :: j

  if (.not.(allocated(total_forces))) then
    call error("write_total_forces: forces array not allocated")
  end if

  write(FILEOUT,*) "**  Total Forces:"
  do i=1, image_n
    write(FILEOUT,'(5X,"Image",1X,I3)') i
    do j=1, geom_len
      write(FILEOUT,'(7X,F16.9)',advance='no') total_forces(i,j)
      if (mod(j,3)==0) then
        write(FILEOUT,*)
      end if
    end do
    write(FILEOUT,*)
  end do

end subroutine write_total_forces

!====================================================================

subroutine write_parall_elastic_forces()

  integer :: i
  integer :: j

  if (.not.(allocated(parall_elastic_forces))) then
    call error("write_parall_elastic_forces: forces array not allocated")
  end if

  write(FILEOUT,*) "**  Parallel Elastic Forces:"
  do i=1, image_n
    write(FILEOUT,'(5X,"Image",1X,I3)') i
    do j=1, geom_len
      write(FILEOUT,'(7X,F16.9)',advance='no') parall_elastic_forces(i,j)
      if (mod(j,3)==0) then
        write(FILEOUT,*)
      end if
    end do
    write(FILEOUT,*)
  end do

end subroutine write_parall_elastic_forces

!====================================================================

subroutine write_perpen_pes_forces()

  integer :: i
  integer :: j

  if (.not.(allocated(perpen_pes_forces))) then
    call error("write_perpen_pes_forces: forces array not allocated")
  end if

  write(FILEOUT,*) "**  Perpendicular Pes Forces:"
  do i=1, image_n
    write(FILEOUT,'(5X,"Image",1X,I3)') i
    do j=1, geom_len
      write(FILEOUT,'(7X,F16.9)',advance='no') perpen_pes_forces(i,j)
      if (mod(j,3)==0) then
        write(FILEOUT,*)
      end if
    end do
    write(FILEOUT,*)
  end do

end subroutine write_perpen_pes_forces

!====================================================================

subroutine write_gnuplot_pes_energy(n,fname)

  integer,      optional, intent(IN) :: n
  character(*), optional, intent(IN) :: fname

  integer, parameter                 :: gp_fnumb = 600
  character(8)                       :: n_str
  integer                            :: i
  character(30)                      :: gp_fname
  logical                            :: gp_file_exist
  integer                            :: err_n
  character(120)                     :: err_msg

  if (.not.(allocated(pes_energy))) then
    call error("write_gnuplot_pes_energy: energy array not allocated")
  end if

  if ((present(fname)).and.(len_trim(fname)/=0)) then
    gp_fname = fname
  else
    gp_fname = "gnuplot.out"
  end if

  inquire(file=gp_fname,exist=gp_file_exist)
  if (gp_file_exist) then
    open(unit=gp_fnumb,file=gp_fname,status='OLD',action='WRITE',&
      &iostat=err_n,iomsg=err_msg,position='APPEND')
    if (err_n/=0) then
      call error("write_gnuplot_pes_energy: "//trim(err_msg))
    end if
  else
    open(unit=gp_fnumb,file=gp_fname,status='NEW',action='WRITE',&
      &iostat=err_n,iomsg=err_msg,position='REWIND')
    if (err_n/=0) then
      call error("write_gnuplot_pes_energy: "//trim(err_msg))
    end if
  end if

  if (present(n)) then
    write(n_str,'(I8)') n
    write(gp_fnumb,'("#",A)') trim(adjustl(n_str))
  end if

  do i=0, image_n+1
    write(gp_fnumb,'(1X,I4,3X,F20.9)') i, pes_energy(i)
  end do

  do i=1, 2
    write(gp_fnumb,*)
  end do

  close(unit=gp_fnumb,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("write_gnuplot_pes_energy: "//trim(err_msg))
  end if

end subroutine write_gnuplot_pes_energy

!====================================================================

subroutine last_geom_bkp(write_elabel,fname)

  logical,                intent(IN) :: write_elabel
  character(*), optional, intent(IN) :: fname

  character(*), parameter            :: my_name   = "last_geom_bkp"
  integer,      parameter            :: bkp_fnumb = 610
  character(30)                      :: bkp_fname
  integer                            :: err_n
  character(120)                     :: err_msg

  if (present(fname).and.(len_trim(fname)>0)) then
    bkp_fname = fname
  else
    bkp_fname = "lastgeom.bkp"
  end if

  open(unit=bkp_fnumb,file=bkp_fname,status='REPLACE',action='WRITE',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  call write_all_images(bkp_fnumb,write_elabel)

  close(unit=bkp_fnumb,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

end subroutine last_geom_bkp

!====================================================================

subroutine write_transition_state()

  integer                            :: i
  integer                            :: n
  integer                            :: tot_ts
  logical, allocatable, dimension(:) :: indx
  real(DBL)                          :: curr
  real(DBL)                          :: prev
  real(DBL)                          :: next
  integer                            :: err_n
  character(120)                     :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_init_images.eqv..false.) then
    call error("write_transition_state: images not initialized")
  end if

  if (.not.allocated(pes_energy)) then
    call error("write_transition_state: energy array not allocated")
  end if

  ! allocation section ------------------------------------
  allocate(indx(image_n),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("write_transition_state: "//trim(err_msg))
  end if

  ! init variables ----------------------------------------
  indx   = .false.
  tot_ts = 0
  do i=1, image_n
    prev = pes_energy(i-1)
    curr = pes_energy(i)
    next = pes_energy(i+1)
    if ((curr>prev).and.(curr>next)) then
      indx(i) = .true.
      tot_ts  = tot_ts+1
    end if
  end do

  ! write transition state(s) -----------------------------
  write(FILEOUT,*) "**  Transition State"
  if (tot_ts==0) then
    write(FILEOUT,'(5X,"No Transition State Founded")')
  else if (tot_ts==1) then
    write(FILEOUT,'(5X,"Founded ",I3," Transition State")') tot_ts
  else
    write(FILEOUT,'(5X,"Founded ",I3," Transition States")') tot_ts
  end if

  n = 0
  do i=1, image_n
    if (indx(i)) then
      n = n+1
      write(FILEOUT,'(5X,"TS ",I3)') n
      call write_image(i,FILEOUT,.false.)
      call write_delta_e(i)
    end if
  end do

  ! deallocation section ----------------------------------
  deallocate(indx,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("write_transition_state: "//trim(err_msg))
  end if

end subroutine write_transition_state

!====================================================================

subroutine write_neb_input_template(verbose)

  logical,     intent(IN) :: verbose

  integer, parameter      :: tmplt_fnumb = 620
  character(*), parameter :: tmplt_fname = "data.in"
  integer                 :: err_n
  character(120)          :: err_msg

  ! open template file ------------------------------------
  open(unit=tmplt_fnumb,file=tmplt_fname,status='REPLACE',action='WRITE',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')

  ! write template ----------------------------------------
  write(tmplt_fnumb,'(A)') "! This is an automatically generated, fully functional input template"
  write(tmplt_fnumb,'(A)') "! that performs a NEB computation by calling the external program Gaussian."
  write(tmplt_fnumb,'(A)') "! Feel free to modify this file for your own purposes."
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! The NEB input file is processed by means of keywords,"
  write(tmplt_fnumb,'(A)') "! that are special words beginning with the ""#"" character."
  write(tmplt_fnumb,'(A)') "! Comment lines begin with the ""!"" character."
  write(tmplt_fnumb,'(A)') "! Comments are useful for obscuring keywords without having to delete them."
  write(tmplt_fnumb,'(A)')
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! If the #ONLYINTERPOLATION keyword is specified,"
  write(tmplt_fnumb,'(A)') "! the program will write only the files containing"
  write(tmplt_fnumb,'(A)') "! the interpolated geometries of the images,"
  write(tmplt_fnumb,'(A)') "! without actually executing the NEB method."
  end if
  write(tmplt_fnumb,'(A)') "!#ONLYINTERPOLATION"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #OPTALGORITHM keyword, followed by a string,"
  write(tmplt_fnumb,'(A)') "! specifies what algorithm to use for the elastic band optimization."
  write(tmplt_fnumb,'(A)') "! Accepted arguments are ""sd"" (for the steepest descent),"
  write(tmplt_fnumb,'(A)') "! and ""fire"" (default, for the fast inertial relaxation engine)."
  end if
  write(tmplt_fnumb,'(A)') "!#OPTALGORITHM fire"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #OPTCYCLE keyword, followed by an integer, specifies"
  write(tmplt_fnumb,'(A)') "! the max optimization cycles that will be performed."
  write(tmplt_fnumb,'(A)') "! A negative argument means infinite cycles until convergence."
  write(tmplt_fnumb,'(A)') "! Despite this keyword is optional, you are strongly advised"
  write(tmplt_fnumb,'(A)') "! to specify it, as optimization cycles number varies between"
  write(tmplt_fnumb,'(A)') "! optimization algorithms and program versions."
  end if
  write(tmplt_fnumb,'(A)') "#OPTCYCLE -1"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #OPTCONV keyword, followed by a real, specifies"
  write(tmplt_fnumb,'(A)') "! the convergence threshold that must be achieved"
  write(tmplt_fnumb,'(A)') "! by all images for a successful execution."
  write(tmplt_fnumb,'(A)') "! Despite this keyword is optional, you are strongly advised"
  write(tmplt_fnumb,'(A)') "! to specify it, as optimization cycles number varies between"
  write(tmplt_fnumb,'(A)') "! optimization algorithms and program versions."
  end if
  write(tmplt_fnumb,'(A)') "#OPTCONV 1.0E-3"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! If the #IDPP keyword is specified,"
  write(tmplt_fnumb,'(A)') "! an interpolation technique called"
  write(tmplt_fnumb,'(A)') "! Image Dependent Pair Potential will be performed."
  write(tmplt_fnumb,'(A)') "! This is way better than the standard linear interpolation,"
  write(tmplt_fnumb,'(A)') "! so we strongly advise to use it."
  end if
  write(tmplt_fnumb,'(A)') "#IDPP"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #IDPPCONV keyword, followed by a real, specifies"
  write(tmplt_fnumb,'(A)') "! the convergence threshold that must be achieved"
  write(tmplt_fnumb,'(A)') "! during the IDPP interpolation technique."
  write(tmplt_fnumb,'(A)') "! Default value is 1.0E-3, and it can be increased"
  write(tmplt_fnumb,'(A)') "! in case the interpolation does not converge the first time."
  end if
  write(tmplt_fnumb,'(A)') "!#IDPPCONV 1.0E-3"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory], conflicts with #GEOMETRIESFILE"
  write(tmplt_fnumb,'(A)') "! The #IMAGES keyword, followed by an integer,"
  write(tmplt_fnumb,'(A)') "! specifies the number of of images that will be interposed"
  write(tmplt_fnumb,'(A)') "! between the initial and final geometries."
  end if
  write(tmplt_fnumb,'(A)') "#IMAGES 6"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory], conflicts with #START, #END, #IMAGES"
  write(tmplt_fnumb,'(A)') "! The #GEOMETRIESFILE keyword, followed by a string,"
  write(tmplt_fnumb,'(A)') "! specifies the name of the file from which the"
  write(tmplt_fnumb,'(A)') "! starting, ending and images geometries will be read."
  write(tmplt_fnumb,'(A)') "! This file must have the same format as the lastgeom.bkp file."
  write(tmplt_fnumb,'(A)') "! This keyword can be used either to restart a NEB calculation,"
  write(tmplt_fnumb,'(A)') "! or to specify some handcrafted geometries."
  end if
  write(tmplt_fnumb,'(A)') "!#GEOMETRIESFILE geometries.input"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory], conflicts with #GEOMETRIESFILE"
  write(tmplt_fnumb,'(A)') "! The #START keyword specifies the initial geometry of the NEB."
  write(tmplt_fnumb,'(A)') "! Geometry is specified similar to a xyz file."
  end if
  write(tmplt_fnumb,'(A)') "#START"
  write(tmplt_fnumb,'(A)') "4"
  write(tmplt_fnumb,'(A)') "O     0.000000    0.000000    0.000000"
  write(tmplt_fnumb,'(A)') "C    -1.203990   -0.000000    0.000000"
  write(tmplt_fnumb,'(A)') "H    -1.805419    0.945297    0.000000"
  write(tmplt_fnumb,'(A)') "H    -1.805419   -0.945297    0.000000"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory], conflicts with #GEOMETRIESFILE"
  write(tmplt_fnumb,'(A)') "! The #END keyword specifies the final geometry of the NEB."
  end if
  write(tmplt_fnumb,'(A)') "#END"
  write(tmplt_fnumb,'(A)') "4"
  write(tmplt_fnumb,'(A)') "O     0.000000    0.000000    0.000000"
  write(tmplt_fnumb,'(A)') "C    -1.317837   -0.000000    0.000000"
  write(tmplt_fnumb,'(A)') "H    -1.542714    1.108433    0.000000"
  write(tmplt_fnumb,'(A)') "H     0.296335   -0.928935    0.000000"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory]"
  write(tmplt_fnumb,'(A)') "! The #STARTENERGY keyword, followed by a real,"
  write(tmplt_fnumb,'(A)') "! specifies the energy of the starting geometry"
  end if
  write(tmplt_fnumb,'(A)') "#STARTENERGY -114.507640701"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory]"
  write(tmplt_fnumb,'(A)') "! The #ENDENERGY keyword, followed by a real,"
  write(tmplt_fnumb,'(A)') "! specifies the energy of the ending geometry"
  end if
  write(tmplt_fnumb,'(A)') "#ENDENERGY -114.423707780"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #CLIMBING keyword specifies that"
  write(tmplt_fnumb,'(A)') "! the Climbing Image NEB method will be applied."
  write(tmplt_fnumb,'(A)') "! The CI-NEB only affects the images that are energy maxima."
  write(tmplt_fnumb,'(A)') "! The #CLIMBING keyword is followed either by an integer,"
  write(tmplt_fnumb,'(A)') "! that specifies the number of maxima on which apply the CI-NEB variant,"
  write(tmplt_fnumb,'(A)') "! or by the ""all"" string (without quotes), that means that the CI-NEB"
  write(tmplt_fnumb,'(A)') "! will be applied to all the images that are energy maxima."
  end if
  write(tmplt_fnumb,'(A)') "#CLIMBING 1"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! If the #CLIMBINGQUICKSTART keyword is specified,"
  write(tmplt_fnumb,'(A)') "! the CI-NEB will be performed from the very first iteration."
  write(tmplt_fnumb,'(A)') "! It's good practice to declare this keyword"
  write(tmplt_fnumb,'(A)') "! if you want to resume a CI-NEB calculation that ended abnormally."
  end if
  write(tmplt_fnumb,'(A)') "!#CLIMBINGQUICKSTART"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #DESCENDING keyword enable the Descending Image NEB."
  write(tmplt_fnumb,'(A)') "! This means that all the images that are energy minima"
  write(tmplt_fnumb,'(A)') "! will be optimized with respect to the real forces acting on them."
  write(tmplt_fnumb,'(A)') "! Contrary to #CLIMBING, #DESCENDING automatically affects all energy minima."
  end if
  write(tmplt_fnumb,'(A)') "!#DESCENDING"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! If the #DESCENDINGQUICKSTART keyword is specified,"
  write(tmplt_fnumb,'(A)') "! the Descending Image NEB will be performed from the very first iteration."
  write(tmplt_fnumb,'(A)') "! It's good practice to declare this keyword"
  write(tmplt_fnumb,'(A)') "! if you want to resume a NEB calculation that ended abnormally."
  end if
  write(tmplt_fnumb,'(A)') "!#DESCENDINGQUICKSTART"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory]"
  write(tmplt_fnumb,'(A)') "! The #PESPROGRAM keyword specifies the kind of external program"
  write(tmplt_fnumb,'(A)') "! that will be used. Accepted arguments are: ""gaussian"" and ""siesta""."
  end if
  write(tmplt_fnumb,'(A)') "#PESPROGRAM gaussian"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #NEWPESPROGRAM keyword must be used if the user specifies,"
  write(tmplt_fnumb,'(A)') "! via the #PESPROGRAM keywork,"
  write(tmplt_fnumb,'(A)') "! an external program other than those currently supported."
  write(tmplt_fnumb,'(A)') "! Its function is to suppress some checks inside the code."
  write(tmplt_fnumb,'(A)') "! Read the documentation to learn how to write an interface for a new program."
  end if
  write(tmplt_fnumb,'(A)') "!#NEWPESPROGRAM"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory]"
  write(tmplt_fnumb,'(A)') "! The #PESEXEC keyword specifies the executable name"
  write(tmplt_fnumb,'(A)') "! of the external program, as it is called on your system."
  write(tmplt_fnumb,'(A)') "! Please note that the executable must be present inside any directory"
  write(tmplt_fnumb,'(A)') "! specified in the PATH environment variable, otherwise the full path"
  write(tmplt_fnumb,'(A)') "! must be provided (e.g. ""/opt/gaussian09/bin/g09"")."
  end if
  write(tmplt_fnumb,'(A)') "#PESEXEC g09"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #PESPROC keyword, followed by an integer,"
  write(tmplt_fnumb,'(A)') "! specifies the number of processors on which the external program will run."
  write(tmplt_fnumb,'(A)') "! This is needed with siesta, because it's run by mpirun."
  end if
  write(tmplt_fnumb,'(A)') "!#PESPROC 2"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #AUXINPUTFILES is followed by an integer N, and then by N strings."
  write(tmplt_fnumb,'(A)') "! This keyword specifies how many and which files are needed"
  write(tmplt_fnumb,'(A)') "! by the external program (such as the siesta pseudo potentials)."
  write(tmplt_fnumb,'(A)') "! These files must be in the master directory where the neb is run,"
  write(tmplt_fnumb,'(A)') "! and are automatically copied to the working directories"
  write(tmplt_fnumb,'(A)') "! from which the external program is run."
  write(tmplt_fnumb,'(A)') "! The following is a shortcut to copy all the "".psml"" files at once:"
  end if
  write(tmplt_fnumb,'(A)') "!#AUXINPUTFILES 1 *.psml"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #AUXOUTPUTFILES is followed by an integer N, and then by N strings."
  write(tmplt_fnumb,'(A)') "! This keyword specifies how many and which files are useful"
  write(tmplt_fnumb,'(A)') "! for speeding up the next calculation (such as the density matrix)."
  write(tmplt_fnumb,'(A)') "! These files are saved at the end of the calculation,"
  write(tmplt_fnumb,'(A)') "! and then copied back into the working directory of the next calculation."
  write(tmplt_fnumb,'(A)') "! These files are image specific, so the files produced by image M"
  write(tmplt_fnumb,'(A)') "! will be reused only for the subsequent calculations on image M."
  write(tmplt_fnumb,'(A)') "! Here we are saying that we want to save all the files ending in "".DM"":"
  end if
  write(tmplt_fnumb,'(A)') "!#AUXOUTPUTFILES 1 *.DM"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Mandatory], numbers of the blocks depend on the selected external program"
  write(tmplt_fnumb,'(A)') "! Any number of #PESINPUTTEMPLATE blocks can be specified in this file,"
  write(tmplt_fnumb,'(A)') "! to ensure flexibility in the case"
  write(tmplt_fnumb,'(A)') "! of implementation of additional external programs."
  write(tmplt_fnumb,'(A)') "! The #PESINPUTTEMPLATE keyword is followed by an integer"
  write(tmplt_fnumb,'(A)') "! that identifies the block."
  write(tmplt_fnumb,'(A)') "! The block is ended with the #ENDPESINPUTTEMPLATE keyword."
  write(tmplt_fnumb,'(A)') "! All text inside the block will be copied and reused as is."
  end if
  write(tmplt_fnumb,'(A)') "#PESINPUTTEMPLATE 1"
  write(tmplt_fnumb,'(A)') "%nproc=2"
  write(tmplt_fnumb,'(A)') "%mem=2GB"
  write(tmplt_fnumb,'(A)') "#p b3lyp/cc-pvdz"
  write(tmplt_fnumb,'(A)') "#ENDPESINPUTTEMPLATE"
  write(tmplt_fnumb,'(A)')
  write(tmplt_fnumb,'(A)') "#PESINPUTTEMPLATE 2"
  write(tmplt_fnumb,'(A)') "gaussian title"
  write(tmplt_fnumb,'(A)')
  write(tmplt_fnumb,'(A)') "0 1"
  write(tmplt_fnumb,'(A)') "#ENDPESINPUTTEMPLATE"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional / Mandatory if #NEWPESPROGRAM is specified]"
  write(tmplt_fnumb,'(A)') "! The #SCFCYCLE keyword, followed by an integer,"
  write(tmplt_fnumb,'(A)') "! specifies the maximum number of SCF cycles that will be performed."
  write(tmplt_fnumb,'(A)') "! It's the NEB program itself that writes this information"
  write(tmplt_fnumb,'(A)') "! in the input file of the external program,"
  write(tmplt_fnumb,'(A)') "! since some older siesta versions need this information"
  write(tmplt_fnumb,'(A)') "! in order to determine if convergence has been reached."
  write(tmplt_fnumb,'(A)') "! For this reason, the end user must make sure to *not* enter this information"
  write(tmplt_fnumb,'(A)') "! in a #PESINPUTTEMPLATE  block under any circumstances."
  end if
  write(tmplt_fnumb,'(A)') "!#SCFCYCLE 400"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional / Mandatory if #NEWPESPROGRAM is specified]"
  write(tmplt_fnumb,'(A)') "! The #SCFCONV keyword, followed by a real,"
  write(tmplt_fnumb,'(A)') "! specifies the SCF convergence threshold."
  write(tmplt_fnumb,'(A)') "! This information must be known by the NEB program"
  write(tmplt_fnumb,'(A)') "! so that, in case of non-convergence of the external calculation,"
  write(tmplt_fnumb,'(A)') "! it can try to rerun the calculation with a lower threshold."
  write(tmplt_fnumb,'(A)') "! For this reason, the end user must make sure to *not* enter this information"
  write(tmplt_fnumb,'(A)') "! in a #PESINPUTTEMPLATE  block under any circumstances."
  end if
  write(tmplt_fnumb,'(A)') "!#SCFCONV 1.0e-6"
  write(tmplt_fnumb,'(A)')
  if (verbose.eqv..true.) then
  write(tmplt_fnumb,'(A)') "! [Optional]"
  write(tmplt_fnumb,'(A)') "! The #SPRINGMODE, followed by a string, specifies"
  write(tmplt_fnumb,'(A)') "! the type of spring force constants to be used to connect the images."
  write(tmplt_fnumb,'(A)') "! Accepted arguments are ""static"" and ""dynamic"" (default)."
  end if
  write(tmplt_fnumb,'(A)') "!#SPRINGMODE dynamic"

  ! close template file -----------------------------------
  close(unit=tmplt_fnumb,iostat=err_n,iomsg=err_msg)

  write(STDOUT,'(A)') "NEB input template has been written in "//tmplt_fname

end subroutine write_neb_input_template

!====================================================================
! Private
!====================================================================

subroutine write_procs_info()

  integer      :: buff
  character(8) :: n

#ifdef USE_MPI
  buff = comm_sz
#else
  buff = 1
#endif

  write(n,'(I8)') buff
  n = adjustl(n)

  write(FILEOUT,*) "**  Procs Info:"
  write(FILEOUT,'(5X,"Executed with: ",A)',advance="NO") trim(n)
  if (buff==1) then
    write(FILEOUT,'(" Process")')
  else
    write(FILEOUT,'(" Processes")')
  end if

end subroutine write_procs_info

!====================================================================

subroutine write_image(n,fnumb,write_elabel)

  integer,     intent(IN) :: n
  integer,     intent(IN) :: fnumb
  logical,     intent(IN) :: write_elabel

  character(*), parameter :: my_name = "write_image"
  integer                 :: i
  character(8)            :: atoms
  character(120)          :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_init_images.eqv..false.) then
    call error(my_name//": images not initialized")
  else if ((n<0).or.(n>(image_n+1))) then
    write(err_msg,'(A,I5)') my_name//": wrong argument ", n
    call error(err_msg)
  end if

  ! write number of atoms and comment line ----------------
  write(atoms,'(I8)') geom_len/3
  write(fnumb,*) trim(adjustl(atoms))
  if (n==0) then
    write(fnumb,'(" Initial Geometry ")',advance="NO")
  else if (n==image_n+1) then
    write(fnumb,'(" Final Geometry   ")',advance="NO")
  else
    write(fnumb,'(" Image ",I3,8X)',advance="NO") n
  end if

  if (flag_init_pes_module) then
    write(fnumb,'(":: Energy ",F15.6)') pes_energy(n)
  else
    write(fnumb,*)
  end if
  
  ! write the actual geometry -----------------------------
  do i=1, geom_len
    if (mod(i,3)==1) then
      write(fnumb,'(4X,A2,1X)',advance='no') element(i/3+1)
    end if

    write(fnumb,'(1X,F13.6)',advance='no') image_geom(n,i)

    if (mod(i,3)==0) then
      if (write_elabel.and.flag_elabel) then
        write(fnumb,'(2X,A)') elabel(i/3)
      else
        write(fnumb,*)
      end if
    end if
  end do

end subroutine write_image

!====================================================================

subroutine write_delta_e(n)

  integer,     intent(IN) :: n

  real(DBL)               :: de_rea
  real(DBL)               :: de_prod
  real(DBL), dimension(2) :: de_au
  real(DBL), dimension(2) :: de_kjmol
  real(DBL), dimension(2) :: de_ev

  de_rea  = pes_energy(n) - pes_energy(0)
  de_prod = pes_energy(n) - pes_energy(image_n+1)

  select case (pes_program)
  case ("gaussian")
    de_au(1) = de_rea
    de_au(2) = de_prod
    de_kjmol = de_au * AU_ON_J * 1.0E-3_DBL * AVOGADRO
    de_ev    = de_au * AU_ON_EV
  case ("siesta")
    de_ev(1) = de_rea
    de_ev(2) = de_prod
    de_au    = de_ev / AU_ON_EV
    de_kjmol = de_ev * EV_ON_J * 1.0E-3_DBL * AVOGADRO
  case default
    if (flag_new_pes_program) then
      return
    else
      call error("write_delta_e: invalid option """//trim(pes_program)//"""")
    end if
  end select

  write(FILEOUT,'(5X,"Energy barrier, direct reaction  : ",F12.6," a.u. ; ",&
    &F12.1," kJ/mol ; ",F12.2," eV")') de_au(1), de_kjmol(1), de_ev(1)
  write(FILEOUT,'(5X,"Energy barrier, reverse reaction : ",F12.6," a.u. ; ",&
    &F12.1," kJ/mol ; ",F12.2," eV")') de_au(2), de_kjmol(2), de_ev(2)

end subroutine write_delta_e

!====================================================================

end module output

