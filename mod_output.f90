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
            write_transition_state

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

subroutine write_all_images(fnumb)

  integer, optional, intent(IN) :: fnumb

  integer                       :: i
  integer                       :: p_fnumb
  logical                       :: is_open

  if (present(fnumb).and.(fnumb>0)) then
    inquire(unit=fnumb,opened=is_open)
    if (is_open) then
      p_fnumb = fnumb
    else
      call error("write_all_image: output stream not opened")
    end if
  else
    p_fnumb = FILEOUT
  end if

  if (p_fnumb==FILEOUT) then
    write(FILEOUT,*) "**  Images Geometry:"
  end if

  do i=0, image_n+1
    call write_image(i,p_fnumb)
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

subroutine last_geom_bkp(fname)

  character(*), optional, intent(IN) :: fname

  integer, parameter                 :: bkp_fnumb = 610
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
    call error("last_geom_bkp: "//trim(err_msg))
  end if

  call write_all_images(bkp_fnumb)

  close(unit=bkp_fnumb,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("last_geom_bkp: "//trim(err_msg))
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
      call write_image(i,FILEOUT)
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

subroutine write_image(n,fnumb)

  integer, intent(IN) :: n
  integer, intent(IN) :: fnumb

  integer             :: i
  character(8)        :: atoms
  character(120)      :: err_msg

  if (flag_init_images.eqv..false.) then
    call error("write_image: images not initialized")
  else if ((n<0).or.(n>(image_n+1))) then
    write(err_msg,'(A,I5)') "write_image: wrong argument ", n
    call error(err_msg)
  end if

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
  
  do i=1, geom_len
    if (mod(i,3)==1) then
      write(fnumb,'(4X,A2,1X)',advance='no') element(i/3+1)
    end if

    write(fnumb,'(1X,F13.6)',advance='no') image_geom(n,i)

    if (mod(i,3)==0) then
      write(fnumb,*)
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
    call error("write_delta_e: invalid option """//trim(pes_program)//"""")
  end select

  write(FILEOUT,'(5X,"DE_rea  : ",F12.6," a.u. ; ",&
    &F12.1," kJ/mol ; ",F12.2," eV")') de_au(1), de_kjmol(1), de_ev(1)
  write(FILEOUT,'(5X,"DE_prod : ",F12.6," a.u. ; ",&
    &F12.1," kJ/mol ; ",F12.2," eV")') de_au(2), de_kjmol(2), de_ev(2)

end subroutine write_delta_e

!====================================================================

end module output

