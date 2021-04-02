module utility

  !==================================================================
  !   Utility Module
  !==================================================================
  !   It contains variables and functions definitions
  ! used by every other module.
  !==================================================================

#ifdef USE_MPI
  use mpi
#endif
  use iso_fortran_env

  implicit none
  save
  private

  ! public parameters -------------------------------------
  public    :: STDOUT,                     &
               FILEOUT,                    &
               SGL,                        &
               DBL,                        &
               QDR,                        &
               LONG,                       &
               AVOGADRO,                   &
               BOHR_ON_ANG,                &
               AU_ON_J,                    &
               EV_ON_J,                    &
               AU_ON_EV,                   &
               PES_MODE,                   &
               IDPP_MODE,                  &
               flag_mpi
#ifdef USE_MPI
  public    :: MMPI_MSG_END_MAIN_EXEC,     &
               MMPI_MSG_INIT_IMAGES,       &
               MMPI_MSG_INIT_PES_MODULE,   &
               MMPI_MSG_COMPUTE_PES_FORCES
#endif
  ! protected variables -----------------------------------
  public    :: main_program_name,          &
               flag_fileout,               &
               start_clock,                &
               clock_rate
  protected :: main_program_name,          &
               flag_fileout,               &
               start_clock,                &
               clock_rate
#ifdef USE_MPI
  public    :: comm_sz,                    &
               proc_id
  protected :: comm_sz,                    &
               proc_id
#endif
  ! public procedures -------------------------------------
  public    :: set_main_program_name,      &
               set_fileout,                &
               close_fileout,              &
               isdigit,                    &
               isalpha,                    &
               isinteger,                  &
               isreal,                     &
               tolower,                    &
               norm,                       &
               alltrue,                    &
               allfalse,                   &
               get_minima,                 &
               get_maxima,                 &
               get_field,                  &
               get_lines,                  &
               triang_numb,                &
               write_date,                 &
               set_start_clock,            &
               human_time,                 &
               end_main_exec,              &
               error
#ifdef USE_MPI
  public    :: update_comm_sz,             &
               update_proc_id
#endif

  !--------------------------------------------------------
  integer, parameter   :: STDOUT                      = OUTPUT_UNIT
  integer, parameter   :: FILEOUT                     = 110
  integer, parameter   :: SGL                         = REAL32
  integer, parameter   :: DBL                         = REAL64
  integer, parameter   :: QDR                         = REAL128
  integer, parameter   :: LONG                        = INT64
  ! physical parameters taken from :
  ! https://physics.nist.gov/cuu/Constants/
  real(DBL), parameter :: AVOGADRO                    = 6.02214076E23_DBL 
  real(DBL), parameter :: BOHR_ON_ANG                 = 0.529177210903_DBL
  real(DBL), parameter :: AU_ON_J                     = 4.3597447222071E-18_DBL
  real(DBL), parameter :: EV_ON_J                     = 1.602176634E-19_DBL
  real(DBL), parameter :: AU_ON_EV                    = AU_ON_J / EV_ON_J

  ! ENUM
  integer, parameter   :: PES_MODE                    = 0
  integer, parameter   :: IDPP_MODE                   = 1

#ifdef USE_MPI
  logical, parameter   :: flag_mpi                    = .true.
  ! ENUM
  integer, parameter   :: MMPI_MSG_END_MAIN_EXEC      = 0
  integer, parameter   :: MMPI_MSG_INIT_IMAGES        = 1
  integer, parameter   :: MMPI_MSG_INIT_PES_MODULE    = 2
  integer, parameter   :: MMPI_MSG_COMPUTE_PES_FORCES = 3
#else
  logical, parameter   :: flag_mpi                    = .false.
#endif

  character(120)       :: main_program_name
  integer(LONG)        :: start_clock
  integer(LONG)        :: clock_rate
  logical              :: flag_start_clock            = .false.
  logical              :: flag_fileout                = .false.
#ifdef USE_MPI
  integer              :: comm_sz
  integer              :: proc_id
#endif

contains

!====================================================================
! Public
!====================================================================

subroutine set_main_program_name(str)

  character(*), intent(IN) :: str

  main_program_name = str

end subroutine set_main_program_name

!====================================================================

#ifdef USE_MPI
subroutine update_comm_sz(sz)

  integer, intent(IN) :: sz

  logical, save       :: first_call = .true.

  if (first_call.eqv..false.) then
    call error("update_comm_sz: subroutine called more than once")
  end if

  comm_sz = sz

  first_call = .false.

end subroutine update_comm_sz
#endif

!====================================================================

#ifdef USE_MPI
subroutine update_proc_id(id)

  integer, intent(IN) :: id

  logical, save       :: first_call = .true.

  if (first_call.eqv..false.) then
    call error("update_proc_id: subroutine called more than once")
  end if

  proc_id = id

  first_call = .false.

end subroutine update_proc_id
#endif

!====================================================================

subroutine set_fileout(fname_in)

  !--------------------------------------------------------
  ! Takes the input file name, generates a name
  ! for the output file, and opens a stream to that file.
  !--------------------------------------------------------

  character(*), intent(IN)  :: fname_in

  character(80)             :: fname_out
  integer                   :: flen
  integer                   :: err_n
  character(120)            :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_fileout) then
    call error("set_fileout: output file already setted")
  end if

  ! working section ---------------------------------------
  flen = len_trim(fname_in)

  if ((flen-3>0).and.(fname_in(flen-3:flen)==".dat")) then
    fname_out = fname_in(:flen-3)//"out"
  else if ((flen-2>0).and.(fname_in(flen-2:flen)==".in")) then
    fname_out = fname_in(:flen-2)//"out"
  else
    fname_out = trim(fname_in)//".out"
  end if

  open(unit=FILEOUT,file=fname_out,status='REPLACE',action='WRITE',&
    &iostat=err_n,iomsg=err_msg,position='REWIND')
  if (err_n/=0) then
    call error("set_fileout: "//trim(err_msg))
  end if

  flag_fileout = .true.

end subroutine set_fileout

!====================================================================

subroutine close_fileout()

  integer        :: err_n
  character(120) :: err_msg

  close(unit=FILEOUT,iostat=err_n,iomsg=err_msg)
  if (err_n/=0) then
    call error("close_fileout: "//trim(err_msg))
  end if

end subroutine close_fileout

!====================================================================

logical function isdigit(c)

  character, intent(IN) :: c

  if ((c>='0').and.(c<='9')) then
    isdigit = .true.
  else
    isdigit = .false.
  end if

end function isdigit

!====================================================================

logical function isalpha(ch)

  character, intent(IN) :: ch

  character             :: c

  c = ch
  call tolower(c)

  if ((c>='a').and.(c<='z')) then
    isalpha = .true.
  else
    isalpha = .false.
  end if

end function isalpha

!====================================================================

logical function isinteger(str)

  character(*), intent(IN) :: str

  integer                  :: i
  logical                  :: first_digit

  if (len_trim(str)==0) then
    isinteger = .false.
    return
  end if

  first_digit = .false.
  isinteger   = .true.

  do i=1, len(str)
    if (i==1) then
      if (isdigit(str(i:i))) then
        first_digit = .true.
      else if (.not.((str(i:i)=="+").or.(str(i:i)=="-"))) then
        isinteger = .false.
        exit
      end if
    else
      if (isdigit(str(i:i))) then
        first_digit = .true.
      else
        isinteger = .false.
        exit
      end if
    end if
  end do

  if (first_digit.eqv..false.) then
    isinteger = .false.
  end if

end function isinteger

!====================================================================

logical function isreal(str)

  character(*), intent(IN) :: str

  integer                  :: i
  integer                  :: exp_pos
  logical                  :: first_exp
  logical                  :: first_dot
  logical                  :: first_digit

  if (len_trim(str)==0) then
    isreal = .false.
    return
  end if

  first_exp   = .false.
  first_dot   = .false.
  first_digit = .false.
  isreal      = .true.

  do i=1, len(str)
    if (i==1) then
      if (isdigit(str(i:i))) then
        first_digit = .true.
      else if (str(i:i)==".") then
        first_dot = .true.
      else if(.not.((str(i:i)=="+").or.(str(i:i)=="-"))) then
        isreal = .false.
        exit
      end if
    else if (i==len(str)) then
      if (isdigit(str(i:i))) then
        first_digit = .true.      
      else if (str(i:i)==".") then
        if (first_dot) then
          isreal = .false.
          exit
        else if (first_exp) then
          isreal = .false.
          exit
        else
          first_dot = .true.
        end if
      else
        isreal = .false.
        exit
      end if
    else
      if (isdigit(str(i:i))) then
        first_digit = .true.
      else if (str(i:i)==".") then
        if (first_dot) then
          isreal = .false.
          exit
        else if (first_exp) then
          isreal = .false.
          exit
        else
          first_dot = .true.
        end if
      else if ((str(i:i)=="e").or.(str(i:i)=="E")) then
        if (first_exp) then
          isreal = .false.
          exit
        else
          exp_pos   = i
          first_exp = .true.
        end if
      else if ((str(i:i)=="+").or.(str(i:i)=="-")) then
        if (first_exp.eqv..false.) then
          isreal = .false.
          exit
        else if (i/=exp_pos+1) then
          isreal = .false.
          exit
        end if
      else
        isreal = .false.
        exit
      end if
    end if
  end do

  if (first_digit.eqv..false.) then
    isreal = .false.
  end if

end function isreal

!====================================================================

subroutine tolower(str)

  !--------------------------------------------------------
  ! Converts all uppercase letters in str in lowercase ones.
  !--------------------------------------------------------

  character(*), intent(INOUT) :: str

  integer                     :: i

  do i=1, len_trim(str)
    if ((str(i:i)>='A').and.(str(i:i)<='Z')) then
      str(i:i) = achar(iachar(str(i:i))+32)
    end if
  end do

end subroutine tolower

!====================================================================

real(DBL) function norm(arr)

  !--------------------------------------------------------
  ! Computes the norm of arr array.
  !--------------------------------------------------------

  real(DBL), dimension(:), intent(IN) :: arr

  norm = sqrt(sum(arr*arr))

end function norm

!====================================================================

logical function alltrue(arr)

  !--------------------------------------------------------
  ! Returns true if all elements in arr are true,
  ! false otherwise.
  !--------------------------------------------------------

  logical, dimension(:), intent(IN) :: arr

  integer                           :: i
  integer                           :: n

  n = size(arr,1)

  do i=1, n
    if (arr(i).eqv..false.) then
      alltrue = .false.
      return
    end if
  end do

  alltrue = .true.

end function alltrue

!====================================================================

logical function allfalse(arr)

  !--------------------------------------------------------
  ! Returns true if all elements in arr are false,
  ! false otherwise.
  !--------------------------------------------------------

  logical, dimension(:), intent(IN) :: arr

  integer                           :: i
  integer                           :: n

  n = size(arr,1)

  do i=1, n
    if (arr(i).eqv..true.) then
      allfalse = .false.
      return
    end if
  end do

  allfalse = .true.

end function allfalse

!====================================================================

subroutine get_minima(vals,flags)

  real(DBL), dimension(:), intent(IN)  :: vals
  logical,   dimension(:), intent(OUT) :: flags

  integer                              :: i
  integer                              :: n
  real(DBL)                            :: prev
  real(DBL)                            :: curr
  real(DBL)                            :: next

  n = size(vals,1)

  ! preliminary checks ------------------------------------
  if (n/=size(flags,1)) then
    call error("get_minima: wrong arguments' size")
  end if

  ! working section ---------------------------------------
  flags = .false.

  do i=2, n-1
    prev = vals(i-1)
    curr = vals(i)
    next = vals(i+1)

    if ((curr<prev).and.(curr<next)) then
      flags(i) = .true.
    end if
  end do

end subroutine get_minima

!====================================================================

subroutine get_maxima(vals,flags)

  real(DBL), dimension(:), intent(IN)  :: vals
  logical,   dimension(:), intent(OUT) :: flags

  integer                              :: i
  integer                              :: n
  real(DBL)                            :: prev
  real(DBL)                            :: curr
  real(DBL)                            :: next

  n = size(vals,1)

  ! preliminary checks ------------------------------------
  if (n/=size(flags,1)) then
    call error("get_maxima: wrong arguments' size")
  end if

  ! working section ---------------------------------------
  flags = .false.

  do i=2, n-1
    prev = vals(i-1)
    curr = vals(i)
    next = vals(i+1)

    if ((curr>prev).and.(curr>next)) then
      flags(i) = .true.
    end if
  end do

end subroutine get_maxima

!====================================================================

subroutine get_field(str_in,str_out,n,err_n,err_msg)

  !--------------------------------------------------------
  ! Takes an input string str_in in which
  ! every field is separated one another by one or more spaces.
  ! Gets the n field and store it in str_out.
  !--------------------------------------------------------

  character(*), intent(IN)  :: str_in
  character(*), intent(OUT) :: str_out
  integer,      intent(IN)  :: n
  integer,      intent(OUT) :: err_n         ! 0 on success, 1 otherwise
  character(*), intent(OUT) :: err_msg       ! message set in case of failure

  integer, parameter        :: SUCCESS       = 0
  integer, parameter        :: FAILURE       = 1
  integer                   :: i
  integer                   :: current_field
  integer                   :: start_field
  integer                   :: end_field
  character(1)              :: ch
  character(3)              :: n_str
  character(8)              :: istr1
  character(8)              :: istr2
  logical                   :: prev_space

  if (n<1) then
    err_n   = FAILURE
    err_msg = "get_field: argument must be a non-zero positive integer"
    return
  end if

  start_field   = -1
  end_field     = -1
  current_field =  0
  prev_space    = .true.
  do i=1, len_trim(str_in)
    ch = str_in(i:i)

    if (ch==" ") then
      if ((prev_space.eqv..false.).and.(current_field==n)) then
        end_field = i-1
        exit
      end if
      prev_space = .true.
    else
      if (prev_space.eqv..true.) then
        current_field = current_field+1
        if (current_field==n) then
          start_field = i
        end if
      end if

      prev_space = .false.
    end if
  end do

  if (start_field/=-1) then
    if (end_field==-1) then
      end_field = len_trim(str_in)
    end if
  else
    write(n_str,'(I3)') n
    n_str   = adjustl(n_str)
    err_n   = FAILURE
    err_msg = "get_field: cannot get field "//trim(n_str)//&
      &" from string """//trim(str_in)//""""
    return
  end if

  if ((end_field-start_field+1)>len(str_out)) then
    err_n   = FAILURE
    write(istr1,'(I8)') len(str_out)
    istr1 = adjustl(istr1)
    write(istr2,'(I8)') end_field-start_field+1
    istr2 = adjustl(istr2)
    err_msg = "get_field: output string too small ("//trim(istr1)//&
      &") to contain the field ("//trim(istr2)//")"
    return
  end if

  str_out = str_in(start_field:end_field)
  err_n   = SUCCESS

end subroutine get_field

!====================================================================

integer function get_lines(fnumb,ending)

  !--------------------------------------------------------
  ! Returns the number of lines until ending string
  ! is encountered.
  ! On exit, the reading buffer points to the ending string.
  !--------------------------------------------------------

  integer,      intent(IN) :: fnumb
  character(*), intent(IN) :: ending

  logical                  :: is_open
  character(120)           :: str
  integer                  :: err_n
  character(120)           :: err_msg

  inquire(unit=fnumb,opened=is_open)
  if (.not.is_open) then
    call error("get_lines: input stream not opened")
  end if

  get_lines = 0
  do
    read(fnumb,'(A120)',iostat=err_n) str
    if (err_n/=0) then
      call error("get_lines: ending string """//ending//""" not found")
    end if

    if (str==ending) then
      backspace(unit=fnumb,iostat=err_n,iomsg=err_msg)
      if (err_n/=0) then
        call error("get_lines: "//trim(err_msg))
      end if
      exit
    end if

    get_lines = get_lines+1
  end do

end function get_lines

!====================================================================

subroutine write_date(str)

  character(*),               intent(IN) :: str

  character(3), dimension(12), parameter :: month = &
    &[ "Jan", "Feb", "Mar", "Apr", "May", "Jun",&
    &  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" ]
  integer, dimension(8)                  :: v
  integer                                :: fnumb
  logical                                :: is_open

  inquire(unit=FILEOUT,opened=is_open)
  if (is_open) then
    fnumb = FILEOUT
  else
    fnumb = STDOUT
  end if

  call date_and_time(values=v)
  if ((v(2)<1).or.(v(2)>12)) then
    call error("write_date: cannot get the date")
  end if

  write(fnumb,'(1X,"DAT",1X,A,1X,I2,1X,A,1X,I4,1X,I2.2,":",I2.2,":",I2.2)') &
    &trim(adjustl(str)),v(3),month(v(2)),v(1),v(5),v(6),v(7)

end subroutine write_date

!====================================================================

integer function triang_numb(i)

  integer, intent(IN) :: i

  if (i<1) then
    call error("triang_numb: argument less than 1")
  end if

  triang_numb = (i*(i+1))/2

end function triang_numb

!====================================================================

subroutine set_start_clock()

  !--------------------------------------------------------
  ! Set start_clock, flag_start_clock and clock_rate
  !--------------------------------------------------------

  if (flag_start_clock) then
    call error("set_start_clock: start clock already setted")
  end if

  call system_clock(start_clock,clock_rate)

  flag_start_clock = .true.

end subroutine set_start_clock

!====================================================================

subroutine human_time(str,scnt,ecnt,cnt_rate)

  character(*),  intent(IN) :: str
  integer(LONG), intent(IN) :: scnt
  integer(LONG), intent(IN) :: ecnt
  integer(LONG), intent(IN) :: cnt_rate

  integer(LONG)             :: msecs
  integer(LONG)             :: secs
  integer(LONG)             :: mins
  integer(LONG)             :: hours
  integer(LONG)             :: days
  integer                   :: fnumb
  logical                   :: is_open

  inquire(unit=FILEOUT,opened=is_open)
  if (is_open) then
    fnumb = FILEOUT
  else
    fnumb = STDOUT
  end if

  msecs = (ecnt-scnt)/(cnt_rate/1000)
  secs  = msecs/1000
  msecs = msecs-(1000*secs)
  mins  = secs/60
  secs  = secs-(60*mins)
  hours = mins/60
  mins  = mins-(60*hours)
  days  = hours/24
  hours = hours-(24*days)

  write(fnumb,'(" CLK ",A,": ",I5,"d ",I2,"h ",I2,"m ",I2,".",I3.3,"s")')&
    &trim(str), days, hours, mins, secs, msecs

end subroutine human_time

!====================================================================

subroutine end_main_exec()

#ifdef USE_MPI
  integer :: cmd
  integer :: err_n

  if ((proc_id==0).and.(comm_sz>1)) then
    cmd = MMPI_MSG_END_MAIN_EXEC
    call mpi_bcast(cmd,1,MPI_INTEGER,0,MPI_COMM_WORLD,err_n)
  end if

  call MPI_Finalize(err_n)
#endif

  stop   ! end_main_exec

end subroutine end_main_exec

!====================================================================

subroutine error(err_msg)

  character(*), intent(IN) :: err_msg

  integer(LONG)            :: end_clock
  integer                  :: fnumb
  logical                  :: is_open

  inquire(unit=FILEOUT,opened=is_open)
  if (is_open) then
    fnumb = FILEOUT
  else
    fnumb = STDOUT
  end if

  write(fnumb,*) "ERR ", trim(err_msg)
  write(fnumb,*) "ERR Aborting Execution"
  if (flag_start_clock) then
    call system_clock(end_clock)
    call human_time("Total Time Before Error",start_clock,end_clock,clock_rate)
  end if
  call write_date("Error Encountered on")
  stop 1 ! error

end subroutine error

!====================================================================

end module utility

