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

program driver

  use utility
  use bfgs
  use lbfgs

  implicit none

  call driver_bfgs_rosenbrock()

contains

!====================================================================

subroutine driver_bfgs_rosenbrock()

  character(*), parameter   :: my_name = "driver_bfgs_rosenbrock"
  character(*), parameter   :: header  = "(2X,""Iter"",5X,""x1"",14X,""x2"",14X,""f(x)"",12X, &
                                          &""d/dx1 f(x)"",6X,""d/dx2 f(x)"",6X,""Norm ("",ES8.1,"")"")"
  character(*), parameter   :: sep     = "(105(""-""))"
  character(*), parameter   :: format1 = "(2X,I6,3X,5(ES13.6,3X))"
  character(*), parameter   :: format2 = "(2X,I6,3X,6(ES13.6,3X))"
  character(*), parameter   :: format3 = "(2X,I6,3X,""skipped"")"
  real(DBL), parameter      :: tol     = 1.0E-5_DBL
  character(120)            :: cmdstr
  character(120)            :: mode

  ! for BFGS
  real(DBL), dimension(2,1) :: old_x
  real(DBL)                 :: old_f
  real(DBL), dimension(2,1) :: old_d
  real(DBL), dimension(2,2) :: old_h
  real(DBL), dimension(2,1) :: new_x
  real(DBL)                 :: new_f
  real(DBL), dimension(2,1) :: new_d
  real(DBL), dimension(2,2) :: new_h
  real(DBL), dimension(2)   :: buff

  ! for L-BFGS
  real(DBL), dimension(2)   :: vold_x
  real(DBL)                 :: vold_f
  real(DBL), dimension(2)   :: vold_d
  real(DBL), dimension(2)   :: vnew_x
  real(DBL)                 :: vnew_f
  real(DBL), dimension(2)   :: vnew_d

  real(DBL)                 :: nrm
  integer, dimension(6)     :: rosenbrock_calls
  integer                   :: i

  rosenbrock_calls = 0

  !--------------------------------------------------------

  write(*,sep)
  write(*,*) my_name//": BFGS with fixed alpha step"
  write(*,header) tol

  mode = "fixed_alpha"
  old_x(1,1) = -1.2_DBL
  old_x(2,1) =  1.0_DBL

  i = 1
  do
    if (i == 1) then
      old_h      = 0.0_DBL
      old_h(1,1) = 1.0_DBL
      old_h(2,2) = 1.0_DBL

      call rosenbrock_f(reshape(old_x,(/2/)),old_f)
      rosenbrock_calls(1) = rosenbrock_calls(1) + 1
      call rosenbrock_d(reshape(old_x,(/2/)),buff)
      old_d = reshape(buff,(/2,1/))
      write(*,format1) i, old_x, old_f, old_d

      cmdstr = "START"
      call bfgs_internal(    &
        cmdstr,              &
        mode,                &
        old_x,               &
        new_x,               &
        old_d,               &
        new_d,               &
        old_h,               &
        new_h,               &
        reset_alpha = .true. &
      )
    else
      cmdstr = "START"
      call bfgs_internal( &
        cmdstr,           &
        mode,             &
        old_x,            &
        new_x,            &
        old_d,            &
        new_d,            &
        old_h,            &
        new_h             &
      )
    end if

    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": expected ""EVALUATE_DF1"", get """//trim(cmdstr)//"""")
    end if

    call rosenbrock_f(reshape(new_x,(/2/)),new_f)
    rosenbrock_calls(1) = rosenbrock_calls(1) + 1
    call rosenbrock_d(reshape(new_x,(/2/)),buff)
    new_d = reshape(buff,(/2,1/))

    cmdstr = "EVALUATED"
    call bfgs_internal( &
      cmdstr,           &
      mode,             &
      old_x,            &
      new_x,            &
      old_d,            &
      new_d,            &
      old_h,            &
      new_h             &
    )

    if (cmdstr /= "DONE") then
      call error(my_name//": expected ""DONE"", get """//trim(cmdstr)//"""")
    end if

    nrm = sqrt(sum(matmul(transpose(new_d),new_d)))
    write(*,format2) i, new_x, new_f, new_d, nrm

    if (nrm < tol) then
      exit
    end if

    old_x = new_x
    old_f = new_f
    old_d = new_d
    old_h = new_h

    i = i+1
  end do

  write(*,*) "total rosenbrock_f calls: ", rosenbrock_calls(1)

  !--------------------------------------------------------

  write(*,sep)
  write(*,*) my_name//": BFGS with accelerated_backtracking_line_search()"
  write(*,header) tol

  mode = "accelerated_bt"
  old_x(1,1) = -1.2_DBL
  old_x(2,1) =  1.0_DBL

  i = 1
  do
    if (i == 1) then
      old_h      = 0.0_DBL
      old_h(1,1) = 1.0_DBL
      old_h(2,2) = 1.0_DBL

      call rosenbrock_f(reshape(old_x,(/2/)),old_f)
      rosenbrock_calls(2) = rosenbrock_calls(2) + 1
      call rosenbrock_d(reshape(old_x,(/2/)),buff)
      old_d = reshape(buff,(/2,1/))
      write(*,format1) i, old_x, old_f, old_d

      cmdstr = "START"
      call bfgs_internal(    &
        cmdstr,              &
        mode,                &
        old_x,               &
        new_x,               &
        old_d,               &
        new_d,               &
        old_h,               &
        new_h,               &
        reset_alpha = .true. &
      )
    else
      cmdstr = "START"
      call bfgs_internal( &
        cmdstr,           &
        mode,             &
        old_x,            &
        new_x,            &
        old_d,            &
        new_d,            &
        old_h,            &
        new_h             &
      )
    end if

    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": expected ""EVALUATE_DF1"", get """//trim(cmdstr)//"""")
    end if

    call rosenbrock_f(reshape(new_x,(/2/)),new_f)
    rosenbrock_calls(2) = rosenbrock_calls(2) + 1
    call rosenbrock_d(reshape(new_x,(/2/)),buff)
    new_d = reshape(buff,(/2,1/))

    cmdstr = "EVALUATED"
    call bfgs_internal( &
      cmdstr,           &
      mode,             &
      old_x,            &
      new_x,            &
      old_d,            &
      new_d,            &
      old_h,            &
      new_h             &
    )

    select case (cmdstr)
    case ("DONE")
      nrm = sqrt(sum(matmul(transpose(new_d),new_d)))
      write(*,format2) i, new_x, new_f, new_d, nrm

      if (nrm < tol) then
        exit
      end if

      old_x = new_x
      old_f = new_f
      old_d = new_d
      old_h = new_h
    case ("SKIPPED")
      write(*,format3) i

!      ! reset old_h to identity matrix
!      old_h      = 0.0_DBL
!      old_h(1,1) = 1.0_DBL
!      old_h(2,2) = 1.0_DBL
    case default
      call error(my_name//": expected ""DONE"" or ""SKIPPED"", get """//trim(cmdstr)//"""")
    end select

    i = i+1
  end do

  write(*,*) "total rosenbrock_f calls: ", rosenbrock_calls(2)

  !--------------------------------------------------------

  write(*,sep)
  write(*,*) my_name//": BFGS with max displacement"
  write(*,header) tol

  mode = "max_displacement"
  old_x(1,1) = -1.2_DBL
  old_x(2,1) =  1.0_DBL

  i = 1
  do
    if (i == 1) then
      old_h      = 0.0_DBL
      old_h(1,1) = 1.0_DBL
      old_h(2,2) = 1.0_DBL

      call rosenbrock_f(reshape(old_x,(/2/)),old_f)
      rosenbrock_calls(3) = rosenbrock_calls(3) + 1
      call rosenbrock_d(reshape(old_x,(/2/)),buff)
      old_d = reshape(buff,(/2,1/))
      write(*,format1) i, old_x, old_f, old_d

      cmdstr = "START"
      call bfgs_internal(    &
        cmdstr,              &
        mode,                &
        old_x,               &
        new_x,               &
        old_d,               &
        new_d,               &
        old_h,               &
        new_h,               &
        reset_alpha = .true. &
      )
    else
      cmdstr = "START"
      call bfgs_internal( &
        cmdstr,           &
        mode,             &
        old_x,            &
        new_x,            &
        old_d,            &
        new_d,            &
        old_h,            &
        new_h             &
      )
    end if

    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": expected ""EVALUATE_DF1"", get """//trim(cmdstr)//"""")
    end if

    call rosenbrock_f(reshape(new_x,(/2/)),new_f)
    rosenbrock_calls(3) = rosenbrock_calls(3) + 1
    call rosenbrock_d(reshape(new_x,(/2/)),buff)
    new_d = reshape(buff,(/2,1/))

    cmdstr = "EVALUATED"
    call bfgs_internal( &
      cmdstr,           &
      mode,             &
      old_x,            &
      new_x,            &
      old_d,            &
      new_d,            &
      old_h,            &
      new_h             &
    )

    if (cmdstr /= "DONE") then
      call error(my_name//": expected ""DONE"", get """//trim(cmdstr)//"""")
    end if

    nrm = sqrt(sum(matmul(transpose(new_d),new_d)))
    write(*,format2) i, new_x, new_f, new_d, nrm

    if (nrm < tol) then
      exit
    end if

    old_x = new_x
    old_f = new_f
    old_d = new_d
    old_h = new_h

    i = i+1
  end do

  write(*,*) "total rosenbrock_f calls: ", rosenbrock_calls(3)

  !--------------------------------------------------------

  write(*,sep)
  write(*,*) my_name//": L-BFGS, memory = 3"
  write(*,header) tol

  vold_x(1) = -1.2_DBL
  vold_x(2) =  1.0_DBL

  call init_lbfgs(3,2)

  i = 1
  do
    if (i == 1) then
      call rosenbrock_f(vold_x,vold_f)
      rosenbrock_calls(4) = rosenbrock_calls(4) + 1
      call rosenbrock_d(vold_x,vold_d)
      write(*,format1) i, vold_x, vold_f, vold_d

      cmdstr = "START"
      call lbfgs_internal(   &
        cmdstr,              &
        vold_x,              &
        vnew_x,              &
        vold_d,              &
        vnew_d,              &
        reset_alpha = .true. &
      )
    else
      cmdstr = "START"
      call lbfgs_internal( &
        cmdstr,            &
        vold_x,            &
        vnew_x,            &
        vold_d,            &
        vnew_d             &
      )
    end if

    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": expected ""EVALUATE_DF1"", get """//trim(cmdstr)//"""")
    end if

    call rosenbrock_f(vnew_x,vnew_f)
    rosenbrock_calls(4) = rosenbrock_calls(4) + 1
    call rosenbrock_d(vnew_x,vnew_d)

    cmdstr = "EVALUATED"
    call lbfgs_internal( &
      cmdstr,            &
      vold_x,            &
      vnew_x,            &
      vold_d,            &
      vnew_d             &
    )

    if (cmdstr /= "DONE") then
      call error(my_name//": expected ""DONE"", get """//trim(cmdstr)//"""")
    end if

    nrm = sqrt(dot_product(vnew_d,vnew_d))
    write(*,format2) i, vnew_x, vnew_f, vnew_d, nrm

    if (nrm < tol) then
      exit
    end if

    vold_x = vnew_x
    vold_f = vnew_f
    vold_d = vnew_d

    i = i+1
  end do

  write(*,*) "total rosenbrock_f calls: ", rosenbrock_calls(4)

  call finalize_lbfgs()

  !--------------------------------------------------------

  write(*,sep)
  write(*,*) my_name//": L-BFGS, memory = 5"
  write(*,header) tol

  vold_x(1) = -1.2_DBL
  vold_x(2) =  1.0_DBL

  call init_lbfgs(5,2)

  i = 1
  do
    if (i == 1) then
      call rosenbrock_f(vold_x,vold_f)
      rosenbrock_calls(5) = rosenbrock_calls(5) + 1
      call rosenbrock_d(vold_x,vold_d)
      write(*,format1) i, vold_x, vold_f, vold_d

      cmdstr = "START"
      call lbfgs_internal(   &
        cmdstr,              &
        vold_x,              &
        vnew_x,              &
        vold_d,              &
        vnew_d,              &
        reset_alpha = .true. &
      )
    else
      cmdstr = "START"
      call lbfgs_internal( &
        cmdstr,            &
        vold_x,            &
        vnew_x,            &
        vold_d,            &
        vnew_d             &
      )
    end if

    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": expected ""EVALUATE_DF1"", get """//trim(cmdstr)//"""")
    end if

    call rosenbrock_f(vnew_x,vnew_f)
    rosenbrock_calls(5) = rosenbrock_calls(5) + 1
    call rosenbrock_d(vnew_x,vnew_d)

    cmdstr = "EVALUATED"
    call lbfgs_internal( &
      cmdstr,            &
      vold_x,            &
      vnew_x,            &
      vold_d,            &
      vnew_d             &
    )

    if (cmdstr /= "DONE") then
      call error(my_name//": expected ""DONE"", get """//trim(cmdstr)//"""")
    end if

    nrm = sqrt(dot_product(vnew_d,vnew_d))
    write(*,format2) i, vnew_x, vnew_f, vnew_d, nrm

    if (nrm < tol) then
      exit
    end if

    vold_x = vnew_x
    vold_f = vnew_f
    vold_d = vnew_d

    i = i+1
  end do

  write(*,*) "total rosenbrock_f calls: ", rosenbrock_calls(5)

  call finalize_lbfgs()

  !--------------------------------------------------------

  write(*,sep)
  write(*,*) my_name//": L-BFGS, memory = 7"
  write(*,header) tol

  vold_x(1) = -1.2_DBL
  vold_x(2) =  1.0_DBL

  call init_lbfgs(7,2)

  i = 1
  do
    if (i == 1) then
      call rosenbrock_f(vold_x,vold_f)
      rosenbrock_calls(6) = rosenbrock_calls(6) + 1
      call rosenbrock_d(vold_x,vold_d)
      write(*,format1) i, vold_x, vold_f, vold_d

      cmdstr = "START"
      call lbfgs_internal(   &
        cmdstr,              &
        vold_x,              &
        vnew_x,              &
        vold_d,              &
        vnew_d,              &
        reset_alpha = .true. &
      )
    else
      cmdstr = "START"
      call lbfgs_internal( &
        cmdstr,            &
        vold_x,            &
        vnew_x,            &
        vold_d,            &
        vnew_d             &
      )
    end if

    if (cmdstr /= "EVALUATE_DF1") then
      call error(my_name//": expected ""EVALUATE_DF1"", get """//trim(cmdstr)//"""")
    end if

    call rosenbrock_f(vnew_x,vnew_f)
    rosenbrock_calls(6) = rosenbrock_calls(6) + 1
    call rosenbrock_d(vnew_x,vnew_d)

    cmdstr = "EVALUATED"
    call lbfgs_internal( &
      cmdstr,            &
      vold_x,            &
      vnew_x,            &
      vold_d,            &
      vnew_d             &
    )

    if (cmdstr /= "DONE") then
      call error(my_name//": expected ""DONE"", get """//trim(cmdstr)//"""")
    end if

    nrm = sqrt(dot_product(vnew_d,vnew_d))
    write(*,format2) i, vnew_x, vnew_f, vnew_d, nrm

    if (nrm < tol) then
      exit
    end if

    vold_x = vnew_x
    vold_f = vnew_f
    vold_d = vnew_d

    i = i+1
  end do

  write(*,*) "total rosenbrock_f calls: ", rosenbrock_calls(6)

  call finalize_lbfgs()

  !--------------------------------------------------------

  write(*,sep)
  write(*,*) "Results, total rosenbrock_f calls:"
  write(*,*) "    Fixed alpha                          : ", rosenbrock_calls(1)
  write(*,*) "    Accelerated backtracking line search : ", rosenbrock_calls(2)
  write(*,*) "    Max displacement                     : ", rosenbrock_calls(3)
  write(*,*) "    L-BFGS (mem = 3)                     : ", rosenbrock_calls(4)
  write(*,*) "    L-BFGS (mem = 5)                     : ", rosenbrock_calls(5)
  write(*,*) "    L-BFGS (mem = 7)                     : ", rosenbrock_calls(6)

  write(*,sep)

end subroutine driver_bfgs_rosenbrock

!====================================================================

subroutine rosenbrock_f(x,f)

  real(DBL), dimension(2), intent(IN)  :: x
  real(DBL),               intent(OUT) :: f

  f = 100.0_DBL*(x(2) - x(1)**2)**2 + (1.0_DBL - x(1))**2

end subroutine rosenbrock_f

!====================================================================

subroutine rosenbrock_d(x,d)

  real(DBL), dimension(2), intent(IN)  :: x
  real(DBL), dimension(2), intent(OUT) :: d

  d(1) = 400.0_DBL*x(1)**3 - 400.0_DBL*x(1)*x(2) - 2.0_DBL + 2.0_DBL*x(1)
  d(2) = 200.0_DBL*x(2) - 200.0_DBL*x(1)**2

end subroutine rosenbrock_d

!====================================================================

end program driver

