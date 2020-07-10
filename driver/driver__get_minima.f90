program driver

  use utility

  integer, parameter :: adim = 8
  integer, parameter :: astart=0
  integer, parameter :: aend=adim-1
  real(DBL), dimension(astart:aend) :: a, b
  logical, dimension(astart:aend) :: f
  integer :: i, mn, mx, indx

  b(0)=-10970.500533_DBL
  b(1)=-10970.414585_DBL
  b(2)=-10970.166030_DBL
  b(3)=-10969.738550_DBL
  b(4)=-10969.540730_DBL
  b(5)=-10969.748298_DBL
  b(6)=-10969.942877_DBL
  b(7)=-10970.001704_DBL

  a = b

  write(*,*) "--- Initial Array ---"
  do i=astart, aend
    write(*,*) i,a(i)
  end do

  write(*,*)
  write(*,*) "--- Get Minima ---"
  call get_minima(a,f)
  mn=0
  do i=astart, aend
    write(*,*) i,f(i)
    if (f(i)) mn=mn+1
  end do
  write(*,*) "tot minima", mn
  do i=1,mn
    indx=minloc(a,1,f)+(astart-1)
    write(*,*) indx
    a(indx)=0.0_DBL
    f(indx)=.false.
  end do

  write(*,*)
  write(*,*) "Result"
  do i=astart, aend
    write(*,*) i,a(i)
  end do

  a = b

  write(*,*)
  write(*,*) "--- Get Maxima ---"
  call get_maxima(a,f)
  mx=0
  do i=astart, aend
    write(*,*) i,f(i)
    if (f(i)) mx=mx+1
  end do
  write(*,*) "tot maxima", mx
  do i=1,mx
    indx=maxloc(a,1,f)+(astart-1)
    write(*,*) indx
    a(indx)=0.0_DBL
    f(indx)=.false.
  end do

  write(*,*)
  write(*,*) "Result"
  do i=astart, aend
    write(*,*) i,a(i)
  end do

end program driver

