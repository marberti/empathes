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

module elastic

  !==================================================================
  !   Elastic Module
  !==================================================================
  ! It deals with elastic forces
  ! and total forces (elastic + potential).
  !==================================================================

  use utility
  use geometry
  use idpp
  use pes_data
  use pes

  implicit none
  save
  private

  ! protected variables -----------------------------------
  public    :: total_forces,                &
               parall_elastic_forces,       &
               perpen_pes_forces,           &
               norm_tangent
  protected :: total_forces,                &
               parall_elastic_forces,       &
               perpen_pes_forces,           &
               norm_tangent
  ! pubic procedures --------------------------------------
  public    :: set_spring_k,                &
               set_spring_mode,             &
               set_mpe_mode,                &
               init_elastic_module,         &
               compute_total_forces,        &
               arbitrary_geom_total_forces, &
               set_total_forces_on_i

  !--------------------------------------------------------
  ! ENUM
  integer, parameter                     :: STATIC_SPRING            = 0
  integer, parameter                     :: DYNAMIC_SPRING           = 1
  integer, parameter                     :: HYBRID_SPRING            = 2
  integer                                :: spring_mode              = DYNAMIC_SPRING

  ! ENUM
  ! to set the minimum path estimator (mpe) mode
  integer, parameter                     :: NEB_MPE                  = 0
  integer, parameter                     :: EB_MPE                   = 1
  integer                                :: mpe_mode                 = NEB_MPE

  logical                                :: flag_spring_k            = .false.
  logical                                :: flag_init_elastic_module = .false.
  real(DBL), allocatable, dimension(:)   :: dynamic_spring_k
  real(DBL), allocatable, dimension(:,:) :: total_forces
  real(DBL), allocatable, dimension(:,:) :: parall_elastic_forces
  real(DBL), allocatable, dimension(:,:) :: perpen_pes_forces
  real(DBL), allocatable, dimension(:,:) :: perpen_idpp_forces
  real(DBL), allocatable, dimension(:,:) :: tangent
  real(DBL), allocatable, dimension(:,:) :: norm_tangent
  real(DBL), allocatable, dimension(:,:) :: half_tangent
  real(DBL)                              :: spring_k

contains

!====================================================================
! Public
!====================================================================

subroutine set_spring_k(k)

  real(DBL), intent(IN) :: k

  if (flag_spring_k) then
    call error("set_spring_k: spring constant already setted")
  end if

  spring_k = k

  flag_spring_k = .true.

end subroutine set_spring_k

!====================================================================

subroutine set_spring_mode(str)

  character(*), intent(INOUT) :: str

  logical, save               :: first_call = .true.

  if (first_call.eqv..false.) then
    call error("set_spring_mode: subroutine called more than once")
  end if

  call tolower(str)

  select case (str)
  case ("static")
    spring_mode=STATIC_SPRING
  case ("dynamic")
    spring_mode=DYNAMIC_SPRING
  case ("hybrid")
    spring_mode=HYBRID_SPRING
    !TODO remove when it's ready
    call error("set_spring_mode: sorry, hybrid spring mode is not implemented yet")
  case default
    call error("set_spring_mode: unknown spring mode """//&
      &trim(str)//"""")
  end select

  first_call = .false.

end subroutine set_spring_mode

!====================================================================

subroutine set_mpe_mode(str)

  character(*), intent(INOUT) :: str

  character(*), parameter     :: my_name    = "set_mpe_mode"
  logical, save               :: first_call = .true.

  if (first_call.eqv..false.) then
    call error(my_name//": subroutine called more than once")
  end if

  call tolower(str)

  select case (str)
  case ("neb")
    mpe_mode = NEB_MPE
  case ("eb")
    mpe_mode = EB_MPE
  case default
    call error(my_name//": unknown MPE mode """//trim(str)//"""")
  end select

  first_call = .false.

end subroutine set_mpe_mode

!====================================================================

subroutine init_elastic_module()

  !--------------------------------------------------------
  ! init_elastic_module performs all the preliminary checks
  ! required to compute elastic forces.
  ! For this reason, all the other procedures inside this module
  ! don't need to check if their arguments are valid.
  !--------------------------------------------------------

  integer        :: err_n
  character(120) :: err_msg

  ! preliminary checks ------------------------------------
  if (flag_init_images.eqv..false.) then
    call error("init_elastic_module: images not initialized")
  end if

  if (flag_init_elastic_module) then
    call error("init_elastic_module: module elastic already initialized")
  end if

  ! allocation section ------------------------------------
  allocate(total_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_elastic_module: "//trim(err_msg))
  end if

  allocate(parall_elastic_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_elastic_module: "//trim(err_msg))
  end if

  allocate(perpen_pes_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_elastic_module: "//trim(err_msg))
  end if

  allocate(perpen_idpp_forces(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_elastic_module: "//trim(err_msg))
  end if

  allocate(tangent(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_elastic_module: "//trim(err_msg))
  end if

  allocate(norm_tangent(image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_elastic_module: "//trim(err_msg))
  end if

  allocate(half_tangent(0:image_n,geom_len),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error("init_elastic_module: "//trim(err_msg))
  end if

  select case (spring_mode)
  case (DYNAMIC_SPRING,HYBRID_SPRING)
    allocate(dynamic_spring_k(image_n+1),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("init_elastic_module: "//trim(err_msg))
    end if
  end select

  ! variables setting -------------------------------------
  if (flag_spring_k.eqv..false.) then
    spring_k      = 1.0_DBL
    flag_spring_k = .true.
  end if

  flag_init_elastic_module = .true.

end subroutine init_elastic_module

!====================================================================

subroutine compute_total_forces(mode,fixed)

  integer, intent(IN) :: mode
  logical, intent(IN) :: fixed

  integer             :: atoms
  character(8)        :: istr

  if (flag_init_elastic_module.eqv..false.) then
    call error("compute_total_forces: module elastic not initialized")
  end if

  select case (mode)
  case (PES_MODE)
    call compute_pes_forces()
    call init_tangents(mode, mpe_mode)
    call compute_parall_elastic_forces(mpe_mode)
    call compute_perpen_pes_forces()
    total_forces = parall_elastic_forces+perpen_pes_forces
  case (IDPP_MODE)
    call compute_idpp_enfo()
    call init_tangents(mode, NEB_MPE)
    call compute_parall_elastic_forces(NEB_MPE)
    call compute_perpen_idpp_forces()
    total_forces = parall_elastic_forces+perpen_idpp_forces
  case default
    write(istr,'(I8)') mode
    istr = adjustl(istr)
    call error("compute_total_forces: mode """//trim(istr)//""" not valid")
  end select

  if (fixed) then
    atoms=geom_len/3
    if (atoms>=1) total_forces(:,1:3) = 0.0_DBL
    if (atoms>=2) total_forces(:,5:6) = 0.0_DBL
    if (atoms>=3) total_forces(:,9:9) = 0.0_DBL
  end if

end subroutine compute_total_forces

!====================================================================

subroutine arbitrary_geom_total_forces(i,ig,f,g)

  integer,                       intent(IN)  :: i      ! image number
  real(DBL), dimension(:),       intent(IN)  :: ig     ! image geometry (Ang)
  real(DBL),                     intent(OUT) :: f      ! pes energy in point ig (Hartree)
  real(DBL), dimension(:),       intent(OUT) :: g      ! total forces in point ig (Hartree/Ang)

  real(DBL), save, allocatable, dimension(:) :: pesg   ! pes forces (Hartree/Ang)
  real(DBL), save, allocatable, dimension(:) :: taup   ! t+
  real(DBL), save, allocatable, dimension(:) :: taum   ! t-
  real(DBL), save, allocatable, dimension(:) :: tg     ! tangent
  real(DBL), save, allocatable, dimension(:) :: ntg    ! normalized tangent
  real(DBL), save, allocatable, dimension(:) :: paelfo ! parallel elastic forces
  real(DBL), save, allocatable, dimension(:) :: pepefo ! perpendicular pes forces
  real(DBL)                                  :: prev
  real(DBL)                                  :: curr
  real(DBL)                                  :: next
  real(DBL)                                  :: c1
  real(DBL)                                  :: c2
  real(DBL)                                  :: coeff
  real(DBL)                                  :: e_max
  real(DBL)                                  :: e_min
  real(DBL)                                  :: dp
  real(DBL)                                  :: conv_threshold
  logical                                    :: converged
  integer                                    :: err_n
  character(120)                             :: err_msg

  if (.not.flag_init_elastic_module) then
    call error("arbitrary_geom_total_forces: module elastic not initialized")
  end if

  ! argument bounds check ---------------------------------
  if ((i<1).or.(i>image_n)) then
    write(err_msg,'(I8)') i
    call error("arbitrary_geom_total_forces: image i out of bounds: "&
      &//trim(adjustl(err_msg)))
  end if

  if (size(ig,1)/=geom_len) then
    call error("arbitrary_geom_total_forces: wrong ig argument size")
  end if

  ! allocation section ------------------------------------
  if (.not.allocated(pesg)) then
    allocate(pesg(geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("arbitrary_geom_total_forces: "//trim(err_msg))
    end if
  end if

  if (.not.allocated(taup)) then
    allocate(taup(geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("arbitrary_geom_total_forces: "//trim(err_msg))
    end if
  end if

  if (.not.allocated(taum)) then
    allocate(taum(geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("arbitrary_geom_total_forces: "//trim(err_msg))
    end if
  end if

  if (.not.allocated(tg)) then
    allocate(tg(geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("arbitrary_geom_total_forces: "//trim(err_msg))
    end if
  end if

  if (.not.allocated(ntg)) then
    allocate(ntg(geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("arbitrary_geom_total_forces: "//trim(err_msg))
    end if
  end if

  if (.not.allocated(paelfo)) then
    allocate(paelfo(geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("arbitrary_geom_total_forces: "//trim(err_msg))
    end if
  end if

  if (.not.allocated(pepefo)) then
    allocate(pepefo(geom_len),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("arbitrary_geom_total_forces: "//trim(err_msg))
    end if
  end if

  ! working section ---------------------------------------
  conv_threshold = get_scfconv()
  call get_pes_forces(i,0,conv_threshold,converged,ig,f,pesg)

  !WARNING Bad Programming Style --------------------------
  ! This code is the same as in:
  !   * compute_half_tangent
  !   * compute_tangent
  !   * normalize_tangent
  !   * compute_parall_elastic_forces
  !   * compute_perpen_pes_forces
  ! but for a single image point.
  ! If a single change is needed, then it must be applied
  ! in different parts of the code.
  !--------------------------------------------------------

  ! compute half-tangents ---------------------------------
  taup = image_geom(i+1,:)-ig
  taum = ig-image_geom(i-1,:)

  ! compute tangent ---------------------------------------
  prev = pes_energy(i-1)
  curr = f
  next = pes_energy(i+1)

  if (prev==next) then
    tg = taup+taum
  else if ((next>curr).and.(curr>prev)) then
    tg = taup
  else if ((next<curr).and.(curr<prev)) then
    tg = taum
  else
    e_max = max(abs(next-curr),abs(prev-curr))
    e_min = min(abs(next-curr),abs(prev-curr))
    if (next>prev) then
      tg = (taup*e_max)+(taum*e_min)
    else
      tg = (taup*e_min)+(taum*e_max)
    end if
  end if

  ! normalize tangent -------------------------------------
  ntg = tg/norm(tg)

  ! compute parallel elastic forces -----------------------
  c1     = norm(taup)
  c2     = norm(taum)
  coeff  = c1-c2
  paelfo = (spring_k*coeff)*ntg

  ! compute perpendicular pes forces ----------------------
  dp     = dot_product(pesg,ntg)
  pepefo = pesg-(dp*ntg)

  ! compute total forces ----------------------------------
  g = paelfo+pepefo

end subroutine arbitrary_geom_total_forces

!====================================================================

subroutine set_total_forces_on_i(i,forces)

  integer,                 intent(IN) :: i
  real(DBL), dimension(:), intent(IN) :: forces

  if ((i<1).or.(i>image_n)) then
    call error("set_total_forces_on_i: image i out of bounds")
  end if

  if (size(forces,1)/=geom_len) then
    call error("set_total_forces_on_i: wrong forces argument size")
  end if

  total_forces(i,:) = forces

end subroutine set_total_forces_on_i

!====================================================================
! Private
!====================================================================

subroutine init_tangents(mode,mpe_mode_arg)

  !--------------------------------------------------------
  ! WARNING: init_tangents *must* be called
  !          after compute_pes_forces or compute_idpp_enfo.
  !          That's because energies are needed
  !          for tangents extimation.
  !--------------------------------------------------------

  integer, intent(IN) :: mode
  integer, intent(IN) :: mpe_mode_arg

  if (flag_init_elastic_module.eqv..false.) then
    call error("init_tangents: module elastic not initialized")
  end if

  call compute_half_tangent()
  call compute_tangent(mode, mpe_mode_arg)
  call normalize_tangent()

end subroutine init_tangents

!====================================================================

subroutine compute_dynamic_spring_k()

  !--------------------------------------------------------
  ! henkelman2000jcp doi: 10.1063/1.1329672
  ! Computes spring force constants dynamically.
  ! Must be used only in PES_MODE.
  !--------------------------------------------------------

  real(DBL), allocatable, dimension(:) :: e
  real(DBL)                            :: e_max
  real(DBL)                            :: e_ref
  real(DBL)                            :: k_max
  real(DBL)                            :: dk
  integer                              :: i
  integer                              :: err_n
  character(120)                       :: err_msg

  ! preliminary checks ------------------------------------
  select case (spring_mode)
  case (DYNAMIC_SPRING,HYBRID_SPRING)
    ! the only cases in which this subroutine call is valid.
    ! do nothing.
  case (STATIC_SPRING)
    call error("compute_dynamic_spring_k: illegal call with static spring mode")
  case default
    call error("compute_dynamic_spring_k: unknown spring mode")
  end select

  ! allocation section ------------------------------------
  if (.not.allocated(e)) then
    allocate(e(image_n+1),stat=err_n,errmsg=err_msg)
    if (err_n/=0) then
      call error("compute_dynamic_spring_k: """//trim(err_msg)//"""")
    end if
  end if

  ! working section ---------------------------------------
  do i=1, image_n+1
    e(i) = max(pes_energy(i),pes_energy(i-1))
  end do

  e_max = maxval(e,1)
  e_ref = max(pes_energy(0),pes_energy(image_n+1))
  ! this two are programmer defined parameter
  k_max = spring_k
  dk    = k_max/2.0_DBL

  do i=1, image_n+1
    if (e(i)>e_ref) then
      dynamic_spring_k(i) = k_max - dk*((e_max-e(i))/(e_max-e_ref))
    else
      dynamic_spring_k(i) = k_max - dk
    end if
  end do

end subroutine compute_dynamic_spring_k

!====================================================================

subroutine compute_parall_elastic_forces(mpe_mode_arg)

  integer,     intent(IN) :: mpe_mode_arg

  character(*), parameter :: my_name = "compute_parall_elastic_forces"
  character(8)            :: istr

  select case (mpe_mode_arg)
  case (NEB_MPE)
    call compute_parall_elastic_forces_neb()
  case (EB_MPE)
    call compute_parall_elastic_forces_eb()
  case default
    write (istr,'(I8)') mpe_mode_arg
    istr = adjustl(istr)
    call error(my_name//": mpe_mode_arg """//trim(istr)//""" not valid")
  end select

end subroutine compute_parall_elastic_forces

!====================================================================

subroutine compute_parall_elastic_forces_neb()

  integer   :: i
  real(DBL) :: coeff
  real(DBL) :: c1
  real(DBL) :: c2

  select case (spring_mode)
  case (STATIC_SPRING)
    do i=1, image_n
      c1    = norm(half_tangent(i,:))
      c2    = norm(half_tangent(i-1,:))
      coeff = c1-c2
      parall_elastic_forces(i,:) = (spring_k*coeff)*norm_tangent(i,:)
    end do
  case (DYNAMIC_SPRING,HYBRID_SPRING)
    call compute_dynamic_spring_k()
    do i=1, image_n
      c1    = norm(half_tangent(i,:))
      c2    = norm(half_tangent(i-1,:))
      coeff = dynamic_spring_k(i+1)*c1 - dynamic_spring_k(i)*c2
      parall_elastic_forces(i,:) = coeff*norm_tangent(i,:)
    end do
  case default
    call error("compute_parall_elastic_forces: unknown spring mode")
  end select

end subroutine compute_parall_elastic_forces_neb

!====================================================================

subroutine compute_parall_elastic_forces_eb()

  ! Implementation of parallel elastic forces
  ! as described in the Elastic Band method (EB).
  ! Equations in appendix C.
  !
  ! source: kolsbjerg2016automated; doi: 10.1063/1.4961868

  character(*),              parameter :: my_name = "compute_parall_elastic_forces_eb"
  real(DBL)                            :: l_eq
  real(DBL), dimension(:), allocatable :: c1
  real(DBL), dimension(:), allocatable :: c2
  real(DBL)                            :: vmax
  real(DBL)                            :: vmin
  integer                              :: i
  logical, dimension(:), allocatable   :: neighbors
  integer                              :: err_n
  character(120)                       :: err_msg

  ! allocation section ------------------------------------
  allocate(c1(geom_len), stat=err_n, errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(c2(geom_len), stat=err_n, errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  allocate(neighbors(0:image_n+1), stat=err_n, errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! compute l_eq ------------------------------------------
  l_eq = norm(image_geom(image_n+1,:) - image_geom(0,:)) / (image_n + 1)

  ! locate the neighbors of the energy maximum ------------
  call get_maxima_neighbors(pes_energy,neighbors,1)

  ! compute parallel elastic forces -----------------------
  do i=1, image_n
    c1 = (norm(half_tangent(i,:))   - l_eq) * (half_tangent(i,:)/norm(half_tangent(i,:)))
    c2 = (norm(half_tangent(i-1,:)) - l_eq) * (half_tangent(i-1,:)/norm(half_tangent(i-1,:)))

    select case (spring_mode)
    case (STATIC_SPRING)
      parall_elastic_forces(i,:) = spring_k * (c1 - c2)
    case (DYNAMIC_SPRING,HYBRID_SPRING)
      call compute_dynamic_spring_k()
      parall_elastic_forces(i,:) = dynamic_spring_k(i+1)*c1 - dynamic_spring_k(i)*c2
    case default
      call error(my_name//": unknown spring mode")
    end select

    if (neighbors(i).eqv..true.) then
      vmax = max(abs(pes_energy(i+1)-pes_energy(i)), abs(pes_energy(i-1)-pes_energy(i)))
      vmin = min(abs(pes_energy(i+1)-pes_energy(i)), abs(pes_energy(i-1)-pes_energy(i)))

      parall_elastic_forces(i,:) = parall_elastic_forces(i,:) * (vmin/vmax)
    end if
  end do

  ! deallocation section ----------------------------------
  deallocate(c1, stat=err_n, errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(c2, stat=err_n, errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

  deallocate(neighbors, stat=err_n, errmsg=err_msg)
  if (err_n /= 0) then
    call error(my_name//": "//trim(err_msg))
  end if

end subroutine compute_parall_elastic_forces_eb

!====================================================================

subroutine compute_perpen_pes_forces()

  integer   :: i
  real(DBL) :: dp

  do i=1, image_n
    dp = dot_product(pes_forces(i,:),norm_tangent(i,:))
    perpen_pes_forces(i,:) = pes_forces(i,:)-(dp*norm_tangent(i,:))
  end do

end subroutine compute_perpen_pes_forces

!====================================================================

subroutine compute_perpen_idpp_forces()

  integer   :: i
  real(DBL) :: dp

  do i=1, image_n
    dp = dot_product(idpp_forces(i,:),norm_tangent(i,:))
    perpen_idpp_forces(i,:) = idpp_forces(i,:)-(dp*norm_tangent(i,:))
  end do

end subroutine compute_perpen_idpp_forces

!====================================================================

subroutine compute_half_tangent()

  !--------------------------------------------------------
  ! Do not call directly, use init_tangents instead
  !--------------------------------------------------------

  integer :: i

  do i=0, image_n
    half_tangent(i,:) = image_geom(i+1,:)-image_geom(i,:)
  end do

end subroutine compute_half_tangent

!====================================================================

subroutine compute_tangent(mode,mpe_mode_arg)

  !--------------------------------------------------------
  ! Do not call directly, use init_tangents instead
  !--------------------------------------------------------

  integer,                  intent(IN) :: mode
  integer,                  intent(IN) :: mpe_mode_arg

  character(*),              parameter :: my_name = "compute_tangent"
  integer                              :: i
  real(DBL)                            :: prev
  real(DBL)                            :: curr
  real(DBL)                            :: next
  real(DBL)                            :: e_max
  real(DBL)                            :: e_min
  real(DBL), allocatable, dimension(:) :: energy
  character(8)                         :: istr
  integer                              :: err_n
  character(120)                       :: err_msg

  ! allocation section ------------------------------------
  allocate(energy(0:image_n+1),stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

  ! mode selection ----------------------------------------
  select case (mode)
  case (PES_MODE)
    energy = pes_energy
  case (IDPP_MODE)
    if (mpe_mode_arg /= NEB_MPE) then
      call error(my_name//": IDPP_MODE used without NEB_MPE")
    end if
    energy = idpp_energy
  case default
    write(istr,'(I8)') mode
    istr = adjustl(istr)
    call error(my_name//": mode """//trim(istr)//""" not valid")
  end select

  ! working section ---------------------------------------
  select case (mpe_mode_arg)
  case (NEB_MPE)

    do i=1, image_n
      prev = energy(i-1)
      curr = energy(i)
      next = energy(i+1)

      if (prev==next) then
        tangent(i,:) = half_tangent(i,:)+half_tangent(i-1,:)
      else if ((next>curr).and.(curr>prev)) then
        tangent(i,:) = half_tangent(i,:)
      else if ((next<curr).and.(curr<prev)) then
        tangent(i,:) = half_tangent(i-1,:)
      else
        e_max = max(abs(next-curr),abs(prev-curr))
        e_min = min(abs(next-curr),abs(prev-curr))
        if (next>prev) then
          tangent(i,:) = (half_tangent(i,:)*e_max)+(half_tangent(i-1,:)*e_min)
        else
          tangent(i,:) = (half_tangent(i,:)*e_min)+(half_tangent(i-1,:)*e_max)
        end if
      end if
    end do

  case (EB_MPE)

    do i=1, image_n
      tangent(i,:) = (half_tangent(i-1,:)/norm(half_tangent(i-1,:))) + &
        &(half_tangent(i,:)/norm(half_tangent(i,:)))
    end do

  case default
    write(istr,'(I8)') mpe_mode_arg
    istr = adjustl(istr)
    call error(my_name//": mpe_mode_arg """//trim(istr)//""" not valid")
  end select

  ! deallocation section ----------------------------------
  deallocate(energy,stat=err_n,errmsg=err_msg)
  if (err_n/=0) then
    call error(my_name//": "//trim(err_msg))
  end if

end subroutine compute_tangent

!====================================================================

subroutine normalize_tangent()

  !--------------------------------------------------------
  ! Do not call directly, use init_tangents instead
  !--------------------------------------------------------

  integer :: i

  do i=1, image_n
    norm_tangent(i,:) = tangent(i,:)/norm(tangent(i,:))
  end do

end subroutine normalize_tangent

!====================================================================

end module elastic

