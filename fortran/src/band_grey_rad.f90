module band_grey_mod

  use netcdf
  use iso_fortran_env
  use io,             only: handle_err
  use utils,          only: bilinear_interp, linear_log_interp, linear_interp, lin_interp
  use params,         only: dp, nf, Finc, Fint, sb, grav, opacity_dir, invert_grid, rdgas, sw_fac, &
                           lw_fac
  use twostr_module,  only: twostr
  use nc_wrap_mod
  use phys_mod,       only: i_h2, i_h2o, i_ch4, &
                            i_nh3, i_co2, i_co, u_vec
  
  implicit none

  integer :: nB
  !integer :: nB, nT, nP   ! Array data
  real(dp), allocatable, dimension(:,:,:) :: h2o_kap, nh3_kap, ch4_kap, co_kap, co2_kap
       
  real(dp), allocatable, dimension(:,:)::  h2h2_cia, h2oh2o_cia
  
  real(dp), allocatable, dimension(:) :: pref_h2o, tref_h2o, &
                                         pref_ch4, tref_ch4, pref_co2,  tref_co2,  &
                                         pref_nh3, tref_nh3, pref_co,   tref_co,   &
                                          tref_h2oh2o, tref_h2h2
  real(dp), allocatable, dimension(:,:) :: bands
  real(dp), allocatable, dimension(:)   :: fbeam
  real(dp), allocatable, dimension(:,:) :: dtau_c, dtau_u, ctau_banded, utau_banded
  
contains

  subroutine read_in_solar(file_name)
    character(len=*), intent(in) :: file_name
    integer :: ncid, status, id

    call open_file(file_name, ncid)
    call get_var(ncid, 'solar_frac', fbeam)
    call close_file(ncid)

    fbeam = Finc*fbeam
  end subroutine read_in_solar

  subroutine read_in_bands(file_name)
    character(*), intent(in) :: file_name

    integer :: ncid
    call open_file(file_name, ncid)

    call get_dimension_length(ncid, 'band', nB)
    call get_var(ncid, 'bands', bands)

    call close_file(ncid)
    
  end subroutine read_in_bands

  subroutine bg_init(npz)
    integer, intent(in) :: npz

    integer :: ncid
    character(80) :: h2o_file, nh3_file, ch4_file, co_file, co2_file, &
                     h2h2_file, h2oh2o_file, solar_file

    real(dp), allocatable :: pref_temp(:)
    
    h2o_file=trim(opacity_dir)//'/'//'h2o_opacity.nc'
    nh3_file=trim(opacity_dir)//'/'//'nh3_opacity.nc'
    co2_file=trim(opacity_dir)//'/'//'co2_opacity.nc'
    !co_file=trim(opacity_dir)//'/'//'co_opacity.nc'
    ch4_file=trim(opacity_dir)//'/'//'ch4_opacity.nc'

    h2h2_file = trim(opacity_dir)//'/'//'h2h2_cia.nc'
    h2oh2o_file = trim(opacity_dir)//'/'//'h2oh2o_cia.nc'
    solar_file = trim(opacity_dir)//'/'//'frac.nc'
    
    ! Read in band data
    call read_in_bands(h2o_file)

    ! Read in solar data
    call read_in_solar(solar_file)

    ! Read in species data
    call get_kap_data(h2o_file, tref_h2o, pref_h2o, h2o_kap)
    call get_kap_data(h2o_file, tref_ch4, pref_ch4, ch4_kap)
    call get_kap_data(co2_file, tref_co2, pref_co2, co2_kap)
    call get_kap_data(nh3_file, tref_nh3, pref_nh3, nh3_kap)
    !call get_kap_data(co_file, tref_co, pref_co, co_kap)

    call get_cia_data(h2h2_file, tref_h2h2, h2h2_cia)
    !call get_cia_data(h2oh2o_file, tref_h2oh2o, h2oh2o_cia)

    ! Allocate tau shape
    allocate(dtau_u(nB, npz))

    if (invert_grid) then
       allocate(dtau_c(nB, npz+1))
    else
       allocate(dtau_c(nB, npz))
    endif

  end subroutine bg_init

  subroutine get_kap_data(file, tref, pref, opacity)
    character(len=*), intent(in) :: file
    real(dp), allocatable, intent(inout) :: pref(:), tref(:), opacity(:,:,:)
    
    integer :: ncid, b
    real(dp), allocatable :: ross_opacity(:,:,:), solar_opacity(:,:,:)
    
    call open_file(file, ncid)
    call get_var(ncid, "pres", pref)
    call get_var(ncid, "temp", tref)
    call get_var(ncid, "opacity", opacity)
    !call get_var(ncid, "ross_opacity", ross_opacity)
    call get_var(ncid, "solar_opacity", solar_opacity)

    opacity(1,:,:) = solar_opacity(1,:,:) ! Make 1st band solar weighted
!!$    do b=5,5
!!$       write(*,*) opacity(5,:,:)/ross_opacity(5,:,:)
!!$       opacity(b,:,:) = ross_opacity(b,:,:)
!!$     enddo
    
    
    call close_file(ncid)
    !deallocate(ross_opacity)
    deallocate(solar_opacity)
  end subroutine get_kap_data

  subroutine get_cia_data(file, tref, opacity)
    character(len=*), intent(in) :: file
    real(dp), allocatable, intent(inout) :: tref(:), opacity(:,:)

    integer :: ncid

    call open_file(file, ncid)
    call get_var(ncid, "temp", tref)
    call get_var(ncid, "opacity", opacity)
    call close_file(ncid)

  end subroutine get_cia_data

  subroutine calc_tau(P, Pe, T, Te, delP, qwat)
    real(dp), dimension(:), intent(in) :: T,P, delP, qwat, Te, Pe

    real(dp) :: kappa(size(T)), k1(size(T)), k2(size(T))
    real(dp) :: kp1,kp2,kp_tot, p1, p2, T1, T2, q1, q2, qe(size(Te))
    real(dp) :: kp(7)
    integer :: npz, k, b, test(1), test2(1), m


    real(dp) :: R, mu, denom, mu_d, alpha(size(T))
    real(dp), parameter :: Rstar = 8314.5
    real(dp) :: xs(5) ! Dry molar mixing ratios (dry mol/total dry mol) - these are constant
    real(dp) :: q_spec(6) ! order:
    ! H2, CH4, NH3, CO2, CO,  H2O
    npz  = size(T)

    denom = 1. 
    !xs = (/1., 3.7e-2, 1.e-2, 0., 0./)
    !xs = (/1., 3.7e-2, 8.28e-4, 8.84e-4, 4.17e-3/)
    !xs = xs/denom
    !xs(i_h2) = 1. - sum(xs(2:))

    !mu_d = sum(u_vec(1:5)*xs)
    !alpha = mu_d/u_vec(i_h2o) * qwat/(1-qwat)
    !alpha= 0.0
    !mu = (mu_d + alpha*u_vec(i_h2o))/(1 + alpha)
    
    !mu = 18.*0.001 + 16*0.001 +2*(1-0.002)
    mu = 18.* 5.67e-2 + 16.*3.70e-2 + 17*8.28e-4 + 44*8.84e-4 + 4.17e-3*28 &
         + (1 - 5.67e-2 - 3.70e-2 - 8.28e-4 - 8.84e-4 - 4.17e-3)*2.

    !mu = 18.*0.001 + 16.*0.001 + (1-0.002)*2
!    R = Rstar/mu

!!$    do k=1,5
!!$       q_spec(k, :) = u_vec(k)*xs(k)/(mu_d + alpha*u_vec(i_h2o))
!!$    enddo
!!$    q_spec(i_h2o,:) = qwat(:)
    
    q_spec(2) = 16./mu * 0.001!3.70e-2
    q_spec(3) = 17./mu * 8.28e-4
    q_spec(4) = 44./mu * 8.84e-4
    q_spec(5) = 28./mu * 4.17e-3
    q_spec(6) = 18./mu * 5.67e-2

    q_spec(1) = 1.
    do k=2,6
       q_spec(1) = q_spec(1) - q_spec(k)
    enddo

    R = rstar/mu

    !write(*,*) 'R=', R, q_spec(1), q_spec(2), q_spec(3)
    dtau_c = 0._dp
    !q_spec(1) = 0.9
    !q_spec(2) = 0.1
    !q_spec(3:6) = 0.0
    
    if(invert_grid) then
       p1 =(pe(1) + p(1))/2
       call linear_interp(p1, Pe(1), P(1), Te(1), T(1), T1)
       p2 = (p(nf) + pe(nf+1))/2
       call linear_interp(p2, P(nf), Pe(nf+1), T(nf), Te(nf+1), T2)

       ! Interpolate water to these layers
!!$       call linear_interp(p1, P(1), P(2), qwat(1), qwat(2), q1)
!!$       call linear_interp(p2, P(nf-1), P(nf), qwat(nf-1), qwat(nf), q2)
!!$
!!$       ! Interpolate water to interfaces
!!$       do k=2,npz
!!$          call linear_log_interp(pe(k), P(k-1), P(k), T(k-1), T(k), qe(k))
!!$       enddo
!!$       call linear_log_interp(pe(1), P(1), P(2), T(1), T(2), qe(1))
!!$       call linear_log_interp(pe(npz+1), P(npz-1), P(npz), T(npz-1), T(npz), qe(npz+1))
       kp(5) = 0.
       do b=1,nB
          
          call bilinear_interp( tref_ch4, log(pref_ch4), ch4_kap(:,:,b), T1, log(p1), kp(2))
          call bilinear_interp( tref_nh3, log(pref_nh3), nh3_kap(:,:,b), T1, log(p1), kp(3))
          call bilinear_interp( tref_co2, log(pref_co2), co2_kap(:,:,b), T1, log(p1), kp(4))
          !call bilinear_interp( tref_co, log(pref_co), co_kap(:,:,b), T1, log(p1), kp(5))
          call bilinear_interp( tref_h2o, log(pref_h2o), h2o_kap(:,:,b), T1, log(p1), kp(6))
          
          call lin_interp(tref_h2h2, h2h2_cia(:,b), T1, kp(1))
          !call lin_interp(tref_h2oh2o, h2oh2o_cia(:,b), T1, kp(7))
          

          dtau_c(b,1) = 0.
          do k=2,6
             dtau_c(b,1) = dtau_c(b,1) + delP(1)/grav * kp(k)*q_spec(k)
          enddo
          dtau_c(b,1) = dtau_c(b,1) + delp(1)/grav * (q_spec(1)**2.*kp(1))&! + q_spec(6)**2.*kp(7)) &
                                                   * p1/R/T1
          !dtau_c(b,1) = (kp1*q1 + kp2*(1-q1)) * delP(1)/grav
          
          !call bilinear_interp( tref_h2o, pref_h2o, h2o_kap(:,:,b), T2, p2, kp1)
          !call bilinear_interp( tref_h2h2, pref_h2h2, h2h2_cia(:,:,b), T2, p2, kp2)

          call bilinear_interp( tref_ch4, log(pref_ch4), ch4_kap(:,:,b), T2, log(p2), kp(2))
          call bilinear_interp( tref_nh3, log(pref_nh3), nh3_kap(:,:,b), T2, log(p2), kp(3))
          call bilinear_interp( tref_co2, log(pref_co2), co2_kap(:,:,b), T2, log(p2), kp(4))
          !call bilinear_interp( tref_co, log(pref_co), co_kap(:,:,b), T2, log(p2), kp(5))
          call bilinear_interp( tref_h2o, log(pref_h2o), h2o_kap(:,:,b), T2, log(p2), kp(6))

          call lin_interp(tref_h2h2, h2h2_cia(:,b), T2, kp(1))
          !call lin_interp(tref_h2oh2o, h2oh2o_cia(:,b), T2, kp(7))
          
          ! Weight CIA coefficient with density
          !kp2 = kp2*p2/T2/rdgas

          dtau_c(b,nf+1) = 0.
          do m=2,6
             dtau_c(b,nf+1) = dtau_c(b,nf+1) + delP(nf+1)/grav * kp(m)*q_spec(m)
          enddo
          dtau_c(b,nf+1) = dtau_c(b,nf+1) + delp(nf+1)/grav * (q_spec(1)**2.*kp(1))&! + q_spec(6)**2.*kp(7))&
               * p2/R/T2

          !dtau_c(b,nf+1) = (kp1*q2 + kp2*(1-q2)) * delp(nf+1)/grav
          
          do k=2,npz

             call bilinear_interp( tref_ch4, log(pref_ch4), ch4_kap(:,:,b), Te(k), log(Pe(k)), kp(2))
             call bilinear_interp( tref_nh3, log(pref_nh3), nh3_kap(:,:,b), Te(k), log(Pe(k)), kp(3))
             call bilinear_interp( tref_co2, log(pref_co2), co2_kap(:,:,b), Te(k), log(Pe(k)), kp(4))
             !call bilinear_interp( tref_co, log(pref_co), co_kap(:,:,b), Te(k), log(Pe(k)), kp(5))
             call bilinear_interp( tref_h2o, log(pref_h2o), h2o_kap(:,:,b), Te(k), log(Pe(k)), kp(6))
             
             call lin_interp(tref_h2h2, h2h2_cia(:,b), Te(k), kp(1))
             !call lin_interp(tref_h2oh2o, h2oh2o_cia(:,b), Te(k), kp(7))

             !call bilinear_interp(tref_h2o, pref_h2o, h2o_kap(:,:,b), Te(k),  Pe(k), kp1)
             !call bilinear_interp(tref_h2o, pref_h2h2, h2h2_cia(:,:,b), Te(k),  Pe(k), kp2)

             ! Weight CIA coefficient by density
             !kp2 = kp2*Pe(k)/rdgas/Te(k)
             dtau_c(b,k) = 0.
             do m=2,6
                dtau_c(b,k) = dtau_c(b,k) + delP(k)/grav * kp(m)*q_spec(m)
             enddo
             dtau_c(b,k) = dtau_c(b,k) + delp(k)/grav * (q_spec(1)**2.*kp(1))&! + q_spec(6)**2.*kp(7))&
                                                       * Pe(k)/R/Te(k)

             if (isnan(dtau_c(b,k)))then
                write(*,*) 'b = ', b, ', k = ', k, ', is NAN'
                write(*,*) 'KP: ', kp
                write(*,*) 'TEMP, PRES: ', Te(k), Pe(k)
                stop
             endif

!!$             write(*,*) 'KAPPAS, BAND ', b, 'LEVEL ', k

!!$             write(*,*) sum(dtau_c(b,1:k))
!!$             write(*,*) q_spec(1)**2.*kp(1)*Pe(k)/R/Te(k), q_spec(2)*kp(2), q_spec(3)*kp(3), q_spec(4)*kp(4)
!!$             write(*,*) q_spec(5)*kp(5), q_spec(6)*kp(6), q_spec(2)**2.*kp(7)*Pe(k)/R/Te(k)
             !dtau_c(b,k) = (kp1*qwat(k) + kp2*(1-qwat(k))) * delp(k)/grav
             if (b .eq. 2) then
                dtau_c(b,k) = dtau_c(b,k)/lw_fac
             else if (b .eq. 1) then
                dtau_c(b,k) = dtau_c(b,k)/sw_fac
             endif
          enddo
       enddo
       allocate(ctau_banded(nB, npz+2))
    endif
    
    allocate(utau_banded(nB, npz+1))
    dtau_u = 0.0_dp
    do b=1, nB
       !writE(*,*) 'Band: ', b
       !write(*,*) '------------------------------------------------------------------------'
       do  k=1,npz
          
          !call bilinear_interp( tref_h2o, pref_h2o, h2o_kap(:,:,b), T(k), P(k), kp1)
          !call bilinear_interp( tref_h2h2, pref_h2h2, h2h2_cia(:,:,b), T(k), P(k), kp2)
!         print*, shape(tref_h2o), shape(pref_h2o),'next',  shape(h2o_kap)
          
          call bilinear_interp( tref_ch4, log(pref_ch4), ch4_kap(:,:,b), T(k), log(P(k)), kp(2))
          call bilinear_interp( tref_nh3, log(pref_nh3), nh3_kap(:,:,b), T(k), log(P(k)), kp(3))
          call bilinear_interp( tref_co2, log(pref_co2), co2_kap(:,:,b), T(k), log(P(k)), kp(4))
!          call bilinear_interp( tref_co, log(pref_co), co_kap(:,:,b), T(k), log(P(k)), kp(5))
          call bilinear_interp( tref_h2o, log(pref_h2o), h2o_kap(:,:,b), T(k), log(P(k)), kp(6))
!!$
         call lin_interp(tref_h2h2, h2h2_cia(:,b), T(k), kp(1))
         !call lin_interp(tref_h2oh2o, h2oh2o_cia(:,b), T(k), kp(7))

          !kp2 = kp2*P(k)/rdgas/T(k)
          ! Is adding the opacities legit?? Ask Ray, understand random overlap approx
          ! in correlated-k approximation
!          write(*,*) kp1, kp2

!!$          if (kp1 .gt. maxval(h2o_kap)) then
!!$             write(*,*) 'kp1 > max(h2o_kap)'
!!$             write(*,*) 'T, P'
!!$             write(*,*) T(k), P(k)
!!$             test = minloc(abs(T(k) - tref_h2o))
!!$             test2= minloc(abs(P(k) - pref_h2o))
!!$             write(*,*) 'NEAREST NEIGHBOUR VS KP1'
!!$             write(*,*) h2o_kap(test(1), test2(1), b), kp1
!!$          endif
         !          write(*,*) kp1*qwat(k)/kp2/(1-qwat(k))

         dtau_u(b,k) = 0.
             do m=2,6
                dtau_u(b,k) = dtau_u(b,k) + delP(k)/grav * kp(m)*q_spec(m)
             enddo
             
             dtau_u(b,k) = dtau_u(b,k) + delp(k)/grav * (q_spec(1)**2.*kp(1))&! + q_spec(6)**2.*kp(7)) &
                                                       * P(k)/R/T(k)
             
            !dtau_u(b,k) = (kp1*qwat(k) + kp2*(1-qwat(k)))*delP(k)/grav
          if (dtau_u(b,k) .lt. 0) then
             write(*,*) 'SUSPICIOUS DTAU', Te(k), Pe(k), kp
             stop
          endif

          if (( b .gt. 1) ) then
             dtau_u(b,k) = dtau_u(b,k)/lw_fac!/10000
          else if (b .eq. 1) then
                dtau_u(b,k) = dtau_u(b,k)/sw_fac
          endif
          
!          write(*,*) b, maxval(kp(2:6)*q_spec(2:6,k)), maxloc(kp(2:6)*q_spec(2:6,k)), &
!               P(k)/R(k)/T(k)*kp(1)*q_spec(1,k)**2., P(k)/R(k)/T(k)*kp(7)*q_spec(2,k)**2. , q_spec(3,k), &
!               q_spec(3,k)*kp(3)
       enddo
       
       
    enddo

!!$    do k=1,npz
!!$       write(*,*) dtau_u(1,k),dtau_u(2,k),dtau_u(3,k),dtau_u(4,k),dtau_u(5,k)
!!$       enddo
    
    if ( .not. invert_grid ) then
       dtau_c = dtau_u
    else
       utau_banded= 0
       ctau_banded = 0.
       do k=1,npz+1
          ctau_banded(:,k+1) = ctau_banded(:,k) + dtau_c(:,k)
       enddo

       do b=1,nB
          do k=1,npz-1
             call linear_interp(P(k), Pe(k), Pe(k+1), ctau_banded(b,k+1), ctau_banded(b,k+2), &
                  utau_banded(b,k+1))
          enddo
          utau_banded(b,npz+1) = ctau_banded(b, npz+2)
       enddo
       
    endif

!!$    do k=1,npz
!!$       write(*,*) ctau_banded(4, k), utau_banded(4,k)
!!$    enddo
    
  end subroutine calc_tau

  subroutine run_twostr(npz,Tf, Te, P, Pe, delP, qwat, flux, olr, tau_V, tau_IR)
    integer, intent(in) :: npz
    real(dp), dimension(:), intent(in)    :: Tf, Te, P, delP, qwat, tau_V, tau_IR, Pe
    real(dp), dimension(:), intent(out) :: flux
    real(dp), intent(out) :: olr

    integer :: k, n
    integer, parameter :: nerr=22
    logical, allocatable :: plnk(:)
    
    ! TWOSTR parameters

    ! Computational layer structure
    real(dp), allocatable :: ssalb(:) 
    real(dp), allocatable :: gg(:) 
    real(dp) :: temp(0:npz) !Level temperature
    real(dp) :: temp2(0:npz+1)
    real(dp), allocatable :: zd(:)
    integer :: maxcly
    integer :: maxulv
    integer :: ntau
!    real(dp), allocatable :: utau_banded(:,:)
    
    
    ! Top/bottom boundary conditions
    real(dp) :: albedo = 0. ! Xianyu hack for internal temp
    real(dp) :: btemp = 1000
    real(dp) :: um0 = 1. ! Cosine of incident beam
    real(dp) :: fisot = 0._dp
    real(dp) :: radius
    real(dp) :: ttemp = 0.0_dp
    real(dp) :: temis = 0.0_dp
    real(dp) :: wvnmlow, wvnmhigh
    
    ! Control variables
    logical :: deltam = .true.
    logical :: planck
    logical :: quiet = .true.
    logical :: spher = .false.
    logical :: usrtau
    logical :: prnt(2) = (/.false.,.false./)
    character(len=127) :: header
    integer :: ierror(nerr)

    ! Output variables
    real(dp), dimension(npz+1) :: rfldir, rfldn, flup, dfdt, uavg, rfldir_test

    ! HII defined variables
    real(dp), dimension(npz+1) :: rfldir_tot, rfldn_tot, flup_tot, tau1, tau2
    
    integer :: ierr, i
    integer:: j(1)

    
    ! Hack, think of better way of doing this?
    allocate(plnk(nB))
    !    allocate(utau_banded(nB, npz+1))
    plnk = .false.
    do i=1,nB
       plnk(i) = .true.
    enddo
    
    btemp = (Fint/sb)**0.25
    !write(*,*) btemp
    ! BTEMP CHANGED
    
!    btemp = 1000 
    ! Initialise dtau
    !write(*,*) qwat
    call calc_tau(P, Pe, Tf, Te, delP, qwat)


    ntau = npz+1
    
    ! Get temperature in correct form for twostr, starting at 0
    if (.not. invert_grid) then
       maxcly = npz
       usrtau = .false.
       do k=0,npz
          temp(k) = Te(k+1)
       enddo
    else
       maxcly = npz+1
       usrtau = .true.
       temp2(0) = Te(1)
       temp2(npz+1) = Te(npz+1)
       do k=1,npz
          temp2(k) = Tf(k)
       enddo
       utau_banded = 0.
    
       do k=1,npz
          utau_banded(:,k+1) = utau_banded(:,k) + dtau_u(:,k)
       enddo
    endif
    
    allocate(ssalb(maxcly), gg(maxcly), zd(0:maxcly))
    ssalb = 0.0_dp
    gg = 0.0_dp

    
    rfldir_tot = 0.0_dp
    rfldn_tot = 0.0_dp
    flup_tot = 0.0_dp

    maxulv = npz + 1
    
    
    ! With dtau set, iterate over bands

!!$    ! HACK - CHANGE DTAU TO CONSTANT VALUE
!!$    !---------------------------------------
!!$    do n=1,npz
!!$       dtau(1,n) = 2*delp(n)/delp(npz)
!!$       dtau(2,n) = 2*delp(n)/delp(npz)
!!$    enddo
!!$    !--------------------------------------
!    write(*,*) 'SHAPES'
    !    write(*,*) shape(bands)
    !write(*,*) utau_banded(3v,:)
    do n=1,nB
       !write(*,*) 'BAND', n, 'DTAUC', dtau_c(n,:)
       wvnmlow = 1./bands(2,n)/100.
       wvnmhigh = 1./bands(1,n)/100.
!       print*, 'BAND:1', dtau_c(n,:)
       if (.not. invert_grid) then

          call twostr(albedo, temp(npz), deltam, dtau_c(n,:), fbeam(n), fisot, gg, header, ierror, &
               maxcly, maxulv, npz, plnk(n), ntau, prnt, quiet, radius, spher, &
               ssalb, temis, temp, ttemp, um0, usrtau, utau_banded(n,:), wvnmlow, wvnmhigh, &
               zd, dfdt, flup, rfldir, rfldn, uavg)
       else
!!$          write(*,*) '--------------------------------'
!!$          do k=1,npz
!!$             write(*,*) dtau_c(n,k), utau_banded(n,k)
!!$         enddo
!!$          write(*,*) '--------------------------------'

          !write(*,*) albedo, temp2(npz+1), fbeam(n)
          call twostr(albedo, temp2(npz+1), deltam, dtau_c(n,:), fbeam(n), fisot, gg, header, ierror, &
               maxcly, maxulv, npz+1, plnk(n), ntau, prnt, quiet, radius, spher, &
               ssalb, temis, temp2, ttemp, um0, usrtau, utau_banded(n,:), wvnmlow, wvnmhigh, &
               zd, dfdt, flup, rfldir, rfldn, uavg)
!!$          write(*,*) 'band, rfldir(bottom), utau(bottom), fbeam , flup(bottom)'
!!$          write(*,*) n, rfldir(maxulv), utau_banded(n,maxulv), fbeam(n), flup(maxulv)
       endif
          do ierr = 1, nerr
             if ( ierror(ierr) .ne. 0 ) then
                WRITE(*,'(/,A,I4,/)')  "TWOSTR REPORTS FATAL ERROR: ", &
                     ierr
             endif
          enddo

       rfldir_tot = rfldir_tot + rfldir
       rfldn_tot = rfldn_tot + rfldn
       flup_tot  = flup_tot + flup
!!$       if (n .eq. 1) then
!!$          write(*,*)'RFLDIR: ', rfldir(npz+1), 'utau_banded: ', utau_banded(1,npz+1)
!!$       endif
       !if (n .eq. 1) rfldir_test = rfldir
    enddo

    ! Extract output net flux
    flux = flup_tot - rfldn_tot - rfldir_tot
    !do k=1,npz+1
    !   write(*,*) flux(k), flup_tot(k), rfldn_tot(k), rfldir_tot(k)
    !enddo
    
    olr = flup_tot(1)
    write(*,*) 'SW TAU at BOA', sum(dtau_u(1,:))
    write(*,*) 'LW TAU at BOA', sum(dtau_u(2,:))
    open(9,file='tau2.out')
    do i=1,nB
       write(9,'(ES13.4)') ( dtau_u(i,k) ,k=1,npz)
    end do
    close(9)

!!$   write(*,*) '-----------------'
!!$   write(*,*) 'flup_tot:, ', (flup_tot(n), n=1,npz+1)
!!$   write(*,*) '-----------------'
!!$    write(*,*) 'rfldn_tot:, ', (rfldn_tot(n), n=1,npz+1)
!!$    write(*,*) '-------------------'
!!$   write(*,*) 'rfldir_tot: ', rfldir_tot
!!$   write(*,*) '-----------------------'
    
    !write(*,*) 'USER DTAU: ', (utau(n+1)-utau(n), n=1,npz)
!!$    write(*,*) '-----------------'
!!$    write(*,*) 'COMP DTAU :, '
!!$    do k=1,npz
!!$      write(*,*) sum(dtau_u(1,:k)), sum(dtau_u(2,:k))!, dtau_c(3,k), dtau_c(4,k), dtau_c(5,k)
!!$    enddo
!!$    
!!$    write(*,*) '-----------------'
!!$
!!$    write    (*,*) 'flux :, ', (flux(n) ,n=1,npz+1)
!!$    write(*,*) '-----------------'
!!$    write(*,*) 'Differenced flux: ', (flux(n+1) - flux(n), n=1,npz)
!!$    write(*,*) '-----------------'
!!$    write(*,*) 'MAXVAL H2OKAP, H2H2CIA'
!!$    write(*,*) maxval(h2o_kap), maxval(h2h2_cia)
    deallocate(plnk, ssalb, gg, zd, utau_banded)
    if (allocated(ctau_banded) ) then
       deallocate(ctau_banded)
    end if
    
  end subroutine run_twostr
end module band_grey_mod
