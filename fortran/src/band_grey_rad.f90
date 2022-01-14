module band_grey_mod

  use netcdf
  use iso_fortran_env
  use io,             only: handle_err
  use utils,          only: bilinear_interp, writeline
  use params,         only: dp, nf, Finc, Fint, sb, grav, opacity_dir
  use twostr_module,  only: twostr
  use nc_wrap_mod
  
  implicit none

  integer :: nB
  !integer :: nB, nT, nP   ! Array data
  real(dp), allocatable, dimension(:,:,:) :: h2o_kap, h2h2_cia
  real(dp), allocatable, dimension(:) :: pref_h2o, pref_h2h2, tref_h2o, tref_h2h2
  real(dp), allocatable, dimension(:,:) :: bands
  real(dp), allocatable, dimension(:)   :: fbeam
  real(dp), allocatable, dimension(:,:) :: dtau
  
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
    character(80) :: water_file, h2_file, solar_file

    water_file=trim(opacity_dir)//'/'//'h2o_opacity.nc'
    h2_file = trim(opacity_dir)//'/'//'h2h2_cia.nc'
    solar_file = trim(opacity_dir)//'/'//'frac.nc'
    
    ! Read in band data
    call read_in_bands(water_file)

    ! Read in solar data
    call read_in_solar(solar_file)

    ! Read in water data
    call open_file(water_file, ncid)
    call get_var(ncid, "pres", pref_h2o)
    call get_var(ncid, "temp", tref_h2o)
    call get_var(ncid, "opacity", h2o_kap)
    call close_file(ncid)

    ! Read in H2-H2 cia data
    call open_file(h2_file, ncid)
    call get_var(ncid, "pres", pref_h2h2)
    call get_var(ncid, "temp", tref_h2h2)
    call get_var(ncid, "opacity", h2h2_cia)
    call close_file(ncid)

    ! Allocate tau shape
    allocate(dtau(nB, npz))
  end subroutine bg_init

  subroutine calc_tau(P, T, delP, qwat)
    real(dp), dimension(:), intent(in) :: T,P, delP, qwat

    real(dp) :: kappa(size(T)), k1(size(T)), k2(size(T))
    real(dp) :: kp1,kp2,kp_tot
    integer :: npz, k, b, test(1), test2(1)
    

    npz  = size(T)
    dtau = 0.0_dp
    do b=1, nB
       !writE(*,*) 'Band: ', b
       !write(*,*) '------------------------------------------------------------------------'
       do  k=1,npz
          
          call bilinear_interp( tref_h2o, pref_h2o, h2o_kap(:,:,b), T(k), P(k), kp1)
          call bilinear_interp( tref_h2h2, pref_h2h2, h2h2_cia(:,:,b), T(k), P(k), kp2)

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
          dtau(b,k) = (kp1*qwat(k) + kp2*(1-qwat(k)))*delP(k)/grav
          
       enddo
    enddo
  end subroutine calc_tau

  subroutine run_twostr(npz,Tf, Te, P, delP, qwat, flux, olr, tau_V, tau_IR)
    integer, intent(in) :: npz
    real(dp), dimension(:), intent(in)    :: Tf, Te, P, delP, qwat, tau_V, tau_IR
    real(dp), dimension(:), intent(out) :: flux
    real(dp), intent(out) :: olr

    integer :: k, nB,n
    integer, parameter :: nerr=22
    logical, allocatable :: plnk(:)
    
    ! TWOSTR parameters

    ! Computational layer structure
    real(dp) :: ssalb(1:npz) 
    real(dp) :: gg(1:npz) 
    real(dp) :: temp(0:npz) !Level temperature
    real(dp) :: zd(0:npz)
    integer :: maxcly
    integer :: maxulv
    integer :: ntau
    real(dp) :: utau(npz+1)
    
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
    logical :: usrtau = .false.
    logical :: prnt(2) = (/.false.,.false./)
    character(len=127) :: header
    integer :: ierror(nerr)

    ! Output variables
    real(dp), dimension(npz+1) :: rfldir, rfldn, flup, dfdt, uavg, rfldir_test

    ! HII defined variables
    real(dp), dimension(npz+1) :: rfldir_tot, rfldn_tot, flup_tot, tau1, tau2
    
    integer :: ierr, i
    integer:: j(1)

    ssalb = 0.0_dp
    gg = 0.0_dp
    
    ! Get temperature in correct form for twostr, starting at 0
    do k=0,npz
       temp(k) = Te(k+1)
    enddo
    
    nB = size(dtau,1)

    ! Hack, think of better way of doing this?
    allocate(plnk(nB))
    plnk = (/.true., .true., .true., .true., .true./)
    btemp = (Fint/sb)**0.25
    ! BTEMP CHANGED
    
!    btemp = 1000 
    ! Initialise dtau
    !write(*,*) qwat
    call calc_tau(P, Tf, delP, qwat)
    
    !write(*,*) 'SHAPE', shape(dtau)
    tau1 = 0.0
    tau2 = 0.0
!!$    do k=1,npz
!!$       dtau(1,k) = tau_V(k+1) - tau_V(k)
!!$       dtau(2,k) = tau_IR(k+1) - tau_IR(k)
!      tau1(k+1) = tau1(k) + dtau(1,k)
!      tau2(k+1) = tau2(k) + dtau(2,k)
!      write(*,*) tau1(k+1), tau2(k+1)
!    enddo
    !write(*,*) tau_V
    !write(*,*) dtau(1,:)
!!$    write(*,*) '---------------'
!!$    write(*,*) (tau(k), k=1,npz+1)
!!$    write(*,*) '---------------'
!!$    write(*,*) 'MINLOC', minloc(abs(tau-1))
!!$    write(*,*) '---------------'

!!$    j = minloc(abs(tau-1))
!!$    i = j(1)
!!$
!!$    write(*,*) '---------------'
!!$    write(*,*) 'PE(i): ', P(i)
!!$    write(*,*) '---------------'
    
    rfldir_tot = 0.0_dp
    rfldn_tot = 0.0_dp
    flup_tot = 0.0_dp

    maxulv = npz + 1
    maxcly = npz 
    
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
    do n=1,nB
       wvnmlow = 1./bands(2,n)/100.
       wvnmhigh = 1./bands(1,n)/100.
       call twostr(albedo, temp(npz), deltam, dtau(n,:), fbeam(n), fisot, gg, header, ierror, &
            maxcly, maxulv, npz, plnk(n), ntau, prnt, quiet, radius, spher, &
            ssalb, temis, temp, ttemp, um0, usrtau, utau, wvnmlow, wvnmhigh, &
            zd, dfdt, flup, rfldir, rfldn, uavg)
       do ierr = 1, nerr
          if ( ierror(ierr) .ne. 0 ) then
             WRITE(*,'(/,A,I4,/)')  "TWOSTR REPORTS FATAL ERROR: ", &
                  ierr
          endif
       enddo

       rfldir_tot = rfldir_tot + rfldir
       rfldn_tot = rfldn_tot + rfldn
       flup_tot  = flup_tot + flup

       !if (n .eq. 1) rfldir_test = rfldir
    enddo

    ! Extract output net flux
    flux = flup_tot - rfldn_tot - rfldir_tot
    olr = flup_tot(1)

!!$    write(*,*) '-----------------'
!!$    write(*,*) 'flup_tot:, ', (flup_tot(n), n=1,npz+1)
!!$    write(*,*) '-----------------'
!!$    write(*,*) 'rfldn_tot:, ', (rfldn_tot(n), n=1,npz+1)
!!$    write(*,*) '-------------------'
!!$    write(*,*) 'rfldir_tot: ', rfldir_tot
!!$    write(*,*) '-----------------------'
    
!!$    write(*,*) 'USER DTAU: ', (utau(n+1)-utau(n), n=1,npz)
!!$    write(*,*) '-----------------'
!!$    write(*,*) 'COMP DTAU :, ', (dtau(2,n) ,n=1,npz)
!!$    write(*,*) '-----------------'
!!$
!!$    write(*,*) 'flux :, ', (flux(n) ,n=1,npz+1)
!!$    write(*,*) '-----------------'
!!$    write(*,*) 'Differenced flux: ', (flux(n+1) - flux(n), n=1,npz)
!!$    write(*,*) '-----------------'
!!$    write(*,*) 'MAXVAL H2OKAP, H2H2CIA'
!!$    write(*,*) maxval(h2o_kap), maxval(h2h2_cia)
    deallocate(plnk)
  end subroutine run_twostr
end module band_grey_mod
