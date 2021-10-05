module radiation_mod

! ==================================================================================
! ==================================================================================
  use params, only : grav, sb, cpair,dp, pi

!==================================================================================
implicit none

! public interfaces

public :: radiation_init, radiation_end, radiation_interface
!==================================================================================


! module variables
! Initalisation
logical :: initialized =.false.

! Flag variables
logical :: do_normal_integration_method = .True.
logical :: do_diffusion_approx = .False.
logical :: constant_bond_albedo = .True.
logical :: sw_scat = .False.
logical :: lw_scat = .False.
logical :: kRoss_scheme_V = .False.
logical :: kRoss_scheme_IR = .False.

! Namelist constants
integer :: nV_b            = 1      ! Number of visual bands
integer :: nIR_b           = 1      ! Number of IR bands
real(dp)    :: solar_constant  = 1360.0 ! Solar constant
real(dp)    :: ir_tau_eq       = 256    ! IR equatorial ref optical depth
real(dp)    :: ir_tau_pole     = 0.0    ! IR polar ref optical depth
real(dp)    :: atm_abs         = 12    ! V band ref optical depth
real(dp)    :: linear_tau      = 1.0    ! Linear lw tau dependence
real(dp)    :: albedo_bond     = 0.0    ! Bond albedo
real(dp)    :: albedo_land     = 0.0    ! Land Bond albedo
real(dp)    :: solar_exponent  = 1.0    ! Optical opacity exponent
real(dp)    :: f1              = 1.0    ! Heng et al. (2011) - p**2 dependence exponent
real(dp)    :: reference_slp   = 1.0e5   ! Reference surface pressure
logical :: tidally_locked  = .false.  ! Tidally locked flag
real(dp)    :: noon_longitude  = 0.0      ! Noon longitude flag

! Flux limited diffusion paramaters
real(dp)    :: k_ir            = 1e-3     ! Semi-grey IR  opacity (for Flux Diffusion)
real(dp), parameter    :: E_lim           = 0.01
real(dp), parameter    :: E_fac           = E_lim**(0.023)
real(dp), parameter    :: E_fac1          = 1.0 - E_fac

! Work variables
real(dp) :: solar
real(dp), dimension(:), allocatable :: sw_down, sw_up
real(dp), dimension(:), allocatable :: solar_tau
real(dp), dimension(:), allocatable :: b
real(dp) :: b_surf
real(dp), dimension(:), allocatable :: lw_down, lw_up
real(dp), dimension(:), allocatable :: lw_tau, dtransV, dtransIR
real(dp), dimension(:), allocatable :: lw_net, sw_net, tot_net
real(dp), dimension(:), allocatable :: t_dt_rad, t_dt_sw, t_dt_lw
real(dp)  :: olr, olr_sw

! Math constants
!real(dp), save :: pi, deg_to_rad , rad_to_deg


namelist/radiation_nml/ do_normal_integration_method, do_diffusion_approx, constant_bond_albedo, sw_scat, lw_scat, &
  & solar_constant, ir_tau_eq, ir_tau_pole, atm_abs, linear_tau, albedo_bond, albedo_land, &
  & solar_exponent, f1, reference_slp, tidally_locked, noon_longitude, &
  & nV_b, nIR_b, kRoss_scheme_V, kRoss_scheme_IR


contains


! ==================================================================================
! ==================================================================================

subroutine radiation_interface(p_half, p_full, t, ts, net_flux, &
  kV_R, kIR_R,  Beta_V, Beta_IR, A_Bond)
  implicit none

  real(dp), intent(in), dimension(:) :: t
  real(dp), intent(in) :: ts
  real(dp), intent(in), dimension(:) :: p_half, p_full
  real(dp), dimension(:), intent(out) :: net_flux
  real(dp), optional, dimension(:,:), intent(in) :: kV_R, kIR_R
  real(dp), optional,  dimension(:), intent(in) :: Beta_V 
  real(dp), optional, dimension(:), intent(in) :: Beta_IR
  real(dp), optional, dimension(:), intent(in) :: A_Bond


  
  integer :: i, j, k, n
  integer :: ie, je, nlev
  real(dp), allocatable, dimension(:) :: band_Ab 
  real(dp), allocatable, dimension(:) :: Beta_Vi, Beta_IRi

  integer :: idx_lim

  real(dp), allocatable, dimension(:) :: sw_down_b, sw_up_b, lw_down_b, lw_up_b, lw_diff
  real(dp) :: solar_b
  real(dp) :: solar_tau_ref, lw_tau_ref

  logical :: used
  nlev = size(t) 

  !allocate(band_Ab(ie-is+1,je-js+1,3))
  allocate(sw_down_b(nlev+1), sw_up_b(nlev+1), lw_down_b(nlev+1), lw_up_b(nlev+1))
  allocate(lw_diff(nlev+1))
  allocate(Beta_Vi(nV_b),Beta_IRi(nIR_b))
  allocate(band_AB(nV_b))

  ! Interface to the grey gas radiation scheme

  ! If we have multiple V bands get fraction of incident flux from argument
  ! Otherwise 1 band with 100%
  if (present(Beta_V)) then
    Beta_Vi(:) = Beta_V(:)
  else
    Beta_Vi(:) = 1.0
  end if

  
  !! First do short wave radiation
  ! Loop over longitude and latitude regime

  ! Find solar irradiation flux
  solar = solar_constant
  ! Set total fluxes to zero
      sw_down(:) = 0.0
      sw_up(:) = 0.0

      ! If solar irradiation flux is zero, no sw radiation, so set outputs = 0 and cycle
      if (present(A_Bond) .eqv. .True.) then
        band_AB(:) = A_Bond(:)
      else
        band_AB(:) = albedo_bond
      end if

        ! n loop does number of short wave bands
      do n = 1, nV_b

        ! Zero the sw_down and sw_up flux for band addition (kRoss_scheme_V)
        sw_down_b(:) = 0.0

        ! Find solar radiation tau at each layer
        if (kRoss_scheme_V .eqv. .True.) then
          call tau_struct(nlev, p_half(:), kV_R(n,:), solar_tau(:))
        else
          solar_tau_ref = atm_abs
          do k = 1, nlev+1
            solar_tau(k) = solar_tau_ref * (p_half(k)/reference_slp)**solar_exponent
          end do
        end if

         ! find dtrans for visual band
         do k = 1, nlev
           dtransV(k) = exp(-(solar_tau(k+1) - solar_tau(k)))
         end do
 
        ! solar_b = flux in band (F_TOA * fraction in band * (1 - band bond albedo))
        ! For now it's a constant passed as A_Bond
        solar_b = solar* Beta_Vi(n) * (1.0 - band_AB(n))
         
        ! Do short wave radiation calculation
        if (sw_scat .eqv. .True.) then
          !call sw_grey_scat()
        else
          call sw_grey_down(nlev, solar_b, solar_tau(:), sw_down_b(:))
        end if

        sw_down(:) = sw_down(:) + sw_down_b(:) 
        !sw_up(i,j,:) = sw_up(i,j,:) + sw_up_b(:)

        !print*, '~~~~~'
        !print*, n,i,j, sw_down(i,j,:), '@@'
        !print*, n,i,j, solar_tau(i,j,:)
        !print*, '----'
      end do
 
  !stop
 
  !! Second, do long wave radiation

      lw_down(:) = 0.0
      lw_up(:) = 0.0


      ! If we have multiple IR bands get fraction of blackbody flux from arugment
      ! Otherwise 1 band with 100%
      if (present(Beta_IR)) then
         Beta_IRi(:) = Beta_IR(:)
      else
        Beta_IRi(:) = 1.0
      end if


      do n = 1, nIR_b

        !Planck function at each layer with picket fence fraction
        b(:) = sb * t(:)**4 * Beta_IRi(n)
        !Planck function at surface with picket fence fraction
        b_surf = sb * ts**4 * Beta_IRi(n)

        ! Zero the lw_down and lw_up flux for band addition (kRoss_scheme_IR)
        lw_down_b(:) = 0.0
        lw_up_b(:) = 0.0

        ! Find IR radiation tau at each layer
        if (kRoss_scheme_IR .eqv. .True.) then
          call tau_struct(nlev, p_half(:), kIR_R(n,:), lw_tau(:))
        else
          ! Heng et al. (2011) p^2 dependence expression here. For no p^2 set f1 = 1.0
          lw_tau_ref = ir_tau_eq
          do k = 1, nlev+1
            lw_tau(k) = lw_tau_ref * (p_half(k)/reference_slp) * f1 &
              & + lw_tau_ref * (p_half(k)/reference_slp)**2 * (1.0 - f1)
          end do
        end if

         ! find dtrans for IR band
         do k = 1, nlev
           dtransIR(k) = exp(-(lw_tau(k+1) - lw_tau(k)))
         end do

        ! Do short wave radiation calculation
        if (lw_scat .eqv. .True.) then
          !call lw_grey_down_scat()
          !call lw_grey_up_scat()
        else
          call lw_grey_down(nlev, b(:), lw_tau(:), dtransIR(:), lw_down_b(:))
          call lw_grey_up(nlev, b(:), b_surf, lw_tau(:), dtransIR(:), lw_up_b(:))
        end if

        if (do_diffusion_approx .eqv. .True.) then
          call flux_lim_diffusion(nlev, p_half(:), p_full(:), t(:), &
            & lw_tau(:), lw_diff(:), idx_lim)
          ! Perform flux adjustment below idx_lim
          do k = 1, nlev+1
            if (k < idx_lim) then
              ! Normal flux calculation
              lw_down(k) = lw_down(k) + lw_down_b(k)
              lw_up(k) = lw_up(k) + lw_up_b(k)
            else
              if (lw_diff(k) > 0.0) then
                ! Add to upward flux - down as normal
                lw_down(k) = lw_down(k) + lw_down_b(k)
                lw_up(k) = lw_up(k) + (E_fac * lw_up_b(k) + E_fac1 * lw_diff(k))
              else if (lw_diff(k) < 0.0) then
                ! Add to downward flux - up as normal
                lw_down(k) = lw_down(k) + (E_fac* lw_down_b(k) + E_fac1 * lw_diff(k))
                lw_up(k) = lw_up(k) + lw_up_b(k)
              else
                ! lw_diff = 0.0 - normal flux calculation
                lw_down(k) = lw_down(k) + lw_down_b(k)
                lw_up(k) = lw_up(k) + lw_up_b(k)
              end if
            end if
          end do
        else
          ! Normal flux calculation
          lw_down(:) = lw_down(:) + lw_down_b(:)
          lw_up(:) = lw_up(:) + lw_up_b(:)
        end if  

        !print*, '~~~~~'
        !print*, i,j, lw_down(i,j,:), '@@'
        !print*, i,j, lw_tau(i,j,:)
        !print*, '----'

     end do
     

  !stop

  !! Third, calculate net fluxes
  do k = 1, nlev+1
    lw_net(k) = lw_up(k) - lw_down(k)
    sw_net(k) = sw_up(k) - sw_down(k)
    tot_net(k) = lw_net(k) + sw_net(k)
 end do

 net_flux = tot_net

  !! Fourth, calculate heating rates from radiation
!  do k = 1, nlev
!    ! Tendency from shortwave flux
!    t_dt_sw(:,:,k) = (sw_net(:,:,k+1) - sw_net(:,:,k)) &
!      & * grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))
!    ! Tendency from longwave flux
!    t_dt_lw(:,:,k) = (lw_net(:,:,k+1) - lw_net(:,:,k)) &
!      & * grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))
    ! Net tendency
!    t_dt_rad(:,:,k) = t_dt_sw(:,:,k) + t_dt_lw(:,:,k)

    !t_dt_rad(:,:,k) = (lw_net(:,:,k+1) - lw_net(:,:,k) - sw_down(:,:,k+1) + sw_down(:,:,k)) &
     ! & * grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))



  !do i = is, ie
    !do j = js, je
      !print*, i,j,lat(i,j),lon(i,j),t_dt_rad(i,j,:)
    !end do
  !end do
  !stop

  olr= lw_up(1)   ! lw OLR
  olr_sw = sw_up(1) ! sw OLR (from scattering)


end subroutine radiation_interface

subroutine sw_grey_down(nlev, solar, solar_tau, sw_down)
  implicit none

  integer, intent(in) :: nlev
  real(dp), intent(in) :: solar
  real(dp), dimension(nlev+1), intent(in) :: solar_tau
  real(dp), dimension(nlev+1), intent(out) :: sw_down

  integer :: k

  do k = 1, nlev+1
    sw_down(k) = solar * exp(-solar_tau(k)) 
  end do
  
end subroutine sw_grey_down

subroutine lw_grey_down(nlev, b, lw_tau, dtrans, lw_down_b)
  implicit none

  integer, intent(in) :: nlev
  real(dp), dimension(nlev+1), intent(in) :: b, lw_tau
  real(dp), dimension(nlev), intent(in) :: dtrans
  real(dp), dimension(nlev+1), intent(out) :: lw_down_b

  integer :: k

  ! boundary condition - zero flux down at TOA
  lw_down_b(1) = 0.0

  if (do_normal_integration_method .eqv. .True.) then
    ! Normal exponential dependence
    do k = 1, nlev
      lw_down_b(k+1) = lw_down_b(k)*dtrans(k) + b(k)*(1.0-dtrans(k))
    end do
  else
    ! 1st order high optical depth fix from ISCA
    do k = 1, nlev
      lw_down_b(k+1) = 2.0 * b(k) * (lw_tau(k+1) - lw_tau(k)) / (2.0 + (lw_tau(k+1) - lw_tau(k))) + &
        & lw_down_b(k) * (2.0 - (lw_tau(k+1) - lw_tau(k)))/(2.0 + (lw_tau(k+1) - lw_tau(k)))
    end do
  end if 

end subroutine lw_grey_down

subroutine lw_grey_up(nlev, b, b_surf, lw_tau, dtrans, lw_up_b)
  implicit none

  integer, intent(in) :: nlev
  real(dp), intent(in) :: b_surf
  real(dp), dimension(nlev+1), intent(in) :: b, lw_tau
  real(dp), dimension(nlev), intent(in) ::  dtrans
  real(dp), dimension(nlev+1), intent(out) :: lw_up_b

  integer :: k

  ! Boundary condition - upward flux from surface
  lw_up_b(nlev+1) = b_surf
  
  if (do_normal_integration_method .eqv. .True.) then
    ! Normal exponential dependence
    do k = nlev, 1, -1
      lw_up_b(k) = lw_up_b(k+1)*dtrans(k) + b(k)*(1.0-dtrans(k))
    end do
  else
    ! 1st order high optical depth fix from ISCA
    do k = nlev, 1,-1
      lw_up_b(k) = 2.0 * b(k) * -1.0 * (lw_tau(k) - lw_tau(k+1)) / (2.0 - (lw_tau(k) - lw_tau(k+1))) + &
        & lw_up_b(k+1) * (2.0 + (lw_tau(k) - lw_tau(k+1)))/(2.0 - (lw_tau(k) - lw_tau(k+1)))
    end do
  end if

end subroutine lw_grey_up

subroutine flux_lim_diffusion(nlev, p_edge, p_full, t_full, lw_tau, F_diff, idx_lim, kIR_R)
  implicit none

  integer, intent(in) :: nlev
  real(dp), dimension(nlev+1), intent(in) :: p_edge
  real(dp), dimension(nlev), intent(in) :: p_full
  real(dp), dimension(nlev), intent(in) :: t_full
  real(dp), dimension(nlev+1), intent(in) :: lw_tau
  real(dp), dimension(nlev+1), intent(out) :: F_diff
  real(dp), dimension(nlev), intent(in), optional :: kIR_R
  integer, intent(out) :: idx_lim

  ! Work Vairiables
  integer :: k
  real(dp) :: tau_lim, A, dt_full, dp_full, dkIR_full
  real(dp) :: t_edge
  real(dp) :: kIR_r_edge

  ! We follow the Rauscher & Menou (2012) scheme
  ! First find the optical depth limit for the diffusion
  A = log10(maxval(p_edge)/minval(p_edge)) / dble(nlev + 1)
  tau_lim = 5.55 / (10.0**(A) - 10.0**(-A))

  ! Perform flux calculation
  do k = 1, nlev
    ! If above optical depth limit, set F_diff = 0 and keep track of idx of limit
    if (lw_tau(k) < tau_lim) then
      F_diff(k) = 0.0
      idx_lim = k+1
      cycle
    end if
    ! If in diffusion regime, perform diffusion calculation
    ! First calculate temperature at cell edge by interpolation
    dt_full = t_full(k-1) - t_full(k)
    dp_full = p_full(k-1) - p_full(k)
    t_edge = t_full(k-1) + (p_edge(k) - p_full(k-1)) * dt_full/dp_full

    if (present(kIR_R) .eqv. .True.) then
      ! Interpolate kIR to edge values
      dKIR_full = kIR_R(k-1) - kIR_R(k)
      kIR_r_edge = kIR_R(k-1) + (p_edge(k) - p_full(k-1)) * dkIR_full/dp_full
    else
      ! use a constant value
      kIR_r_edge = k_ir
    endif

    F_diff(k) = ((16.0 * grav * sb * t_edge**3) / (3.0 * kIR_r_edge)) &
      &  * dt_full/dp_full  

  end do
  F_diff(nlev+1) = 0.0

end subroutine flux_lim_diffusion


subroutine tau_struct(nlev, p_half, kRoss, tau_struc)
  implicit none

  integer, intent(in) :: nlev
  real(dp), dimension(nlev+1), intent(in) :: p_half
  real(dp), dimension(nlev), intent(in) :: kRoss !, ssa, asy_g
  real(dp), dimension(nlev+1), intent(out) :: tau_struc

  real(dp) :: tau_sum, tau_lay, dP
  integer :: k
  
  ! running sum of optical depth
  tau_sum = 0.0

  ! Upper most tau_struc is given by some low pressure value (here 1e-9 bar = 1e-4 pa)
  !dP = (p_half(1) - 1e-4)
  !tau_lay = (kRoss(1) * dP) / grav
  !tau_sum = tau_sum + tau_lay
  tau_struc(1) = tau_sum

  ! Integrate from top to bottom
  do k = 1, nlev
 
    ! Pressure difference between layer edges
    dP = (p_half(k+1) - p_half(k))

    ! Optical depth of layer assuming hydrostatic equilibirum
    tau_lay = (kRoss(k) * dP) / grav
 
    ! Add to running sum
    tau_sum = tau_sum + tau_lay
  
    ! Optical depth structure is running sum
    tau_struc(k+1) = tau_sum

  end do

end subroutine tau_struct

! ==================================================================================

subroutine radiation_init(num_levels)

!-------------------------------------------------------------------------------------
integer, intent(in)               ::  num_levels
!-------------------------------------------------------------------------------------
integer :: ierr, io, unit, i

logical :: water_file_exists

!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile
unit  = 15
open(unit, file='input.nml')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=radiation_nml, iostat=io, end=10)
enddo
10 close(unit)

!pi    = 4.0 * atan(1.0)
!deg_to_rad = 2.0*pi/360.0
!rad_to_deg = 360.0/2.0/pi

initialized = .true.

! Allocate all work variables
allocate(sw_down(num_levels+1))
allocate(sw_up(num_levels+1))
allocate(solar_tau(num_levels+1))
allocate(dtransV(num_levels))

allocate(b(num_levels))

allocate(lw_down(num_levels+1))
allocate(lw_up(num_levels+1))
allocate(lw_tau(num_levels+1))
allocate(dtransIR(num_levels))

allocate(lw_net(num_levels+1))
allocate(sw_net(num_levels+1))
allocate(tot_net(num_levels+1))


end subroutine radiation_init


! ==================================================================================


subroutine radiation_end

deallocate(sw_down, sw_up)
deallocate(solar_tau)
deallocate(b)
deallocate(lw_down, lw_up)
deallocate(lw_tau, dtransV, dtransIR)
deallocate(lw_net,sw_net, tot_net)

end subroutine radiation_end



! ==================================================================================

end module radiation_mod

