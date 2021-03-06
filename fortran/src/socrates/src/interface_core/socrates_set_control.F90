! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates control type
!
!------------------------------------------------------------------------------
module socrates_set_control
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_CONTROL'
contains

subroutine set_control(control, spectrum, l_set_defaults, &
  l_rayleigh, l_gas, l_continuum, l_cont_gen, l_orog, l_solvar, &
  l_rescale, l_ir_source_quad, l_mixing_ratio, &
  l_aerosol, l_aerosol_mode, l_aerosol_ccn, &
  l_tile, l_clear, &
  l_flux_up_band, l_flux_down_band, &
  l_flux_up_clear_band, l_flux_down_clear_band, &
  l_blue_flux_surf, &
  n_tile, n_cloud_layer, n_aer_mode, &
  isolir, i_cloud_representation, i_overlap, i_inhom, i_mcica_sampling, &
  i_st_water, i_cnv_water, i_st_ice, i_cnv_ice )

use def_control,  only: StrCtrl, allocate_control
use def_spectrum, only: StrSpecData
use realtype_rd,  only: RealK
use ereport_mod,  only: ereport
use errormessagelength_mod, only: errormessagelength
use missing_data_mod, only: imdi
use rad_pcf, only: &
  ip_solar, ip_pifm80, ip_scatter_full, ip_infra_red, ip_elsasser, &
  ip_two_stream, ip_ir_gauss, ip_spherical_harmonic, ip_overlap_k_eqv_scl, &
  ip_cloud_off, ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat, &
  ip_cloud_csiw, ip_max_rand, ip_rand, ip_exponential_rand, ip_homogeneous, &
  ip_scaling, ip_mcica, ip_cairns, ip_cloud_ice_water, ip_cloud_mcica, &
  ip_no_scatter_abs, ip_no_scatter_ext, ip_solver_no_scat, &
  ip_solver_homogen_direct, ip_scatter_approx, ip_solver_mix_app_scat, &
  ip_solver_homogen_direct, ip_solver_mix_direct_hogan, ip_cloud_mix_max, &
  ip_cloud_part_corr, ip_cloud_mix_random, ip_solver_triple_app_scat, &
  ip_solver_triple_hogan, ip_cloud_triple, ip_cloud_part_corr_cnv, &
  ip_cloud_clear, ip_scale_ses2, ip_overlap_mix_ses2, &
  i_normal, i_err_fatal
use def_mcica, only: ip_mcica_optimal_sampling

implicit none

! Control options:
type(StrCtrl), intent(inout) :: control

! Spectral data:
type(StrSpecData), intent(in), optional :: spectrum

logical, intent(in), optional :: l_set_defaults, &
  l_rayleigh, l_gas, l_continuum, l_cont_gen, l_orog, l_solvar, &
  l_rescale, l_ir_source_quad, l_mixing_ratio, &
  l_aerosol, l_aerosol_mode, l_aerosol_ccn, &
  l_tile, l_clear, &
  l_flux_up_band, l_flux_down_band, &
  l_flux_up_clear_band, l_flux_down_clear_band, &
  l_blue_flux_surf

integer, intent(in), optional :: n_tile, n_cloud_layer, n_aer_mode

integer, intent(in), optional :: isolir, &
  i_cloud_representation, i_overlap, i_inhom, i_mcica_sampling, &
  i_st_water, i_cnv_water, i_st_ice, i_cnv_ice

! Local variables
integer :: i
integer :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_CONTROL_DEFAULTS'
character (len=errormessagelength) :: cmessage


! Logical options
if (present(l_rayleigh)) control%l_rayleigh = l_rayleigh
if (present(l_gas)) control%l_gas = l_gas
if (present(l_continuum)) control%l_continuum = l_continuum
if (present(l_cont_gen)) control%l_cont_gen = l_cont_gen
if (present(l_orog)) control%l_orog = l_orog
if (present(l_solvar)) control%l_solvar = l_solvar
if (present(l_rescale)) control%l_rescale = l_rescale
if (present(l_ir_source_quad)) control%l_ir_source_quad = l_ir_source_quad
if (present(l_mixing_ratio)) control%l_mixing_ratio = l_mixing_ratio
if (present(l_aerosol)) control%l_aerosol = l_aerosol
if (present(l_aerosol_mode)) control%l_aerosol_mode = l_aerosol_mode
if (present(l_aerosol_ccn)) control%l_aerosol_ccn = l_aerosol_ccn
if (present(l_tile)) control%l_tile = l_tile


! Integer options
if (present(isolir)) control%isolir = isolir
if (present(i_cloud_representation)) &
  control%i_cloud_representation = i_cloud_representation
if (present(i_overlap)) control%i_overlap = i_overlap
if (present(i_inhom)) control%i_inhom = i_inhom
if (present(i_mcica_sampling)) control%i_mcica_sampling = i_mcica_sampling
if (present(i_st_water)) control%i_st_water = i_st_water
if (present(i_cnv_water)) control%i_cnv_water = i_cnv_water
if (present(i_st_ice)) control%i_st_ice = i_st_ice
if (present(i_cnv_ice)) control%i_cnv_ice = i_cnv_ice


! Diagnostic options
if (present(l_flux_up_band)) control%l_flux_up_band = l_flux_up_band
if (present(l_flux_down_band)) control%l_flux_down_band = l_flux_down_band
if (present(l_flux_up_clear_band)) &
  control%l_flux_up_clear_band = l_flux_up_clear_band
if (present(l_flux_down_clear_band)) &
  control%l_flux_down_clear_band = l_flux_down_clear_band
if (present(l_blue_flux_surf)) control%l_blue_flux_surf = l_blue_flux_surf

if (present(l_clear)) control%l_clear = l_clear
control%l_clear = control%l_clear &
             .or. control%l_flux_up_clear_band &
             .or. control%l_flux_down_clear_band


! Defaults and checking
if (present(l_set_defaults)) then
  if (l_set_defaults) then

    ! Two-stream is the default
    call set_int_default(control%i_angular_integration, ip_two_stream)

    select case (control%i_angular_integration)
    case(ip_two_stream)

      ! Source-specific options
      select case(control%isolir)
      case(ip_solar)
        call set_int_default(control%i_2stream,        ip_pifm80)
        call set_int_default(control%i_scatter_method, ip_scatter_full)
      case(ip_infra_red)
        call set_int_default(control%i_2stream,        ip_elsasser)
        call set_int_default(control%i_scatter_method, ip_scatter_full)
        if (.not.present(l_ir_source_quad)) control%l_ir_source_quad = .true.
      case default
        cmessage = 'Radiative source option unrecognised.'
        ierr=i_err_fatal
        CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      end select

      ! Set general default options
      if (.not.present(l_rescale)) control%l_rescale = .true.
      if (.not.present(l_gas)) control%l_gas = .true.
      if (.not.present(l_continuum)) control%l_continuum = .true.
      if (.not.present(l_cont_gen)) control%l_cont_gen = .true.
      call set_int_default(control%first_band, 1)
      if (present(spectrum)) then
        call set_int_default(control%last_band, spectrum%basic%n_band)
      end if
      call set_int_default(control%n_channel, 1)
      call set_int_default(control%n_order_forward, 2)
      call set_int_default(control%i_gas_overlap, ip_overlap_k_eqv_scl)

      ! Consistent tiling options
      if (present(n_tile)) then
        if (n_tile < 1) control%l_tile = .false.
      end if

      ! Consistent cloud options and defaults
      call set_int_default(control%i_cloud_representation, ip_cloud_off)
      call set_int_default(control%i_overlap, ip_max_rand)
      call set_int_default(control%i_inhom, ip_homogeneous)
      call set_int_default(control%i_mcica_sampling, ip_mcica_optimal_sampling)
      if (present(n_cloud_layer)) then
        if (n_cloud_layer < 1) control%i_cloud_representation = ip_cloud_off
      end if
      select case(control%i_cloud_representation)
      case(ip_cloud_homogen, ip_cloud_ice_water)
        control%l_cloud = .true.
        if (.not.control%l_drop .and. .not.control%l_ice) then
          control%l_drop = .true.
          control%l_ice  = .true.
        end if
        if (control%i_inhom == ip_mcica) then
          control%i_cloud = ip_cloud_mcica
          if ( (control%i_scatter_method == ip_no_scatter_abs) .or. &
               (control%i_scatter_method == ip_no_scatter_ext) ) then
            control%i_solver       = ip_solver_no_scat
            control%i_solver_clear = ip_solver_no_scat
          else
            control%i_solver       = ip_solver_homogen_direct
            control%i_solver_clear = ip_solver_homogen_direct
          end if
        else
          if (control%i_scatter_method == ip_scatter_approx) then
            control%i_solver       = ip_solver_mix_app_scat
            control%i_solver_clear = ip_solver_homogen_direct
          else
            control%i_solver       = ip_solver_mix_direct_hogan
            control%i_solver_clear = ip_solver_homogen_direct
          end if
          if (control%i_overlap == ip_max_rand) then
            control%i_cloud = ip_cloud_mix_max
          else if (control%i_overlap == ip_exponential_rand) then
            control%i_cloud = ip_cloud_part_corr
          else
            control%i_cloud = ip_cloud_mix_random
          end if
        end if
      case(ip_cloud_conv_strat, ip_cloud_csiw)
        ! Not compatible with control%i_inhom == ip_mcica
        if (control%i_inhom == ip_mcica) control%i_inhom = ip_homogeneous
        control%l_cloud = .true.
        if (.not.control%l_drop .and. .not.control%l_ice) then
          control%l_drop = .true.
          control%l_ice  = .true.
        end if
        if (control%i_scatter_method == ip_scatter_approx) then
          control%i_solver       = ip_solver_triple_app_scat
          control%i_solver_clear = ip_solver_homogen_direct
        else
          control%i_solver       = ip_solver_triple_hogan
          control%i_solver_clear = ip_solver_homogen_direct
        end if
        if (control%i_overlap == ip_max_rand) then
          control%i_cloud = ip_cloud_triple
        else if (control%i_overlap == ip_exponential_rand) then
          control%i_cloud = ip_cloud_part_corr_cnv
        else
          control%i_cloud = ip_cloud_mix_random
        end if
      case(ip_cloud_off)
        ! No treatment of cloud
        control%l_cloud        = .false.
        control%l_drop         = .false.
        control%l_ice          = .false.
        control%l_microphysics = .false.
        control%i_cloud        = ip_cloud_clear
        if ( (control%i_scatter_method == ip_no_scatter_abs) .or. &
             (control%i_scatter_method == ip_no_scatter_ext) ) then
          control%i_solver       = ip_solver_no_scat
          control%i_solver_clear = ip_solver_no_scat
        else
          control%i_solver       = ip_solver_homogen_direct
          control%i_solver_clear = ip_solver_homogen_direct
        end if
      end select

      ! Consistent aerosol options
      if (present(n_aer_mode)) then
        if (n_aer_mode < 1) control%l_aerosol_mode = .false.
      end if

    case(ip_ir_gauss, ip_spherical_harmonic)

      ! Currently no checking done or defaults set

    case default

      cmessage = 'Method of angular integration unrecognised.'
      ierr=i_err_fatal
      CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)

    end select

  end if
end if


if (present(spectrum)) then
  ! Allocate band-by-band control options
  call allocate_control(control, spectrum)
  
  ! Set properties for individual bands.
  if (control%n_channel == control%last_band - control%first_band + 1) then
    do i = 1, control%n_channel
      control%map_channel(control%first_band + i-1) = i
    end do
  else
    do i = 1, spectrum%basic%n_band
      control%map_channel(i) = 1
    end do
  end if
  do i = 1, spectrum%basic%n_band
    control%weight_band(i)           = 1.0_RealK
    control%i_scatter_method_band(i) = control%i_scatter_method
    control%i_gas_overlap_band(i)    = control%i_gas_overlap
    if (any(spectrum%gas%i_scale_fnc(i,:) == ip_scale_ses2)) then
      ! If SES2 scaling is used in the band the overlap must also use SES2:
      control%i_gas_overlap_band(i)  = ip_overlap_mix_ses2
    end if
  end do
end if

contains

  subroutine set_int_default(i_option, i_default)
    implicit none
    integer, intent(inout) :: i_option
    integer, intent(in) :: i_default
    if (i_option == imdi) then
      i_option = i_default
    end if
  end subroutine set_int_default

end subroutine set_control
end module socrates_set_control
