! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Run the Socrates radiative transfer code

module socrates_runes

use rad_pcf, only: &
  ip_source_illuminate              => ip_solar, &
  ip_source_thermal                 => ip_infra_red, &
  ip_cloud_representation_off       => ip_cloud_off, &
  ip_cloud_representation_ice_water => ip_cloud_ice_water, &
  ip_cloud_representation_csiw      => ip_cloud_csiw, &
  ip_overlap_max_random             => ip_max_rand, &
  ip_overlap_random                 => ip_rand, &
  ip_overlap_exponential_random     => ip_exponential_rand, &
  ip_inhom_homogeneous              => ip_homogeneous, &
  ip_inhom_scaling                  => ip_scaling, &
  ip_inhom_mcica                    => ip_mcica, &
  ip_inhom_cairns                   => ip_cairns

implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_RUNES'
contains

subroutine runes(n_profile, n_layer, spectrum, spectrum_name, mcica_data, &
  n_cloud_layer, n_aer_mode, n_tile, &
  p_layer, t_layer, t_level, mass, density, &
  h2o, o3, &
  p_layer_1d, t_layer_1d, t_level_1d, mass_1d, density_1d, &
  h2o_1d, o3_1d, &
  co2_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, &
  o2_mix_ratio, so2_mix_ratio, cfc11_mix_ratio, cfc12_mix_ratio, &
  cfc113_mix_ratio, hcfc22_mix_ratio, hfc134a_mix_ratio, &
  t_ground, cos_zenith_angle, solar_irrad, orog_corr, &
  l_grey_albedo, grey_albedo, albedo_diff, albedo_dir, &
  l_tile, frac_tile, t_tile, albedo_diff_tile, albedo_dir_tile, &
  cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim, &
  cloud_frac_1d, conv_frac_1d, &
  liq_frac_1d, ice_frac_1d, liq_conv_frac_1d, ice_conv_frac_1d, &
  liq_mmr_1d, ice_mmr_1d, liq_conv_mmr_1d, ice_conv_mmr_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d, &
  cloud_vertical_decorr, conv_vertical_decorr, &
  cloud_horizontal_rsd, &
  layer_heat_capacity, layer_heat_capacity_1d, &
  i_source, i_cloud_representation, i_overlap, i_inhom, &
  i_mcica_sampling, i_st_water, i_st_ice, &
  rand_seed, &
  l_rayleigh, l_mixing_ratio, l_aerosol_mode, &
  l_invert, l_debug, i_profile_debug, &
  flux_direct, flux_down, flux_up, heating_rate, &
  flux_up_tile, flux_up_blue_tile, flux_direct_blue_surf, flux_down_blue_surf, &
  flux_direct_1d, flux_down_1d, flux_up_1d, heating_rate_1d, &
  flux_up_tile_1d, flux_up_blue_tile_1d, &
  total_cloud_cover)

use def_spectrum, only: StrSpecData
use def_mcica,    only: StrMcica
use def_control,  only: StrCtrl,  deallocate_control
use def_dimen,    only: StrDim
use def_atm,      only: StrAtm,   deallocate_atm
use def_bound,    only: StrBound, deallocate_bound
use def_cld,      only: StrCld,   deallocate_cld, deallocate_cld_prsc, &
                                  deallocate_cld_mcica
use def_aer,      only: StrAer,   deallocate_aer, deallocate_aer_prsc
use def_out,      only: StrOut,   deallocate_out

use socrates_set_spectrum, only: spectrum_array_name, spectrum_array, &
                                 mcica_spectrum_name, mcica_data_array

use socrates_set_control,   only: set_control
use socrates_set_dimen,     only: set_dimen
use socrates_set_atm,       only: set_atm
use socrates_set_bound,     only: set_bound
use socrates_set_cld,       only: set_cld
use socrates_set_cld_mcica, only: set_cld_mcica
use socrates_set_aer,       only: set_aer

use realtype_rd, only: RealK
use ereport_mod,  only: ereport
use errormessagelength_mod, only: errormessagelength
use rad_pcf, only: i_normal, i_err_fatal

implicit none

! Spectral data:
type (StrSpecData), intent(in), target, optional :: spectrum
character(len=*), intent(in), optional :: spectrum_name

! Mcica data
type (StrMcica), intent(in), target, optional :: mcica_data

integer, intent(in) :: n_profile
!   Number of columns to operate on
integer, intent(in) :: n_layer
!   Number of layers for radiation
integer, intent(in), optional :: n_tile
!   Number of surface tiles
integer, intent(in), optional :: n_cloud_layer
!   Number of potentially cloudy layers
integer, intent(in), optional :: n_aer_mode
!   Number of aerosol modes

real(RealK), intent(in), optional :: p_layer(n_profile, n_layer)
real(RealK), intent(in), optional :: p_layer_1d(n_layer)
!   Pressure at layer centres
real(RealK), intent(in), optional :: t_layer(n_profile, n_layer)
real(RealK), intent(in), optional :: t_layer_1d(n_layer)
!   Temperature at layer centres
real(RealK), intent(in), optional :: t_level(n_profile, 0:n_layer)
real(RealK), intent(in), optional :: t_level_1d(0:n_layer)
!   Temperature at layer boundaries
real(RealK), intent(in), optional :: mass(n_profile, n_layer)
real(RealK), intent(in), optional :: mass_1d(n_layer)
!   Mass of layer (kg m-2)
real(RealK), intent(in), optional :: density(n_profile, n_layer)
real(RealK), intent(in), optional :: density_1d(n_layer)
!   Density of layer (kg m-3)
real(RealK), intent(in), optional :: h2o(n_profile, n_layer)
real(RealK), intent(in), optional :: h2o_1d(n_layer)
!   Mass mixing ratio of water vapour
real(RealK), intent(in), optional :: o3(n_profile, n_layer)
real(RealK), intent(in), optional :: o3_1d(n_layer)
!   Mass mixing ratio of ozone

real(RealK), intent(in), optional :: &
  co2_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, &
  o2_mix_ratio, so2_mix_ratio, cfc11_mix_ratio, cfc12_mix_ratio, &
  cfc113_mix_ratio, hcfc22_mix_ratio, hfc134a_mix_ratio
!   Trace gas mass mixing ratios

real(RealK), intent(in), optional :: t_ground(n_profile)
!   Effective radiative temperature over whole grid-box
real(RealK), intent(in), optional :: cos_zenith_angle(n_profile)
!   Cosine of solar zenith angle
real(RealK), intent(in), optional :: solar_irrad(n_profile)
!   Solar irradiance at top-of-atmosphere (mean over timestep)
real(RealK), intent(in), optional :: orog_corr(n_profile)
!   Orographic correction factor

logical, intent(in), optional :: l_grey_albedo
!   Set a single grey albedo / emissivity for the surface
real(RealK), intent(in), optional :: grey_albedo
!   Grey surface albedo

real(RealK), intent(in), optional :: albedo_diff(:, :)
!   Spectral diffuse albedo (n_profile, n_band)
real(RealK), intent(in), optional :: albedo_dir(:, :)
!   Spectral direct albedo (n_profile, n_band)

logical, intent(in), optional :: l_tile
!   Use tiled surface properties
real(RealK), intent(in), optional :: frac_tile(:, :)
!   Tile fractions (n_profile, n_tile)
real(RealK), intent(in), optional :: albedo_diff_tile(:, :, :)
!   Diffuse tile albedo (n_profile, n_tile, n_band)
real(RealK), intent(in), optional :: albedo_dir_tile(:, :, :)
!   Direct tile albedo (n_profile, n_tile, n_band)
real(RealK), intent(in), optional :: t_tile(:, :)
!   Tile temperatures (n_profile, n_tile)

real(RealK), intent(in), dimension (n_profile, n_layer), optional :: &
  cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim
real(RealK), intent(in), dimension (n_layer), optional :: &
  cloud_frac_1d, conv_frac_1d, &
  liq_frac_1d, ice_frac_1d, liq_conv_frac_1d, ice_conv_frac_1d, &
  liq_mmr_1d, ice_mmr_1d, liq_conv_mmr_1d, ice_conv_mmr_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d
!   Liquid and ice cloud fractions, gridbox mean mixing ratios, and
!   effective dimensions

real(RealK), intent(in), optional :: cloud_vertical_decorr
!   Decorrelation pressure scale for cloud vertical overlap
real(RealK), intent(in), optional :: conv_vertical_decorr
!   Decorrelation pressure scale for convective cloud vertical overlap
real(RealK), intent(in), optional :: cloud_horizontal_rsd
!   Relative standard deviation of sub-grid cloud condensate

real(RealK), intent(in), optional :: layer_heat_capacity(n_profile, n_layer)
real(RealK), intent(in), optional :: layer_heat_capacity_1d(n_layer)
!   Heat capacity of layer

integer, intent(in), optional :: i_source
!   Select source of radiation
integer, intent(in), optional :: &
  i_cloud_representation, i_overlap, i_inhom, &
  i_mcica_sampling, i_st_water, i_st_ice
!   Select treatment of cloud
integer, intent(in), optional :: rand_seed(n_profile)
!   Random seed for cloud generator

logical, intent(in), optional :: l_rayleigh
!   Include Rayleigh scattering
logical, intent(in), optional :: l_mixing_ratio
!   Assume mass mixing ratios are with respect to dry mass
logical, intent(in), optional :: l_aerosol_mode
!   Include aerosol optical properties specified by mode

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information


! Output fields:
real(RealK), intent(out), optional :: flux_direct(n_profile, 0:n_layer)
real(RealK), intent(out), optional :: flux_direct_1d(0:n_layer)
!   Direct (unscattered) downwards flux (Wm-2)
real(RealK), intent(out), optional :: flux_down(n_profile, 0:n_layer)
real(RealK), intent(out), optional :: flux_down_1d(0:n_layer)
!   Downwards flux (Wm-2)
real(RealK), intent(out), optional :: flux_up(n_profile, 0:n_layer)
real(RealK), intent(out), optional :: flux_up_1d(0:n_layer)
!   Upwards flux (Wm-2)
real(RealK), intent(out), optional :: heating_rate(n_profile, n_layer)
real(RealK), intent(out), optional :: heating_rate_1d(n_layer)
!   Heating rate (Ks-1)
real(RealK), intent(out), optional :: flux_up_tile(:, :) ! (n_profile, n_tile)
real(RealK), intent(out), optional :: flux_up_tile_1d(:) ! (n_tile)
!   Upwards flux on tiles (Wm-2)
real(RealK), intent(out), optional :: flux_up_blue_tile(:, :)
real(RealK), intent(out), optional :: flux_up_blue_tile_1d(:)
!   Upwards blue flux on tiles (Wm-2)
real(RealK), intent(out), optional :: flux_direct_blue_surf(n_profile)
!   Direct blue flux at the surface
real(RealK), intent(out), optional :: flux_down_blue_surf(n_profile)
!   Total downward blue flux at the surface
real(RealK), intent(out), optional :: total_cloud_cover(n_profile)
!   Total cloud cover

! Spectral data:
type(StrSpecData), pointer :: spec => null()

! Mcica data:
type(StrMcica), pointer :: mcica => null()
type(StrMcica), target :: mcica_dummy

! Controlling options:
type(StrCtrl) :: control

! Dimensions:
type(StrDim) :: dimen

! Atmospheric properties:
type(StrAtm) :: atm

! Boundary conditions:
type(StrBound) :: bound

! Cloud properties:
type(StrCld) :: cld

! Aerosol properties:
type(StrAer) :: aer

! Output fields from core radiation code:
type(StrOut) :: radout

integer :: l, i, id_spec, id_mcica
!   Loop variables
logical :: l_inv
!   Local logical for field inversion
logical :: l_blue_flux_surf
!   Output blue fluxes if requested
real(RealK) :: flux_divergence(n_profile, n_layer)
!   Flux divergence across layer (Wm-2)

integer :: ierr = i_normal
character (len=errormessagelength) :: cmessage
character (len=*), parameter :: RoutineName = 'RUNES'


if (present(spectrum_name)) then
  do id_spec=1, size(spectrum_array)
    if (spectrum_array_name(id_spec) == spectrum_name) exit
    if (id_spec == size(spectrum_array)) then
      cmessage = 'Spectrum name not found.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if
  end do
  spec => spectrum_array(id_spec)
  if ( (i_cloud_representation /= ip_cloud_representation_off) .and. &
       (i_inhom == ip_inhom_mcica) ) then
    if (allocated(mcica_data_array).and.allocated(mcica_spectrum_name)) then
      do id_mcica=1, size(mcica_data_array)
        if (mcica_spectrum_name(i_source, id_mcica) == spectrum_name) exit
        if (id_mcica == size(mcica_data_array)) then
          cmessage = 'Spectrum name not associated with MCICA data.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if
      end do
      mcica => mcica_data_array(id_mcica)
    else
      cmessage = 'MCICA data has not been read in correctly.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if
  else
    mcica => mcica_dummy
  end if
else if (present(spectrum)) then
  spec => spectrum
  if ( (i_cloud_representation /= ip_cloud_representation_off) .and. &
       (i_inhom == ip_inhom_mcica) ) then
    if (present(mcica_data)) then
      mcica => mcica_data
    else
      cmessage = 'No mcica_data object has been passed to socrates_runes.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if
  else
    mcica => mcica_dummy
  end if
else
  cmessage = 'No spectrum name or object supplied.'
  ierr=i_err_fatal
  call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
end if

if (present(flux_up_blue_tile)     .or. &
    present(flux_up_blue_tile_1d)  .or. &
    present(flux_direct_blue_surf) .or. &
    present(flux_down_blue_surf) ) then
  l_blue_flux_surf = .true.
else
  l_blue_flux_surf = .false.
end if

call set_control(control, spec, &
  isolir                 = i_source, &
  l_rayleigh             = l_rayleigh, &
  l_mixing_ratio         = l_mixing_ratio, &
  l_aerosol_mode         = l_aerosol_mode, &
  l_tile                 = l_tile, &
  l_blue_flux_surf       = l_blue_flux_surf, &
  n_tile                 = n_tile, &
  n_cloud_layer          = n_cloud_layer, &
  n_aer_mode             = n_aer_mode, &
  i_cloud_representation = i_cloud_representation, &
  i_overlap              = i_overlap, &
  i_inhom                = i_inhom, &
  i_mcica_sampling       = i_mcica_sampling, &
  i_st_water             = i_st_water, &
  i_st_ice               = i_st_ice, &
  l_set_defaults         = .true.)

call set_dimen(dimen, control, n_profile, n_layer, &
  mcica_data    = mcica, &
  n_tile        = n_tile, &
  n_cloud_layer = n_cloud_layer, &
  n_aer_mode    = n_aer_mode )

call set_atm(atm, dimen, spec, n_profile, n_layer, &
  p_layer           = p_layer,           &
  t_layer           = t_layer,           &
  mass              = mass,              &
  density           = density,           &
  t_level           = t_level,           &
  h2o               = h2o,               &
  o3                = o3,                &
  p_layer_1d        = p_layer_1d,        &
  t_layer_1d        = t_layer_1d,        &
  mass_1d           = mass_1d,           &
  density_1d        = density_1d,        &
  t_level_1d        = t_level_1d,        &
  h2o_1d            = h2o_1d,            &
  o3_1d             = o3_1d,             &
  co2_mix_ratio     = co2_mix_ratio,     &
  n2o_mix_ratio     = n2o_mix_ratio,     &
  ch4_mix_ratio     = ch4_mix_ratio,     &
  o2_mix_ratio      = o2_mix_ratio,      &
  so2_mix_ratio     = so2_mix_ratio,     &
  cfc11_mix_ratio   = cfc11_mix_ratio,   &
  cfc12_mix_ratio   = cfc12_mix_ratio,   &
  cfc113_mix_ratio  = cfc113_mix_ratio,  &
  hcfc22_mix_ratio  = hcfc22_mix_ratio,  &
  hfc134a_mix_ratio = hfc134a_mix_ratio, &
  l_invert          = l_invert,          &
  l_debug           = l_debug,           &
  i_profile_debug   = i_profile_debug )

call set_bound(bound, control, dimen, spec, n_profile, &
  n_tile           = n_tile, &
  t_ground         = t_ground, &
  cos_zenith_angle = cos_zenith_angle, &
  solar_irrad      = solar_irrad, &
  orog_corr        = orog_corr, &
  l_grey_albedo    = l_grey_albedo, &
  grey_albedo      = grey_albedo, &
  albedo_diff      = albedo_diff, &
  albedo_dir       = albedo_dir, &
  frac_tile        = frac_tile, &
  t_tile           = t_tile, &
  albedo_diff_tile = albedo_diff_tile, &
  albedo_dir_tile  = albedo_dir_tile, &
  l_debug          = l_debug, &
  i_profile_debug  = i_profile_debug )

call set_cld(cld, control, dimen, spec, atm, &
  cloud_frac            = cloud_frac, &
  liq_frac              = liq_frac, &
  ice_frac              = ice_frac, &
  liq_mmr               = liq_mmr, &
  ice_mmr               = ice_mmr, &
  liq_dim               = liq_dim, &
  ice_dim               = ice_dim, &
  cloud_frac_1d         = cloud_frac_1d, &
  liq_frac_1d           = liq_frac_1d, &
  ice_frac_1d           = ice_frac_1d, &
  liq_mmr_1d            = liq_mmr_1d, &
  ice_mmr_1d            = ice_mmr_1d, &
  liq_dim_1d            = liq_dim_1d, &
  ice_dim_1d            = ice_dim_1d, &
  cloud_vertical_decorr = cloud_vertical_decorr, &
  l_invert              = l_invert, &
  l_debug               = l_debug, &
  i_profile_debug       = i_profile_debug )

call set_cld_mcica(cld, mcica, control, dimen, spec, atm, &
  rand_seed            = rand_seed, &
  cloud_horizontal_rsd = cloud_horizontal_rsd)

call set_aer(aer, control, dimen, spec, n_profile, n_layer)

! DEPENDS ON: radiance_calc
call radiance_calc(control, dimen, spec, atm, cld, aer, bound, radout)

if (present(l_invert)) then
  l_inv = l_invert
else
  l_inv = .false.
end if

! set heating rates and diagnostics
if (present(heating_rate).or.present(heating_rate_1d)) then
  if (l_inv) then
    do i=1, n_layer
      do l=1, n_profile
        flux_divergence(l, n_layer-i+1) = &
          sum(radout%flux_down(l, i-1, 1:control%n_channel)) - &
          sum(radout%flux_down(l, i,   1:control%n_channel)) + &
          sum(radout%flux_up(  l, i,   1:control%n_channel)) - &
          sum(radout%flux_up(  l, i-1, 1:control%n_channel))
      end do
    end do
  else
    do i=1, n_layer
      do l=1, n_profile
        flux_divergence(l, i) = &
          sum(radout%flux_down(l, i-1, 1:control%n_channel)) - &
          sum(radout%flux_down(l, i,   1:control%n_channel)) + &
          sum(radout%flux_up(  l, i,   1:control%n_channel)) - &
          sum(radout%flux_up(  l, i-1, 1:control%n_channel))
      end do
    end do
  end if
  if (present(heating_rate)) then
    if (present(layer_heat_capacity)) then
      heating_rate = flux_divergence / layer_heat_capacity
    else if (present(layer_heat_capacity_1d)) then
      do i=1, n_layer
        heating_rate(1:n_profile, i) = &
          flux_divergence(1:n_profile, i) / layer_heat_capacity_1d(i)
      end do
    else
      ! Just return the flux_divergence if no heat capacity supplied
      heating_rate = flux_divergence
    end if
  end if
  if (present(heating_rate_1d)) then
    if (n_profile == 1) then
      if (present(layer_heat_capacity)) then
        heating_rate_1d(1:n_layer) = &
          flux_divergence(1, 1:n_layer) / layer_heat_capacity(1, 1:n_layer)
      else if (present(layer_heat_capacity_1d)) then
        heating_rate_1d(1:n_layer) = &
          flux_divergence(1, 1:n_layer) / layer_heat_capacity_1d(1:n_layer)
      else
        heating_rate_1d(1:n_layer) = flux_divergence(1, 1:n_layer)
      end if
    else
      if (present(layer_heat_capacity)) then
        heating_rate_1d = &
          sum(flux_divergence / layer_heat_capacity, 1) &
          / real(n_profile, RealK)
      else if (present(layer_heat_capacity_1d)) then
        heating_rate_1d = &
          (sum(flux_divergence, 1) / layer_heat_capacity_1d) &
          / real(n_profile, RealK)
      else
        heating_rate_1d = sum(flux_divergence, 1) / real(n_profile, RealK)
      end if
    end if
  end if
end if

call sum_flux_channels(flux_direct, flux_direct_1d, radout%flux_direct)
call sum_flux_channels(flux_down, flux_down_1d, radout%flux_down)
call sum_flux_channels(flux_up, flux_up_1d, radout%flux_up)
call sum_tile_channels(flux_up_tile, flux_up_tile_1d, radout%flux_up_tile)
call sum_tile_channels(flux_up_blue_tile, flux_up_blue_tile_1d, &
                       radout%flux_up_blue_tile)
if (present(flux_direct_blue_surf)) &
  flux_direct_blue_surf = radout%flux_direct_blue_surf(1:n_profile)
if (present(flux_down_blue_surf)) &
  flux_down_blue_surf = radout%flux_down_blue_surf(1:n_profile)
if (present(total_cloud_cover)) then
  if (control%i_cloud_representation == ip_cloud_representation_off) then
    total_cloud_cover = 0.0_RealK
  else
    if (control%i_inhom == ip_inhom_mcica) then
      total_cloud_cover = real(cld%n_subcol_cld(1:n_profile), RealK) &
                        / real(mcica%n_subcol_gen, RealK)
    else
      total_cloud_cover = radout%tot_cloud_cover(1:n_profile)
    end if
  end if
end if

call deallocate_out(radout)
call deallocate_aer_prsc(aer)
call deallocate_aer(aer)
call deallocate_cld_mcica(cld)
call deallocate_cld_prsc(cld)
call deallocate_cld(cld)
call deallocate_bound(bound)
call deallocate_atm(atm)
call deallocate_control(control)

contains

  subroutine sum_flux_channels(field, field_1d, field_channels)
  
  implicit none
  
  real(RealK), intent(out), optional :: field(n_profile, 0:n_layer)
  real(RealK), intent(out), optional :: field_1d(0:n_layer)
  real(RealK), intent(in) :: field_channels(:, 0:, :)
  
  if (present(field)) then
    if (l_inv) then
      do i=0, n_layer
        do l=1, n_profile
          field(l, n_layer-i) = &
            sum(field_channels(l, i, 1:control%n_channel))
        end do
      end do
    else
      do i=0, n_layer
        do l=1, n_profile
          field(l, i) = &
            sum(field_channels(l, i, 1:control%n_channel))
        end do
      end do
    end if
  end if
  if (present(field_1d)) then
    if (n_profile == 1) then
      if (l_inv) then
        do i=0, n_layer
          field_1d(n_layer-i) = &
            sum(field_channels(1, i, 1:control%n_channel))
        end do
      else
        do i=0, n_layer
          field_1d(i) = &
            sum(field_channels(1, i, 1:control%n_channel))
        end do
      end if
    else
      if (l_inv) then
        do i=0, n_layer
          field_1d(n_layer-i) = &
            sum(field_channels(1:n_profile, i, 1:control%n_channel)) &
            / real(n_profile, RealK)
        end do
      else
        do i=0, n_layer
          field_1d(i) = &
            sum(field_channels(1:n_profile, i, 1:control%n_channel)) &
            / real(n_profile, RealK)
        end do
      end if
    end if
  end if
  
  end subroutine sum_flux_channels
  
  
  subroutine sum_tile_channels(field, field_1d, field_channels)
  
  implicit none
  
  real(RealK), intent(out), optional :: field(:, :)
  real(RealK), intent(out), optional :: field_1d(:)
  real(RealK), intent(in) :: field_channels(:, :, :)
  
  integer :: ll
  
  if (present(field)) then
    field(:, :) = 0.0_RealK
    if (present(n_tile).and.control%l_tile) then
      do i=1, n_tile
        do ll=1, bound%n_point_tile
          l = bound%list_tile(ll)
          field(l, i) = sum(field_channels(ll, i, 1:control%n_channel))
        end do
      end do
    end if
  end if
  if (present(field_1d)) then
    field_1d(:) = 0.0_RealK
    if (present(n_tile).and.control%l_tile) then
      if (bound%n_point_tile == 1) then
        do i=1, n_tile
          field_1d(i) = sum(field_channels(1, i, 1:control%n_channel))
        end do
      else if (bound%n_point_tile > 1) then
        do i=1, n_tile
          field_1d(i) = &
            sum(field_channels(1:bound%n_point_tile, i, 1:control%n_channel)) &
            / real(bound%n_point_tile, RealK)
        end do
      end if
    end if
  end if
  
  end subroutine sum_tile_channels

end subroutine runes
end module socrates_runes
