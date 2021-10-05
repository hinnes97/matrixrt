! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Calculate the parameters for external illumination of the atmosphere

module socrates_illuminate

use def_orbit, only: &
  ip_elements_user, &
  ip_elements_earth_fixed, &
  ip_elements_earth_secular_variation, &
  ip_spin_user, &
  ip_spin_earth_day, &
  ip_spin_fixed_sun

implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_ILLUMINATE'
contains

subroutine illuminate( &
    l_calendar_360, l_observer, l_stellar_position, l_stellar_angle, &
    n_profile, i_elements, i_spin, &
    year, month, day_of_month, day_of_year, second_of_day, length_of_timestep, &
    epoch, eccentricity, eccentricity_inc, arg_periapsis, arg_periapsis_inc, &
    obliquity, obliquity_inc, semimajor_axis, semimajor_axis_inc, &
    mean_anomaly, mean_anomaly_inc, hour_angle, hour_angle_inc, &
    fixed_zenith_angle, fixed_azimuth_angle, observer_lon, observer_lat, &
    latitude, longitude, stellar_constant, &
    sin_stellar_declination, stellar_eqn_of_time, stellar_constant_scaling, &
    sin_observer_declination, observer_eqn_of_time, observed_orbital_phase, &
    cos_zenith_angle, lit_fraction, stellar_irradiance, &
    observer_cos_zen_ang, observed_fraction)

  use realtype_rd, only: RealK
  use def_orbit, only: StrOrbit, set_orbit
  use solpos_mod, only: solpos
  use solang_mod, only: solang

  implicit none

  logical, intent(in), optional :: l_calendar_360, l_observer
  logical, intent(in), optional :: l_stellar_position, l_stellar_angle
  integer, intent(in), optional :: n_profile
  integer, intent(in), optional :: i_elements, i_spin
  integer, intent(in), optional :: year, month, day_of_month, day_of_year
  real(RealK), intent(in), optional :: second_of_day, length_of_timestep
  real(RealK), intent(in), optional :: &
    epoch, eccentricity, eccentricity_inc, arg_periapsis, arg_periapsis_inc, &
    obliquity, obliquity_inc, semimajor_axis, semimajor_axis_inc, &
    mean_anomaly, mean_anomaly_inc, hour_angle, hour_angle_inc, &
    fixed_zenith_angle, fixed_azimuth_angle, observer_lon, observer_lat
  real(RealK), intent(in), optional :: latitude(:)
  real(RealK), intent(in), optional :: longitude(:)
  real(RealK), intent(in), optional :: stellar_constant

  real(RealK), intent(inout), optional :: sin_stellar_declination
  real(RealK), intent(inout), optional :: stellar_eqn_of_time
  real(RealK), intent(inout), optional :: stellar_constant_scaling
  real(RealK), intent(inout), optional :: sin_observer_declination
  real(RealK), intent(inout), optional :: observer_eqn_of_time
  real(RealK), intent(inout), optional :: observed_orbital_phase

  real(RealK), intent(out), optional :: cos_zenith_angle(:)
  real(RealK), intent(out), optional :: lit_fraction(:)
  real(RealK), intent(out), optional :: stellar_irradiance(:)
  real(RealK), intent(out), optional :: observer_cos_zen_ang(:)
  real(RealK), intent(out), optional :: observed_fraction(:)

  ! Local variables
  type(StrOrbit) :: orbit
  logical :: l_solpos, l_solang
  real(RealK) :: eqt, sindec, scs, sindec_obs, eqt_obs, phase_obs


  call set_orbit(orbit, &
    i_elements          = i_elements, &
    i_spin              = i_spin, &
    epoch               = epoch, &
    eccentricity        = eccentricity, &
    eccentricity_inc    = eccentricity_inc, &
    arg_periapsis       = arg_periapsis, &
    arg_periapsis_inc   = arg_periapsis_inc, &
    obliquity           = obliquity, &
    obliquity_inc       = obliquity_inc, &
    semimajor_axis      = semimajor_axis, &
    semimajor_axis_inc  = semimajor_axis_inc, &
    mean_anomaly        = mean_anomaly, &
    mean_anomaly_inc    = mean_anomaly_inc, &
    hour_angle          = hour_angle, &
    hour_angle_inc      = hour_angle_inc, &
    fixed_zenith_angle  = fixed_zenith_angle, &
    fixed_azimuth_angle = fixed_azimuth_angle, &
    observer_lon        = observer_lon, &
    observer_lat        = observer_lat )

  if (present(l_stellar_position)) then
    l_solpos = l_stellar_position
  else
    l_solpos = .false.
  end if

  if (l_solpos) then
    ! Find stellar position from the planet surface
    call solpos(orbit, year, day_of_year, second_of_day, length_of_timestep, &
      sindec, eqt, scs, l_observer, l_calendar_360, &
      sindec_obs, eqt_obs, phase_obs)
    if (present(sin_stellar_declination))  sin_stellar_declination  = sindec
    if (present(stellar_eqn_of_time))      stellar_eqn_of_time      = eqt
    if (present(stellar_constant_scaling)) stellar_constant_scaling = scs
    if (present(sin_observer_declination)) sin_observer_declination = sindec_obs
    if (present(observer_eqn_of_time))     observer_eqn_of_time     = eqt_obs
    if (present(observed_orbital_phase))   observed_orbital_phase   = phase_obs
  end if

  if (present(l_stellar_angle)) then
    l_solang = l_stellar_angle
  else
    l_solang = .false.
  end if

  if (l_solang) then
    if (present(sin_stellar_declination))  sindec = sin_stellar_declination
    if (present(stellar_eqn_of_time))      eqt    = stellar_eqn_of_time
    if (present(stellar_constant_scaling)) scs    = stellar_constant_scaling
    ! Find stellar zenith angle and lit fraction of timestep for each gridpoint
    call solang(orbit, second_of_day, length_of_timestep, &
      sindec, eqt, latitude, longitude, n_profile, &
      lit_fraction, cos_zenith_angle)
    if (present(stellar_irradiance).and.present(stellar_constant)) then
      stellar_irradiance = stellar_constant*scs*lit_fraction
    end if
  end if

  if (present(l_observer)) then
    if (l_observer) then
      if (present(sin_observer_declination)) &
                                         sindec_obs = sin_observer_declination
      if (present(observer_eqn_of_time)) eqt_obs    = observer_eqn_of_time
      ! Find observer zenith angle and observed fraction of timestep
      call solang(orbit, second_of_day, length_of_timestep, &
        sindec_obs, eqt_obs, latitude, longitude, n_profile, &
        observed_fraction, observer_cos_zen_ang)
    end if
  end if

end subroutine illuminate
end module socrates_illuminate
