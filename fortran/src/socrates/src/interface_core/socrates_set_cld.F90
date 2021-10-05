! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates cloud type
!
!------------------------------------------------------------------------------
module socrates_set_cld
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_CLD'
contains

subroutine set_cld(cld, control, dimen, spectrum, atm, &
  cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim, &
  cloud_frac_1d, conv_frac_1d, &
  liq_frac_1d, ice_frac_1d, liq_conv_frac_1d, ice_conv_frac_1d, &
  liq_mmr_1d, ice_mmr_1d, liq_conv_mmr_1d, ice_conv_mmr_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d, &
  cloud_vertical_decorr, conv_vertical_decorr, &
  l_invert, l_debug, i_profile_debug)

use def_cld,      only: StrCld, allocate_cld, allocate_cld_prsc
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use def_atm,      only: StrAtm
use realtype_rd,  only: RealK
use rad_pcf,      only: &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat, ip_cloud_csiw, &
  ip_clcmp_st_water, ip_clcmp_st_ice, ip_clcmp_cnv_water, ip_clcmp_cnv_ice, &
  ip_phase_water, ip_phase_ice, ip_cloud_type_homogen, &
  ip_cloud_type_water, ip_cloud_type_ice, &
  ip_cloud_type_strat, ip_cloud_type_conv, &
  ip_cloud_type_sw, ip_cloud_type_si, ip_cloud_type_cw, ip_cloud_type_ci, &
  ip_drop_unparametrized, ip_ice_unparametrized, &
  i_normal, i_err_fatal
use ereport_mod,  only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none


! Cloud properties:
type(StrCld),      intent(out) :: cld

! Control options:
type(StrCtrl),     intent(in)  :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

! Atmospheric properties:
type(StrAtm),      intent(in)  :: atm

real(RealK), intent(in), dimension(:, :), optional :: &
  cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim
!   Liquid and ice cloud fractions, gridbox mean mixing ratios, and
!   effective dimensions

real(RealK), intent(in), dimension(:), optional :: &
  cloud_frac_1d, conv_frac_1d, &
  liq_frac_1d, ice_frac_1d, liq_conv_frac_1d, ice_conv_frac_1d, &
  liq_mmr_1d, ice_mmr_1d, liq_conv_mmr_1d, ice_conv_mmr_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d
!   Liquid and ice cloud fractions, gridbox mean mixing ratios, and
!   effective dimensions input as 1d fields

real(RealK), intent(in), optional :: cloud_vertical_decorr
!   Decorrelation pressure scale for cloud vertical overlap
real(RealK), intent(in), optional :: conv_vertical_decorr
!   Decorrelation pressure scale for convective cloud vertical overlap

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information


! Local variables
integer :: i, j, k, l
!   Loop variables
integer :: i_phase, i_param_type, n_cloud_parameter
!   Working variables
integer :: i_cloud_type(dimen%nd_cloud_component)
!   Types of cloud to which each component contributes
real(RealK), dimension(dimen%nd_profile, dimen%id_cloud_top:dimen%nd_layer) :: &
  frac, frac_liq, frac_ice
!   Cloud fraction working arrays

logical :: l_inv
!   Local flag to invert fields in the vertical

real(RealK) :: eps = EPSILON(1.0)
real(RealK) :: min_cloud_fraction = 0.001

integer                      :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_CLD'
character (len=errormessagelength) :: cmessage

! Functions called
integer, external :: set_n_cloud_parameter


! Allocate structure for the core radiation code interface
call allocate_cld(cld, dimen, spectrum)
call allocate_cld_prsc(cld, dimen, spectrum)

if (.not.control%l_cloud) then
  return
end if

if (present(l_invert)) then
  l_inv = l_invert
else
  l_inv = .false.
end if

!------------------------------------------------------------------------------
! Set properties of condensed components
!------------------------------------------------------------------------------

if (control%l_ice .and. control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_st_ice
  case (ip_cloud_conv_strat, ip_cloud_csiw)
    cld%n_condensed = 4
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_st_ice
    cld%type_condensed(3) = ip_clcmp_cnv_water
    cld%type_condensed(4) = ip_clcmp_cnv_ice
  end select
else if (control%l_ice .and. .not.control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water)
    cld%n_condensed = 1
    cld%type_condensed(1) = ip_clcmp_st_ice
  case (ip_cloud_conv_strat, ip_cloud_csiw)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_ice
    cld%type_condensed(2) = ip_clcmp_cnv_ice
  end select
else if (.not.control%l_ice .and. control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water)
    cld%n_condensed = 1
    cld%type_condensed(1) = ip_clcmp_st_water
  case (ip_cloud_conv_strat, ip_cloud_csiw)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_cnv_water
  end select
else
  cmessage = 'Cloud on, but no condensed components included.'
  ierr=i_err_fatal
  call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
end if

do i=1, cld%n_condensed
  select case (cld%type_condensed(i))
  case (ip_clcmp_st_water)
    i_phase = ip_phase_water
    i_param_type = control%i_st_water
  case (ip_clcmp_st_ice)
    i_phase = ip_phase_ice
    i_param_type = control%i_st_ice
  case (ip_clcmp_cnv_water)
    i_phase = ip_phase_water
    i_param_type = control%i_cnv_water
  case (ip_clcmp_cnv_ice)
    i_phase = ip_phase_ice
    i_param_type = control%i_cnv_ice
  end select

  select case (i_phase)
  case (ip_phase_water)
    if (i_param_type <= 0) then
      cld%i_condensed_param(i) = ip_drop_unparametrized
      cmessage = 'Prescribed liquid cloud not yet implemented.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (i_param_type > spectrum%dim%nd_drop_type) then
      cmessage = 'Liquid cloud type outside allowed range.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (spectrum%drop%l_drop_type(i_param_type)) then
      ! Take parametrisation from spectral file
      cld%i_condensed_param(i) = spectrum%drop%i_drop_parm(i_param_type)
      cld%condensed_n_phf(i) = spectrum%drop%n_phf(i_param_type)
      ! DEPENDS ON: set_n_cloud_parameter
      n_cloud_parameter = set_n_cloud_parameter( cld%i_condensed_param(i), &
        cld%type_condensed(i), cld%condensed_n_phf(i) )
      do j=1, spectrum%basic%n_band
        do k=1, n_cloud_parameter
          cld%condensed_param_list(k, i, j) &
            = spectrum%drop%parm_list(k, j, i_param_type)
        end do
      end do
      ! Assign droplet mass mixing ratio and effective radius
      select case (cld%type_condensed(i))
      case (ip_clcmp_st_water)
        if ((present(liq_mmr).or.present(liq_mmr_1d)) .and. &
            (present(liq_dim).or.present(liq_dim_1d))) then
          call set_cld_field(cld%condensed_mix_ratio(:, :, i), &
                             liq_mmr, liq_mmr_1d)
          call set_cld_field(cld%condensed_dim_char(:, :, i), &
                             liq_dim, liq_dim_1d)
        else
          cmessage = 'Liquid MMR and effective radius not provided.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if
      case (ip_clcmp_cnv_water)
        if ((present(liq_conv_mmr).or.present(liq_conv_mmr_1d)) .and. &
            (present(liq_conv_dim).or.present(liq_conv_dim_1d))) then
          call set_cld_field(cld%condensed_mix_ratio(:, :, i), &
                             liq_conv_mmr, liq_conv_mmr_1d)
          call set_cld_field(cld%condensed_dim_char(:, :, i), &
                             liq_conv_dim, liq_conv_dim_1d)
        else
          cmessage = 'Convective liquid MMR and effective radius not provided.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if        
      end select
      ! Constrain effective radius to be within parametrisation bounds
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          cld%condensed_dim_char(l, k, i) = min( max( &
            cld%condensed_dim_char(l, k, i), &
            spectrum%drop%parm_min_dim(i_param_type) ), &
            spectrum%drop%parm_max_dim(i_param_type) )
        end do
      end do
    else
      cmessage = 'Liquid cloud type not in spectral file.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if
  case (ip_phase_ice)
    if (i_param_type <= 0) then
      cld%i_condensed_param(i) = ip_ice_unparametrized
      cmessage = 'Prescribed ice cloud not yet implemented.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (i_param_type > spectrum%dim%nd_ice_type) then
      cmessage = 'Ice cloud type outside allowed range.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (spectrum%ice%l_ice_type(i_param_type)) then
      ! Take parametrisation from spectral file
      cld%i_condensed_param(i) = spectrum%ice%i_ice_parm(i_param_type)
      cld%condensed_n_phf(i) = spectrum%ice%n_phf(i_param_type)
      n_cloud_parameter = set_n_cloud_parameter( cld%i_condensed_param(i), &
        cld%type_condensed(i), cld%condensed_n_phf(i) )
      do j=1, spectrum%basic%n_band
        do k=1, n_cloud_parameter
          cld%condensed_param_list(k, i, j) &
            = spectrum%ice%parm_list(k, j, i_param_type)
        end do
      end do
      ! Assign ice mass mixing ratio and effective dimension
      select case (cld%type_condensed(i))
      case (ip_clcmp_st_ice)
        if (present(ice_mmr).or.present(ice_mmr_1d)) then
          call set_cld_field(cld%condensed_mix_ratio(:, :, i), &
                             ice_mmr, ice_mmr_1d)
        else
          cmessage = 'Ice MMR not provided.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if
        if (present(ice_dim).or.present(ice_dim_1d)) then
          call set_cld_field(cld%condensed_dim_char(:, :, i), &
                             ice_dim, ice_dim_1d)
        else
          call set_ice_dim()
        end if
      case (ip_clcmp_cnv_ice)
        if (present(ice_conv_mmr).or.present(ice_conv_mmr_1d)) then
          call set_cld_field(cld%condensed_mix_ratio(:, :, i), &
                             ice_conv_mmr, ice_conv_mmr_1d)
        else
          cmessage = 'Convective ice MMR not provided.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if
        if (present(ice_conv_dim).or.present(ice_conv_dim_1d)) then
          call set_cld_field(cld%condensed_dim_char(:, :, i), &
                             ice_conv_dim, ice_conv_dim_1d)
        else
          call set_ice_dim()
        end if
      end select
      ! Constrain effective dimension to be within parametrisation bounds
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          cld%condensed_dim_char(l, k, i) = min( max( &
            cld%condensed_dim_char(l, k, i), &
            spectrum%ice%parm_min_dim(i_param_type) ), &
            spectrum%ice%parm_max_dim(i_param_type) )
        end do
      end do
    else
      cmessage = 'Ice cloud type not in spectral file.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if
  end select
end do

! Set the decorrelation scalings for cloud vertical overlap
if (present(cloud_vertical_decorr)) then
  cld%dp_corr_strat = cloud_vertical_decorr
else
  cld%dp_corr_strat = 0.0_RealK
end if
if (present(conv_vertical_decorr)) then
  cld%dp_corr_conv  = conv_vertical_decorr
else
  cld%dp_corr_conv = 0.0_RealK
end if


!------------------------------------------------------------------------------
! Set cloud amounts and convert mixing ratios to in-cloud values
!------------------------------------------------------------------------------

! Set cloud fractions
select case (control%i_cloud_representation)
case (ip_cloud_homogen)
  cld%n_cloud_type = 1
  do i = 1, cld%n_condensed
    i_cloud_type(i) = ip_cloud_type_homogen
  end do
  if (present(cloud_frac).or.present(cloud_frac_1d)) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_homogen), &
                       cloud_frac, cloud_frac_1d)
  else
    cmessage = 'Cloud fraction not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if
case (ip_cloud_ice_water)
  cld%n_cloud_type = 2
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      i_cloud_type(i) = ip_cloud_type_water
    case (ip_clcmp_st_ice)
      i_cloud_type(i) = ip_cloud_type_ice
    end select
  end do
  if ((present(liq_frac).or.present(liq_frac_1d)) .and. &
      (present(ice_frac).or.present(ice_frac_1d)) .and. &
      (present(cloud_frac).or.present(cloud_frac_1d))) then
    call set_cld_field(frac_liq, liq_frac, liq_frac_1d)
    call set_cld_field(frac_ice, ice_frac, ice_frac_1d)
    call set_cld_field(frac, cloud_frac, cloud_frac_1d)
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        if (frac_liq(l, k) + frac_ice(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          cld%frac_cloud(l, k, ip_cloud_type_water) = &
            frac(l, k)*frac_liq(l, k) / (frac_liq(l, k)+frac_ice(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_ice) = &
            frac(l, k)*frac_ice(l, k) / (frac_liq(l, k)+frac_ice(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_water) = 0.0_RealK
          cld%frac_cloud(l, k, ip_cloud_type_ice) = 0.0_RealK
        end if
      end do
    end do
  else if ((present(liq_frac).or.present(liq_frac_1d)) .and. &
           (present(ice_frac).or.present(ice_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_water), &
                       liq_frac, liq_frac_1d)
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_ice), &
                       ice_frac, ice_frac_1d)
  else
    cmessage = 'Liquid and ice cloud fractions not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if
case (ip_cloud_conv_strat)
  cld%n_cloud_type = 2
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      i_cloud_type(i) = ip_cloud_type_strat
    case (ip_clcmp_st_ice)
      i_cloud_type(i) = ip_cloud_type_strat
    case (ip_clcmp_cnv_water)
      i_cloud_type(i) = ip_cloud_type_conv
    case (ip_clcmp_cnv_ice)
      i_cloud_type(i) = ip_cloud_type_conv
    end select
  end do
  if ((present(cloud_frac).or.present(cloud_frac_1d)) .and. &
      (present(conv_frac).or.present(conv_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_strat), &
                       cloud_frac, cloud_frac_1d)
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_conv), &
                       conv_frac, conv_frac_1d)
  else
    cmessage = 'Cloud fraction and convective cloud fraction not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if
case (ip_cloud_csiw)
  cld%n_cloud_type = 4
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      i_cloud_type(i) = ip_cloud_type_sw
    case (ip_clcmp_st_ice)
      i_cloud_type(i) = ip_cloud_type_si
    case (ip_clcmp_cnv_water)
      i_cloud_type(i) = ip_cloud_type_cw
    case (ip_clcmp_cnv_ice)
      i_cloud_type(i) = ip_cloud_type_ci
    end select
  end do
  if ((present(liq_frac).or.present(liq_frac_1d)) .and. &
      (present(ice_frac).or.present(ice_frac_1d)) .and. &
      (present(cloud_frac).or.present(cloud_frac_1d))) then
    call set_cld_field(frac_liq, liq_frac, liq_frac_1d)
    call set_cld_field(frac_ice, ice_frac, ice_frac_1d)
    call set_cld_field(frac, cloud_frac, cloud_frac_1d)
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        if (frac_liq(l, k) + frac_ice(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          cld%frac_cloud(l, k, ip_cloud_type_sw) = &
            frac(l, k)*frac_liq(l, k) / (frac_liq(l, k)+frac_ice(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_si) = &
            frac(l, k)*frac_ice(l, k) / (frac_liq(l, k)+frac_ice(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_sw) = 0.0_RealK
          cld%frac_cloud(l, k, ip_cloud_type_si) = 0.0_RealK
        end if
      end do
    end do
  else if ((present(liq_frac).or.present(liq_frac_1d)) .and. &
           (present(ice_frac).or.present(ice_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_sw), &
                       liq_frac, liq_frac_1d)
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_si), &
                       ice_frac, ice_frac_1d)
  else
    cmessage = 'Liquid and ice cloud fractions not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if
  if ((present(liq_conv_frac).or.present(liq_conv_frac_1d)) .and. &
      (present(ice_conv_frac).or.present(ice_conv_frac_1d)) .and. &
      (present(conv_frac).or.present(conv_frac_1d))) then
    call set_cld_field(frac_liq, liq_conv_frac, liq_conv_frac_1d)
    call set_cld_field(frac_ice, ice_conv_frac, ice_conv_frac_1d)
    call set_cld_field(frac, conv_frac, conv_frac_1d)
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        if (frac_liq(l, k) + frac_ice(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          cld%frac_cloud(l, k, ip_cloud_type_cw) = &
            frac(l, k)*frac_liq(l, k) / (frac_liq(l, k)+frac_ice(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_ci) = &
            frac(l, k)*frac_ice(l, k) / (frac_liq(l, k)+frac_ice(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_cw) = 0.0_RealK
          cld%frac_cloud(l, k, ip_cloud_type_ci) = 0.0_RealK
        end if
      end do
    end do
  else if ((present(liq_conv_frac).or.present(liq_conv_frac_1d)) .and. &
           (present(ice_conv_frac).or.present(ice_conv_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_cw), &
                       liq_conv_frac, liq_conv_frac_1d)
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_ci), &
                       ice_conv_frac, ice_conv_frac_1d)
  else
    cmessage = 'Liquid and ice convective cloud fractions not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if
end select

! Convert mass mixing ratios to in-cloud values
do i = 1, cld%n_condensed
  do k = dimen%id_cloud_top, atm%n_layer
    do l = 1, atm%n_profile
      cld%condensed_mix_ratio(l, k, i) = cld%condensed_mix_ratio(l, k, i) &
        / max(cld%frac_cloud(l, k, i_cloud_type(i)), eps)
    end do
  end do
end do

! Normalise the cloud fractions
do k = dimen%id_cloud_top, atm%n_layer
  do l = 1, atm%n_profile
    cld%w_cloud(l, k) = sum(cld%frac_cloud(l, k, 1:cld%n_cloud_type))
    if (cld%w_cloud(l, k) > min_cloud_fraction) then
      do j=1, cld%n_cloud_type
        cld%frac_cloud(l, k, j) = cld%frac_cloud(l, k, j) / cld%w_cloud(l, k)
      end do
    else
      cld%w_cloud(l, k) = 0.0_RealK
      cld%frac_cloud(l, k, 1:cld%n_cloud_type) = 0.0_RealK
    end if
    if (cld%w_cloud(l, k) > 1.0_RealK + min_cloud_fraction) then
      cmessage = 'Cloud fraction greater than 1.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (cld%w_cloud(l, k) > 1.0_RealK) then
      cld%w_cloud(l, k) = 1.0_RealK
    end if
  end do
end do


contains

  subroutine set_cld_field(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, dimen%id_cloud_top:)
!     Output field
    real(RealK), intent(in), optional :: full_field(:, :)
!     Full field variable
    real(RealK), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    if (present(full_field)) then
      if (l_inv) then
        do k = dimen%id_cloud_top, atm%n_layer
          do l=1, atm%n_profile
            out_field(l, k) = full_field(l, atm%n_layer+1-k)
          end do
        end do
      else
        do k = dimen%id_cloud_top, atm%n_layer
          do l=1, atm%n_profile
            out_field(l, k) = full_field(l, k)
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_inv) then
        do k = dimen%id_cloud_top, atm%n_layer
          do l=1, atm%n_profile
            out_field(l, k) = oned_field(atm%n_layer+1-k)
          end do
        end do
      else
        do k = dimen%id_cloud_top, atm%n_layer
          do l=1, atm%n_profile
            out_field(l, k) = oned_field(k)
          end do
        end do
      end if
    else
      cmessage = 'The required cloud fields have not been provided.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)      
    end if
  end subroutine set_cld_field


  subroutine set_ice_dim()

    use rad_pcf, only: &
      ip_ice_adt, ip_ice_agg_de, ip_ice_agg_de_sun
    implicit none

    ! Parameters for the aggregate parametrization.
    real (RealK), parameter :: a0_agg_cold = 7.5094588e-04_RealK
    real (RealK), parameter :: b0_agg_cold = 5.0830326e-07_RealK
    real (RealK), parameter :: a0_agg_warm = 1.3505403e-04_RealK
    real (RealK), parameter :: b0_agg_warm = 2.6517429e-05_RealK
    real (RealK), parameter :: t_switch    = 216.208_RealK
    real (RealK), parameter :: t0_agg      = 279.5_RealK
    real (RealK), parameter :: s0_agg      = 0.05_RealK
    real (RealK), parameter :: dge2de      = &
      (3.0_RealK/2.0_RealK)*(3.0_RealK/(2.0_RealK*sqrt(3.0_RealK)))

    select case (cld%i_condensed_param(i))

    case (ip_ice_adt)
      ! This parametrization is based on the mean maximum
      ! dimension of the crystal, determined as a function of
      ! the local temperature. The size is limited to its value
      ! at the freezing level.
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          cld%condensed_dim_char(l, k, i) = min( 7.198755e-04_RealK, &
            exp(5.522e-02_RealK * (atm%t(l, k) - 2.7965e+02_RealK)) &
            / 9.702e+02_RealK )
        end do
      end do

    case (ip_ice_agg_de, ip_ice_agg_de_sun)
      ! Aggregate parametrization based on effective dimension.
      ! The fit provided here is based on Stephan Havemann's fit of
      ! Dge with temperature, consistent with David Mitchell's treatment
      ! of the variation of the size distribution with temperature. The
      ! parametrization of the optical properties is based on De
      ! (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
      ! (=(2*SQRT(3)/3)*volume/projected area), which explains the
      ! conversion factor. The fit to Dge is in two sections, because
      ! Mitchell's relationship predicts a cusp at 216.208 K. Limits
      ! of 8 and 124 microns are imposed on Dge: these are based on this
      ! relationship and should be reviewed if it is changed. Note also
      ! that the relationship given here is for polycrystals only.
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          ! Preliminary calculation of Dge.
          if (atm%t(l, k) < t_switch) then
            cld%condensed_dim_char(l, k, i) = &
              a0_agg_cold*exp(s0_agg*(atm%t(l, k)-t0_agg))+b0_agg_cold
          else
            cld%condensed_dim_char(l, k, i) = &
              a0_agg_warm*exp(s0_agg*(atm%t(l, k)-t0_agg))+b0_agg_warm
          end if
          ! Limit and convert to De.
          cld%condensed_dim_char(l, k, i) = dge2de * min( 1.24e-04_RealK, &
            max(8.0e-06_RealK, cld%condensed_dim_char(l, k, i)) )
        end do
      end do

      case default
        ! A default value of 30-microns is assumed.
        cld%condensed_dim_char(:, :, i) = 30.0e-6_RealK

    end select

  end subroutine set_ice_dim

end subroutine set_cld
end module socrates_set_cld
