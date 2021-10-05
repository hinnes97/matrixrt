! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the MCICA variables in the Socrates cloud type
!
!------------------------------------------------------------------------------
module socrates_set_cld_mcica
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_CLD_MCICA'
contains

subroutine set_cld_mcica(cld, mcica_data, control, dimen, spectrum, atm, &
  rand_seed, cloud_horizontal_rsd)

use def_cld,      only: StrCld, allocate_cld_mcica
use def_mcica,    only: StrMcica, ip_mcica_full_sampling, &
  ip_mcica_single_sampling, ip_mcica_optimal_sampling
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use def_atm,      only: StrAtm
use realtype_rd,  only: RealK
use rad_pcf,      only: &
  ip_cloud_off, ip_cloud_mcica, ip_solar, &
  i_normal, i_err_fatal
use ereport_mod,  only: ereport
use errormessagelength_mod, only: errormessagelength

use socrates_cloud_gen, only: cloud_gen

implicit none


! Cloud properties:
type(StrCld),      intent(inout) :: cld

! Mcica data:
type(StrMcica),    intent(in) :: mcica_data

! Control options:
type(StrCtrl),     intent(in) :: control

! Dimensions:
type(StrDim),      intent(in) :: dimen

! Spectral data:
type(StrSpecData), intent(in) :: spectrum

! Atmospheric properties:
type(StrAtm),      intent(in) :: atm

integer,     intent(in), optional :: rand_seed(:)
!   Random seed for cloud generator
real(RealK), intent(in), optional :: cloud_horizontal_rsd
!   Relative standard deviation of sub-grid cloud condensate

real(RealK), dimension(dimen%nd_profile, dimen%id_cloud_top:dimen%nd_layer) :: &
  c_cloud, &
!   Amount of convective cloud
  c_ratio, &
!   Ratio of convective cloud condensate to mean condensate
  ls_ratio, &
!   Ratio of large-scale cloud condensate to mean condensate
  dp_corr_cloud, &
!   Cloud fraction decorrelation length
  dp_corr_cond, &
!   Cloud condensate decorrelation length
  cond_rsd
!   Relative standard deviation of sub-grid cloud condensate
integer :: rnd_seed(dimen%nd_profile)
!   Random seed

integer :: l, i, j, k
integer :: i_k, i_band, n_k
integer :: n_subcol_fill

integer                      :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_CLD_MCICA'
character (len=errormessagelength) :: cmessage


if (control%i_cloud_representation /= ip_cloud_off .and. &
    control%i_cloud == ip_cloud_mcica) then

  ! Allocate MCICA data arrays
  call allocate_cld_mcica(cld, dimen, spectrum)

  ! Set the number of sub-columns to be sampled by each k-term
  select case (control%i_mcica_sampling)
  case (ip_mcica_full_sampling)
    cld%subcol_k = mcica_data%n_subcol_gen
  case (ip_mcica_single_sampling)
    cld%subcol_k = 1
  case (ip_mcica_optimal_sampling)
    if (control%isolir == ip_solar) then
      do i_k=1, spectrum%dim%nd_k_term
        do i_band=1, spectrum%basic%n_band
          cld%subcol_k(i_band, i_k) = mcica_data%sw_subcol_k(i_band, i_k)
        end do
      end do
    else
      do i_k=1, spectrum%dim%nd_k_term
        do i_band=1, spectrum%basic%n_band
          cld%subcol_k(i_band, i_k) = mcica_data%lw_subcol_k(i_band, i_k)
        end do
      end do
    end if
  case default
    cmessage = 'The value of i_mcica_sampling is not valid.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)           
  end select

  ! Set the first sub-column to be sampled by each k-term
  cld%first_subcol_k(control%first_band, 1) = 1
  do i_band = control%first_band, control%last_band
    n_k = spectrum%gas%i_band_k(i_band, spectrum%gas%index_absorb(1, i_band))
    do i_k = 1, n_k
      cld%first_subcol_k(i_band, i_k+1) &
        = cld%first_subcol_k(i_band, i_k) + cld%subcol_k(i_band, i_k)
    end do
    if (i_band < control%last_band) &
      cld%first_subcol_k(i_band+1, 1) = cld%first_subcol_k(i_band, n_k+1)
  end do

  ! Set order of sub-columns to balance SW and LW effect
  if (control%isolir == ip_solar) then
    do l=1, dimen%nd_subcol_req
      cld%subcol_reorder(l) = mod(l, mcica_data%n_subcol_gen)
    end do
  else
    select case (control%i_mcica_sampling)
    case (ip_mcica_full_sampling)
      do l=1, dimen%nd_subcol_req
        cld%subcol_reorder(l) = l
      end do
    case (ip_mcica_single_sampling)
      do l=1, dimen%nd_subcol_req
        cld%subcol_reorder(l) = mod(mcica_data%lw_subcol_reorder_single(l), &
                                    mcica_data%n_subcol_gen)
      end do
    case (ip_mcica_optimal_sampling)
      do l=1, dimen%nd_subcol_req
        cld%subcol_reorder(l) = mod(mcica_data%lw_subcol_reorder_optimal(l), &
                                    mcica_data%n_subcol_gen)
      end do
    end select
  end if
  where (cld%subcol_reorder == 0) cld%subcol_reorder=mcica_data%n_subcol_gen

  if (present(rand_seed)) then
    rnd_seed(1:atm%n_profile) = rand_seed(1:atm%n_profile)
  else
    do i=1, atm%n_profile
      rnd_seed(i) = 10 + i
    end do
  end if
  if (present(cloud_horizontal_rsd)) then
    cond_rsd(:, :) = cloud_horizontal_rsd
  else
    cond_rsd(:, :) = 0.75_RealK
  end if

  ! Currently a single value of the decorrelation scale is used
  dp_corr_cloud(:, :) = cld%dp_corr_strat
  dp_corr_cond(:, :) = dp_corr_cloud(:, :) * 0.5_RealK

  ! Convective cloud amounts not yet implemented
  c_cloud(:, :) = 0.0_RealK
  c_ratio(:, :) = 0.0_RealK
  ls_ratio(:, :) = 1.0_RealK

  ! Call the cloud generator
  call cloud_gen(dimen%nd_layer, dimen%id_cloud_top, atm%n_layer, &
    dimen%nd_profile, 1, atm%n_profile, &
    mcica_data%n_subcol_gen, mcica_data%n1, mcica_data%n2, &
    mcica_data%ipph, control%i_overlap, rnd_seed, &
    dp_corr_cloud, dp_corr_cond, &
    cond_rsd, cld%w_cloud, c_cloud, c_ratio, ls_ratio, atm%p, &
    mcica_data%xcw, &
    cld%n_subcol_cld, cld%c_sub(:,:,:,1))
  
  ! All cloud types assumed to have the same sub-grid distribution:
  do i=2, cld%n_cloud_type
    cld%c_sub(:,:,:,i) = cld%c_sub(:,:,:,1)
  end do

  select case (control%i_mcica_sampling)
  case (ip_mcica_full_sampling)
    ! In this case we treat the clear sub-columns as cloudy sub-columns
    ! so the clear-sky fraction is implicit in the summing of the 
    ! "cloudy" sub-columns
    cld%frac_cloudy = 1.0_RealK
  case default
    ! Otherwise, where there are less cloudy subcolumns than required,
    ! the cloudy sub-columns are copied to all sampled sub-columns.
    n_subcol_fill = min(dimen%nd_subcol_req, mcica_data%n_subcol_gen)
    do i=1, atm%n_profile
      if (cld%n_subcol_cld(i) < n_subcol_fill .and. &
          cld%n_subcol_cld(i) > 0) then
        do j=cld%n_subcol_cld(i)+1, n_subcol_fill
          do k=dimen%id_cloud_top, atm%n_layer
            l=j-cld%n_subcol_cld(i)
            cld%c_sub(i,k,j,:)=cld%c_sub(i,k,l,1)
          end do
        end do
      end if
    end do
    ! The total cloud cover is set to the correct fraction.
    cld%frac_cloudy = real(cld%n_subcol_cld, RealK) &
                    / real(mcica_data%n_subcol_gen, RealK)
  end select

end if

end subroutine set_cld_mcica
end module socrates_set_cld_mcica
