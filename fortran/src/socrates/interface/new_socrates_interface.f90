MODULE socrates_interface_mod

  ! Socrates calculation interface module
  ! Takes FMS time, spectra, temperature, and pressure
  ! Outputs FMS heating rate, and downwards surface LW and SW
  ! MDH 30/01/18

  ! Socrates interface significantly updated in order to allow for
  ! seasonally-varying insolation, time-varying CO2, setting of radiation
  ! parameters using namelists, specified Ozone, radiation timestep!=dt_atmos,
  ! and correcting errors. Some of the additional code was adapted from Isca's
  ! RRTM_RADIATION, which was originally written by Martin Jucker and adapted by
  ! the Isca development team.
  ! SIT 06/08/18


  !----------

!#ifdef INTERNAL_FILE_NML
!  use mpp_mod, only: input_nml_file
!#else
  ! Socrates modules
  USE read_control_mod
  USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
  USE def_spectrum
  !USE constants_mod, only: grav, rdgas, rvgas, cp_air
  USE soc_constants_mod  
  USE socrates_config_mod
  USE params, ONLY: dp

  IMPLICIT NONE


  ! Input spectra
  TYPE (StrSpecData) :: spectrum_calc
  TYPE (StrSpecData) :: spectrum_lw, spectrum_lw_hires
  TYPE (StrSpecData) :: spectrum_sw, spectrum_sw_hires

  ! Control options:
  TYPE(StrCtrl) :: control_calc
  TYPE(StrCtrl) :: control_sw, control_sw_hires
  TYPE(StrCtrl) :: control_lw, control_lw_hires

  INTEGER :: n_soc_bands_lw, n_soc_bands_sw
  INTEGER :: n_soc_bands_lw_hires, n_soc_bands_sw_hires


  REAL :: dt_last !Time of last radiation calculation - used to tell whether it is time to recompute radiation or not
  REAL(r_def), allocatable, dimension(:,:,:) :: tdt_soc_sw_store, tdt_soc_lw_store, thd_sw_flux_down_store, thd_lw_flux_up_store
  REAL(r_def), allocatable, dimension(:,:) :: net_surf_sw_down_store, surf_lw_down_store, coszen_store
  REAL(r_def), allocatable, dimension(:) :: outputted_soc_spectral_olr, spectral_olr_store
  REAL(r_def), allocatable, dimension(:)     :: soc_bins_lw, soc_bins_sw


!ExoFMS namelist
!  real :: stellar_constant
!  character :: lw_spectral_filename
!  character :: sw_spectral_filename
  real :: grav,rdgas,rvgas,cp_air

!  namelist /socrates_rad_nml/ stellar_constant, lw_spectral_filename, sw_spectral_filename

CONTAINS

  SUBROUTINE socrates_init()
    !! Initialises Socrates spectra, arrays, and constants

 !   USE astronomy_mod, only: astronomy_init
    !    USE interpolator_mod, only: interpolate_type, interpolator_init, ZERO
    use socrates_config_mod
    use rad_ccf, only: grav_acc, r_gas_dry, cp_air_dry

    ! Arguments
    !INTEGER, INTENT(in), DIMENSION(4) :: axes
    !! NB axes refers to the handles of the axes defined in fv_diagnostics
    !TYPE(time_type), INTENT(in)       :: Time, delta_t_atmos
    !INTEGER, INTENT(in)               :: num_levels
    !REAL, INTENT(in) , DIMENSION(:,:)   :: lat
    !REAL, INTENT(in) , DIMENSION(:,:)   :: lonb, latb
        
    integer :: io, stdlog_unit, nml_unit
    integer :: res, time_step_seconds
    real    :: day_in_s_check
    logical :: file_exist
    character(len=60) :: filename = "input.nml"
    !-------------------------------------------------------------------------------------
    
!#ifdef INTERNAL_FILE_NML
!   read (input_nml_file, nml=socrates_rad_nml, iostat=io)
    !#else
    grav=12.43!grav_acc
    rdgas=3779!r_gas_dry
   rvgas=3779.*2.2/18.
    cp_air=rdgas*7./2.!cp_air_dry

  inquire(file='input.nml', exist=file_exist)
  if ( file_exist ) then
     write(*,*) 'FILE EXISTS'
     open(99, file='input.nml', iostat=io)
     if (io .ne. 0) write(*,*) 'OPEN ERROR: IOSTAT=',io
     rewind(99)
     read (99, socrates_rad_nml, iostat=io)
     if (io .ne. 0) write(*,*) 'READ ERROR: IOSTAT=',io
     close(99)
     if (io .ne. 0) write(*,*) 'CLOSE ERROR: IOSTAT=',io
   endif
   write(*,*) 'STELLAR CONSTANT', stellar_constant
!#endif
!stdlog_unit = stdlog()
!write(stdlog_unit, socrates_rad_nml)

      if (lw_spectral_filename .eq. 'unset') then
       stop 'lw_spectral_filename is unset, and must point to a valid spectral file'
      endif

      if (sw_spectral_filename .eq. 'unset') then
       stop 'sw_spectral_filename is unset, and must point to a valid spectral file'
      endif

      if (lw_hires_spectral_filename .eq. 'unset') then
         write(*,*) 'WARNING: lw_hires_spectral_filename is unset, making equal to lw_spectral_filename'
        lw_hires_spectral_filename = lw_spectral_filename
    endif

    if (sw_hires_spectral_filename .eq. 'unset') then
       write(*,*) 'WARNING: sw_hires_spectral_filename is unset, making equal to sw_spectral_filename'
           sw_hires_spectral_filename = sw_spectral_filename
     endif
     
     write(*,*) 'IN INIT FUNCTION'
    ! Socrates spectral files -- should be set by namelist
    control_lw%spectral_file = lw_spectral_filename
    control_lw_hires%spectral_file = lw_hires_spectral_filename

    control_sw%spectral_file = sw_spectral_filename
    control_sw_hires%spectral_file = sw_hires_spectral_filename

    ! Read in spectral files
    CALL read_spectrum(control_lw%spectral_file,spectrum_lw)
    CALL read_spectrum(control_lw_hires%spectral_file,spectrum_lw_hires)
    CALL read_spectrum(control_sw%spectral_file,spectrum_sw)
    CALL read_spectrum(control_sw_hires%spectral_file,spectrum_sw_hires)

    ! Set Socrates configuration
    CALL read_control(control_lw,spectrum_lw)
    CALL read_control(control_lw_hires,spectrum_lw_hires)
    CALL read_control(control_sw,spectrum_sw)
    CALL read_control(control_sw_hires,spectrum_sw_hires)

    ! Specify LW and SW setups
    control_sw%isolir=1
    control_sw_hires%isolir=1
    control_lw%isolir=2
    control_lw_hires%isolir=2

! New hires part
if(socrates_hires_mode) then
        allocate(soc_bins_lw(spectrum_lw_hires%dim%nd_band))
        allocate(soc_bins_sw(spectrum_sw_hires%dim%nd_band))    
        soc_bins_lw = spectrum_lw_hires%Basic%wavelength_long            
        soc_bins_sw = spectrum_sw_hires%Basic%wavelength_short       
    else
        allocate(soc_bins_lw(spectrum_lw%dim%nd_band))
        allocate(soc_bins_sw(spectrum_sw%dim%nd_band))    
        soc_bins_lw = spectrum_lw%Basic%wavelength_long            
        soc_bins_sw = spectrum_sw%Basic%wavelength_short              
    endif    
    
    !Need to actually give bins arrays values    
    
!    id_soc_bins_lw = diag_axis_init('soc_bins_lw', soc_bins_lw, 'cm^-1', 'n', 'socrates lw spectral bin centers', set_name='socrates_lw_bins')
!    id_soc_bins_sw = diag_axis_init('soc_bins_sw', soc_bins_sw, 'cm^-1', 'n', 'socrates sw spectral bin centers', set_name='socrates_sw_bins')

! End of new hires part



    ! Number of bands
    n_soc_bands_lw = spectrum_lw%dim%nd_band
    n_soc_bands_lw_hires = spectrum_lw_hires%dim%nd_band
    n_soc_bands_sw = spectrum_sw%dim%nd_band
    n_soc_bands_sw_hires = spectrum_sw_hires%dim%nd_band

    if (socrates_hires_mode .eqv. .True.) then
        allocate(outputted_soc_spectral_olr(n_soc_bands_lw_hires))
    else
        allocate(outputted_soc_spectral_olr(n_soc_bands_lw ))
    endif

! new hires
!        if (id_soc_spectral_olr > 0) then 
!            if (socrates_hires_mode == .True.) then
!                allocate(spectral_olr_store(b,1)-1, size(latb,2)-1, n_soc_bands_lw_hires))
!            else
!                allocate(spectral_olr_store(size(lonb,1)-1, size(latb,2)-1, n_soc_bands_lw ))
!            endif
!        endif 

    ! Print Socrates init data from one processor
       PRINT*, 'Initialised Socrates v17.03'
       PRINT*, 'Stellar constant = ', stellar_constant
       PRINT*, 'Longwave spectral file = ', TRIM(control_lw%spectral_file), ' WITH ', n_soc_bands_lw, ' bands'
       PRINT*, 'Longwave hires spectral file = ', TRIM(control_lw_hires%spectral_file), ' WITH ', n_soc_bands_lw_hires, ' bands'
       PRINT*, 'Shortwave spectral file = ', TRIM(control_sw%spectral_file), ' WITH ', n_soc_bands_sw, ' bands'
       PRINT*, 'Shortwave hires spectral file = ', TRIM(control_sw_hires%spectral_file), ' WITH ', n_soc_bands_sw_hires, ' bands'
       PRINT*, ' '
       PRINT*, '-----------------------------------'
       PRINT*, ' '

    return
  end subroutine socrates_init
  ! ==================================================================================


  ! Set up the call to the Socrates radiation scheme
  ! -----------------------------------------------------------------------------
  subroutine socrates_interface(rlat, rlon, soc_lw_mode,  &
       fms_temp, fms_spec_hum, fms_ozone, fms_co2, fms_h2, fms_t_surf, fms_p_full, fms_p_half, &
       fms_z_full, fms_z_half, fms_albedo, fms_coszen, n_profile, n_layer,        &
       output_heating_rate, output_flux_down, output_flux_up,  output_soc_spectral_olr,&
       output_flux_direct, t_half_level_out )

    use realtype_rd
    use read_control_mod
    use socrates_calc_mod
    use compress_spectrum_mod
    use def_spectrum
    use def_dimen,   only: StrDim
    use def_control, only: StrCtrl,  allocate_control,   deallocate_control
    use def_atm,     only: StrAtm,   allocate_atm,       deallocate_atm
    use def_cld,     only: StrCld,   allocate_cld,       deallocate_cld, &
         allocate_cld_prsc,  deallocate_cld_prsc, &
         allocate_cld_mcica, deallocate_cld_mcica
    use def_aer,     only: StrAer,   allocate_aer,       deallocate_aer, &
         allocate_aer_prsc,  deallocate_aer_prsc
    use def_bound,   only: StrBound, allocate_bound,     deallocate_bound
    use def_out,     only: StrOut,                       deallocate_out

    !-----------------------------------------------------------------------
    implicit none


    INTEGER(i_def), intent(in) :: n_profile, n_layer
    logical, intent(in) :: soc_lw_mode
    INTEGER(i_def) :: nlat

    ! Input arrays
    real(r_def), intent(in) :: fms_temp(:), fms_spec_hum(:), fms_ozone(:), fms_co2(:), fms_h2(:)
    real(r_def), intent(in) :: fms_p_full(:)
    real(r_def), intent(in) :: fms_p_half(:)
    real(r_def), intent(in) :: fms_t_surf, fms_albedo
    real(r_def), intent(in) :: fms_coszen
    real(r_def), intent(in) :: rlon
    real(r_def), intent(in) :: rlat
    real(r_def), intent(in) :: fms_z_full(:), fms_z_half(:)


    ! Output arrays
    real(r_def), intent(out) :: output_heating_rate(:)
    real(r_def), intent(out) :: output_flux_up(:)
    real(r_def), intent(out) :: output_flux_down(:)
    real(r_def), intent(out), optional :: output_flux_direct(:)    
    real(r_def), intent(out), optional :: output_soc_spectral_olr(:)
    real(r_def), intent(out), optional :: t_half_level_out(:)

    ! Hi-res output
    INTEGER, PARAMETER :: out_unit=20
    CHARACTER(len=200) :: file_name
    REAL(r_def) :: soc_spectral_olr(n_profile, size(outputted_soc_spectral_olr,1))

    ! Arrays to send to Socrates
    real(r_def), dimension(n_profile, n_layer) :: input_p, input_t, input_mixing_ratio, &
         input_d_mass, input_density, input_layer_heat_capacity, &
         soc_heating_rate, input_o3_mixing_ratio, &
          input_co2_mixing_ratio,z_full_reshaped, input_h2_mixing_ratio
    real(r_def), dimension(n_profile,0:n_layer) :: input_p_level, input_t_level, soc_flux_direct, &
         soc_flux_down, soc_flux_up, z_half_reshaped
    real(r_def), dimension(n_profile) :: input_t_surf, input_cos_zenith_angle, input_solar_irrad, &
         input_orog_corr, input_planet_albedo


    ! Socrates options
    integer(i_def) :: input_n_cloud_layer
    integer(i_def) :: input_n_aer_mode
    integer(i_def) :: input_cld_subcol_gen
    integer(i_def) :: input_cld_subcol_req



    ! Dimensions:
    type(StrDim) :: dimen
    type(StrAtm) :: atm_input

    ! Loop variables
    integer(i_def) :: i, j, k, l, n, lon, si, sj, sk, sb

    ! chunking loop variable
    integer(i_def) :: n_chunk_loop, idx_chunk_start, idx_chunk_end, i_chunk, n_profile_chunk

    !DIAG Diagnostic
    logical :: used


    !----------------------------i
    ! Set array sizes
    input_n_cloud_layer = n_layer
    input_n_aer_mode = 1!n_layer
    input_cld_subcol_gen = n_layer
    input_cld_subcol_req = n_layer
    !si = size(fms_temp,1)
    !sj = size(fms_temp,2)
    sk = size(fms_temp)

    sb=n_soc_bands_lw





      !Set input T, p, p_level, and mixing ratio profiles
          input_t(n_profile,:) = fms_temp!reshape(fms_temp(:,:,:),(/si*sj,sk /))
          input_p(n_profile,:) = fms_p_full!reshape(fms_p_full(:,:,:),(/si*sj,sk /))
          input_p_level(n_profile,:) = fms_p_half!reshape(fms_p_half(:,:,:),(/si*sj,sk+1 /))
          
          if (account_for_effect_of_water .eqv. .true.) then
              input_mixing_ratio(n_profile,:) = fms_spec_hum!reshape(fms_spec_hum(:,:,:) / (1. - fms_spec_hum(:,:,:)),(/si*sj,sk /)) !Mass mixing ratio = q / (1-q)
          else
              input_mixing_ratio = 0.0
          endif
          
          if (account_for_effect_of_ozone .eqv. .true.) then
            input_o3_mixing_ratio(n_profile,:) = fms_ozone!reshape(fms_ozone(:),(/n_profile,sk /))
          else         
            input_o3_mixing_ratio = 0.0
          endif

          input_co2_mixing_ratio(n_profile,:) = fms_co2!reshape(fms_co2(:),(/n_profile,sk/))
          input_h2_mixing_ratio(n_profile,:) = fms_h2

          !-------------

          !Default parameters
          input_cos_zenith_angle(n_profile) = fms_coszen!reshape((fms_coszen(:,:)),(/si*sj /))
          input_orog_corr = 0.0
          input_planet_albedo(n_profile) = fms_albedo!reshape(fms_albedo(:,:),(/n_profile /))

          !Set tide-locked flux - should be set by namelist eventually!
          input_solar_irrad(n_profile) = stellar_constant
          input_t_surf(n_profile) = fms_t_surf!reshape(fms_t_surf(:,:),(/si*sj /))
          z_full_reshaped(n_profile,:) = fms_z_full!reshape(fms_z_full(:,:,:), (/si*sj, sk/))
          z_half_reshaped(n_profile,:) = fms_z_half!reshape(fms_z_half(:,:,:), (/si*sj, sk+1/))

          !--------------
          !Set input t_level by scaling t
       if (use_pressure_interp_for_half_levels) then
             DO k = 1,n_layer-1
                input_t_level(:,k) = input_t(:,k) + (input_t(:,k+1)-input_t(:,k)) * &
                     ((input_p_level(:,k)-input_p(:,k))/(input_p(:,k+1)-input_p(:,k)))
             END DO
!              input_t_level(:,n_layer) = input_t(:,n_layer) + input_t(:,n_layer) - input_t_level(:,n_layer-1)
             input_t_level(:,n_layer) = input_t(:,n_layer) + (input_t(:,n_layer)-input_t(:,n_layer-1)) * &
                  ((input_p_level(:,n_layer)-input_p(:,n_layer))/(input_p(:,n_layer)-input_p(:,n_layer-1)))
             
!              input_t_level(:,0) = input_t(:,1) - (input_t_level(:,1) - input_t(:,1))
             input_t_level(:,0) = input_t(:,1) + (input_t(:,2)-input_t(:,1)) * &
                  ((input_p_level(:,0)-input_p(:,1))/(input_p(:,2)-input_p(:,1)))
             
       else

          call interp_temp(z_full_reshaped,z_half_reshaped,input_t, input_t_level)

       endif

         if (present(t_half_level_out)) then
             t_half_level_out(:) = input_t_level(n_profile,:)!reshape(input_t_level,(/si,sj,sk+1 /))
         endif

          !Set input dry mass, density, and heat capacity profiles
          DO i=n_layer, 1, -1
             input_d_mass(:,i) = (input_p_level(:,i)-input_p_level(:,i-1))/grav
             input_density(:,i) = input_p(:,i)/(rdgas*input_t(:,i))
             input_layer_heat_capacity(:,i) = input_d_mass(:,i)*cp_air
          END DO


          ! Zero heating rate
          soc_heating_rate = 0.0

          ! Test if LW or SW mode
          if (soc_lw_mode .eqv. .TRUE.) then
             control_lw%isolir = 2
             CALL read_control(control_lw, spectrum_lw)
             if (socrates_hires_mode .eqv. .FALSE.) then
                control_calc = control_lw
                spectrum_calc = spectrum_lw
             else
                control_calc = control_lw_hires
                spectrum_calc = spectrum_lw_hires
             end if

          else
             control_sw%isolir = 1
             CALL read_control(control_sw, spectrum_sw)
             if(socrates_hires_mode .eqv. .FALSE.) then
                control_calc = control_sw
                spectrum_calc = spectrum_sw
             else
                control_calc = control_sw_hires
                spectrum_calc = spectrum_sw_hires
             end if

          end if


          ! Do calculation
          CALL read_control(control_calc, spectrum_calc)

          if (soc_lw_mode .eqv. .TRUE.) then
            CALL socrates_calc(control_calc, spectrum_calc,                      &
               n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,             &
               input_cld_subcol_gen, input_cld_subcol_req,                                  &
               input_p,                                    & 
               input_t,                                    &
               input_t_level,                              & 
               input_d_mass,                               &
               input_density,                              &
               input_mixing_ratio,                         &
               input_o3_mixing_ratio,                      &
               input_co2_mixing_ratio,                     &
               input_h2_mixing_ratio,                      &
               input_t_surf,                                 &
               input_cos_zenith_angle,                       &
               input_solar_irrad,                            &
               input_orog_corr,                              &
               l_planet_grey_surface,                                                       &
               input_planet_albedo,                          &
               input_planet_emissivity,                                                     &
               input_layer_heat_capacity,                  &
               soc_flux_direct,                            &
               soc_flux_down,                              &
               soc_flux_up,                                &
               soc_heating_rate,                           &
               soc_spectral_olr)

          else
            CALL socrates_calc(control_calc, spectrum_calc,                      &
               n_profile, n_layer, input_n_cloud_layer, input_n_aer_mode,             &
               input_cld_subcol_gen, input_cld_subcol_req,                                  &
               input_p,                                    &
               input_t,                                    &
               input_t_level,                              &
               input_d_mass,                               &
               input_density,                              &
               input_mixing_ratio,                         &
               input_o3_mixing_ratio,                      &
               input_co2_mixing_ratio,                     &
               input_h2_mixing_ratio,                      &
               input_t_surf,                                 &
               input_cos_zenith_angle,                       &
               input_solar_irrad,                            &
               input_orog_corr,                              &
               l_planet_grey_surface,                                                       &
               input_planet_albedo,                          &
               input_planet_emissivity,                                                     &
               input_layer_heat_capacity,                  &
               soc_flux_direct,                            &
               soc_flux_down,                           &
               soc_flux_up,                             &
               soc_heating_rate)
          endif


          ! Set output arrays
          output_flux_up = soc_flux_up(n_profile,:)
          output_flux_down = soc_flux_down(n_profile,:)!reshape(soc_flux_down(:,:),(/si,sj,sk+1 /))          

          if(present(output_flux_direct)) then
              output_flux_direct = soc_flux_direct(n_profile,:)!reshape(soc_flux_direct(:,:),(/si,sj,sk+1 /))          
          endif
          
          output_heating_rate = soc_heating_rate(n_profile,:)!reshape(soc_heating_rate(:,:),(/si,sj,sk /))

          if (soc_lw_mode .eqv. .TRUE.) then
              output_soc_spectral_olr = soc_spectral_olr(n_profile,:) !reshape(soc_spectral_olr(:,:),(/si,sj,int(n_soc_bands_lw,i_def) /))
          endif

  end subroutine socrates_interface

subroutine run_socrates(rad_lat, rad_lon, temp_in, q_in, h2_in, t_surf_in, p_full_in, p_half_in, z_full_in, z_half_in, albedo_in, &
       temp_tend, net_surf_sw_down, surf_lw_down, net_flux)  

!    use astronomy_mod, only: diurnal_solar
    !use constants_mod,         only: pi, wtmco2, wtmozone, rdgas, gas_constant
  !    use interpolator_mod,only: interpolator

    use socrates_config_mod    

    real(dp), intent(in)       :: t_surf_in, albedo_in
    real(dp), intent(in), dimension(:)   :: temp_in, p_full_in, q_in, z_full_in, h2_in
    real(dp), intent(in), dimension(:)  :: p_half_in, z_half_in
    real(dp), intent(inout), dimension(:) :: temp_tend
    real(dp), intent(out)   :: net_surf_sw_down, surf_lw_down
    !real, intent(in) :: delta_t
    real(dp), intent(in) :: rad_lat, rad_lon
    real(dp), intent(out) :: net_flux(:)

    integer(i_def) :: n_profile, n_layer

    real(r_def) :: t_surf_for_soc, rad_lat_soc, rad_lon_soc, albedo_soc, coszen_soc
    real(r_def), dimension(size(temp_in)) :: tg_tmp_soc, q_soc, ozone_soc, co2_soc, h2_soc, &
    p_full_soc, output_heating_rate_sw, output_heating_rate_lw, output_heating_rate_total, z_full_soc
    real(r_def), dimension(size(temp_in)+1) :: p_half_soc, t_half_out, z_half_soc, output_soc_flux_sw_down, output_soc_flux_sw_up, &
         output_soc_flux_lw_down, output_soc_flux_lw_up

    logical :: soc_lw_mode, used
    integer :: seconds, days, year_in_s
    real(dp) :: r_seconds, r_days, r_total_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, &
         dt_rad_radians, day_in_s, r_solday, r_dt_rad_avg
    real(dp) :: coszen, fracsun   
    real(dp), dimension(size(temp_in)) :: ozone_in, co2_in
    real(dp), dimension(size(temp_in)) :: thd_sw_flux_down, thd_lw_flux_up
    !type(time_type) :: Time_loc

    real :: pi, wtmco2, wtmozone, rdgas, gas_constant
    pi = 4._dp*atan(1._dp)
    wtmco2 = 44._dp
    wtmozone=48._dp
    rdgas=3779._dp
    gas_constant=8314._dp

       !make sure we run perpetual when solday > 0)
 !         if(solday > 0)then
 !            Time_loc = set_time(seconds,solday)
 !         else
 !            Time_loc = Time
 !         endif

       !Set tide-locked flux if tidally-locked = .true. Else use diurnal-solar
       !to calculate insolation from orbit!
       if (tidally_locked .eqv. .TRUE.) then
           coszen = COS(rad_lat)*COS(rad_lon)
           if (coszen < 0.0_dp) coszen = 0.0_dp
           
       endif

       ozone_in = 0.00_dp
       co2_in = 0.5_dp

       
      !get ozone 
!       if(do_read_ozone)then
!         call interpolator( o3_interp, Time_diag, p_half_in, ozone_in, trim(ozone_field_name))
!         if (input_o3_file_is_mmr==.false.) then
!             ozone_in = ozone_in * wtmozone / (1000. * gas_constant / rdgas ) !Socrates expects all abundances to be mass mixing ratio. So if input file is volume mixing ratio, it must be converted to mass mixing ratio using the molar masses of dry air and ozone
             ! Molar mass of dry air calculated from gas_constant / rdgas, and converted into g/mol from kg/mol by multiplying by 1000. This conversion is necessary because wtmozone is in g/mol.
             
!         endif 
!       endif

!       if (input_co2_mmr==.false.) then
!           co2_in = co2_ppmv * 1.e-6 * wtmco2 / (1000. * gas_constant / rdgas )!Convert co2_ppmv to a mass mixing ratio, as required by socrates
             ! Molar mass of dry air calculated from gas_constant / rdgas, and converted into g/mol from kg/mol by multiplying by 1000. This conversion is necessary because wtmco2 is in g/mol.           
!       else
!           co2_in = co2_ppmv * 1.e-6 !No need to convert if it is already a mmr
!       endif
      
       !get co2 
!       if(do_read_co2)then
!         call interpolator( co2_interp, Time_diag, p_half_in, co2_in, trim(co2_field_name))
!         if (input_co2_mmr==.false.) then
!             co2_in = co2_in * 1.e-6 * wtmco2 / (1000. * gas_constant / rdgas )
             ! Molar mass of dry air calculated from gas_constant / rdgas, and converted into g/mol from kg/mol by multiplying by 1000. This conversion is necessary because wtmco2 is in g/mol.
             
!         endif
!       endif

       n_profile = INT(1, kind(i_def))!INT(size(temp_in,2)*size(temp_in,1), kind(i_def))
       n_layer   = INT(size(temp_in), kind(i_def))
       t_surf_for_soc = REAL(t_surf_in, kind(r_def))
       
       ! LW calculation
       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes
       soc_lw_mode = .TRUE.
       rad_lat_soc = REAL(rad_lat, kind(r_def))
       rad_lon_soc = REAL(rad_lon, kind(r_def))
       tg_tmp_soc =  REAL(temp_in, kind(r_def))
       q_soc      =  REAL(q_in, kind(r_def))
       h2_soc     = REAL(h2_in, kind(r_def))
       ozone_soc  =  REAL(ozone_in, kind(r_def)) 
       co2_soc    =  REAL(co2_in, kind(r_def))      
       p_full_soc = REAL(p_full_in, kind(r_def))
       p_half_soc = REAL(p_half_in, kind(r_def))
       albedo_soc = REAL(albedo_in, kind(r_def))
       z_full_soc = REAL(z_full_in, kind(r_def))
       z_half_soc = REAL(z_half_in, kind(r_def))
       coszen_soc = REAL(coszen, kind(r_def))

       

       CALL socrates_interface(rad_lat_soc, rad_lon_soc, soc_lw_mode,  &
            tg_tmp_soc, q_soc, ozone_soc, co2_soc, h2_soc, t_surf_for_soc, p_full_soc, p_half_soc, &
            z_full_soc, z_half_soc, albedo_soc, coszen_soc, n_profile, n_layer,     &
            output_heating_rate_lw, output_soc_flux_lw_down, output_soc_flux_lw_up, &
            output_soc_spectral_olr = outputted_soc_spectral_olr, t_half_level_out = t_half_out)

       !tg_tmp_soc = tg_tmp_soc + output_heating_rate_lw*delta_t !Output heating rate in K/s, so is a temperature tendency
       surf_lw_down = REAL(output_soc_flux_lw_down(n_layer+1), kind(dp))
       !thd_lw_flux_up = REAL(output_soc_flux_lw_up, kind(dp))

       temp_tend(:) = temp_tend(:) + real(output_heating_rate_lw, kind(dp))
       
       ! SW calculation
       ! Retrieve output_heating_rate, and downward surface SW and LW fluxes

       soc_lw_mode = .FALSE.
       CALL socrates_interface(rad_lat_soc, rad_lon_soc, soc_lw_mode,  &
            tg_tmp_soc, q_soc, ozone_soc, co2_soc, h2_soc, t_surf_for_soc, p_full_soc, p_half_soc, &
            z_full_soc, z_half_soc, albedo_soc, coszen_soc, n_profile, n_layer,     &
            output_heating_rate_sw, output_soc_flux_sw_down, output_soc_flux_sw_up)

       !tg_tmp_soc = tg_tmp_soc + output_heating_rate_sw*delta_t !Output heating rate in K/s, so is a temperature tendency
! MDH added ST fix
       net_surf_sw_down = REAL(output_soc_flux_sw_down(n_layer+1)-output_soc_flux_sw_up(n_layer+1) , kind(dp))
       !thd_sw_flux_down = REAL(output_soc_flux_sw_down, kind(dp))


       temp_tend(:) = temp_tend(:) + real(output_heating_rate_sw, kind(dp))
       
       output_heating_rate_total = output_heating_rate_lw + output_heating_rate_sw
       net_flux = output_soc_flux_lw_up - output_soc_flux_lw_down + &
                  output_soc_flux_sw_up - output_soc_flux_sw_down

!hires
!            if(id_soc_spectral_olr > 0) then 
!                used = send_data ( id_soc_spectral_olr, outputted_soc_spectral_olr, Time_diag)
!            endif      

end subroutine run_socrates  

subroutine run_socrates_end

 !   use interpolator_mod, only: interpolator_end
    USE socrates_config_mod    

!    if(do_read_ozone) call interpolator_end(o3_interp)
!    if(do_read_co2)   call interpolator_end(co2_interp)
    

end subroutine run_socrates_end

!*****************************************************************************************
        subroutine interp_temp(z_full,z_half,temp_in, t_half)
          implicit none

          real(r_def),dimension(:,:),intent(in)  :: z_full,z_half,temp_in
          real(r_def),dimension(size(z_half,1), size(z_half,2)),intent(out) :: t_half

          integer i,k,kend
          real dzk,dzk1,dzk2
          
! note: z_full(kend) = z_half(kend), so there's something fishy
! also, for some reason, z_half(k=1)=0. so we need to deal with k=1 separately
          kend=size(z_full,2)
          do k=2,kend
                do i=1,size(temp_in,1)
                   dzk2 = 1./( z_full(i,k-1)   - z_full(i,k) )
                   dzk  = ( z_half(i,k  )   - z_full(i,k) )*dzk2
                   dzk1 = ( z_full(i,k-1)   - z_half(i,k) )*dzk2 
                   t_half(i,k) = temp_in(i,k)*dzk1 + temp_in(i,k-1)*dzk
                enddo
          enddo
! top of the atmosphere: need to extrapolate. z_half(1)=0, so need to use values
! on full grid
             do i=1,size(temp_in,1)
                !standard linear extrapolation
                !top: use full points, and distance is 1.5 from k=2
                t_half(i,1) = 0.5*(3*temp_in(i,1)-temp_in(i,2))
                !bottom: z=0 => distance is
                !-z_full(kend-1)/(z_full(kend)-z_full(kend-1))
                t_half(i,kend+1) = temp_in(i,kend-1) &
                     + (z_half(i,kend+1) - z_full(i,kend-1))&
                     * (temp_in(i,kend ) - temp_in(i,kend-1))&
                     / (z_full(i,kend  ) - z_full(i,kend-1))
             enddo


        end subroutine interp_temp
!*****************************************************************************************
  
end module socrates_interface_mod
