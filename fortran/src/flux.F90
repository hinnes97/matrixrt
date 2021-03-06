module flux_mod

  use params, only: dp, sb, invert_grid, moist_rad, surface
!  use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
  use condense, only : rain_out
  
#ifdef SOC
  use socrates_interface_mod, only: run_socrates
#elif defined SHORT_CHAR
  !use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
  use toon_mod, only : toon_driver
#elif defined PICKET
  use k_Rosseland_mod, only: k_func_Freedman_local, gam_func_Parmentier, AB_func_Parmentier
  use short_char_ross, only: short_char_ross_driver
  use radiation_mod, only : radiation_interface
#elif defined TWOSTR
  use band_grey_mod, only : run_twostr
#endif
  
  implicit none
contains
  subroutine get_fluxes(nf, ne, Tf, pf, Te, pe, &
       net_F, mu_s, Finc, Fint, olr, q, Ts, fup, fdn)
    !! Input variables
    integer, intent(in) :: nf, ne                         ! Number of layers, levels (lev = lay + 1)
    real(dp), dimension(:), intent(in) :: Tf, pf           ! Temperature [K], pressure [pa] at layers
    real(dp), dimension(:), intent(in) :: Te, pe           ! pressure [pa] at levels
    real(dp), intent(in) :: Finc, mu_s                        ! Incident flux [W m-2] and cosine zenith angle
    real(dp), intent(in) :: Fint                              ! Internal flux [W m-2]
    real(dp), intent(in) :: q(:)
    real(dp), intent(in) :: Ts                            ! Surface temperature
    
    !! Output variables
    real(dp), dimension(ne), intent(out) :: net_F, fup, fdn
    real(dp), intent(out) :: olr

    integer :: i
    
#ifdef SOC
    
    real(dp) :: rad_lat, rad_lon, t_surf_in, albedo_in, net_surf_sw_down, surf_lw_down
    real(dp), dimension(size(Tf)) :: q_in, temp_tend, h2_in, ch4_in
    
#elif defined PICKET
    real(dp), dimension(3) :: gam_V, Beta_V, A_Bond
    real(dp), dimension(2) :: beta
    real(dp) :: gam_1, gam_2, Tint, Tirr
    real(dp) :: met
    real(dp), dimension(2,nf) :: kIR_Ross
    real(dp), dimension(3,nf) :: kV_Ross        

#elif defined TWOSTR
    real(dp), allocatable :: delp(:)
#endif

#ifdef SOC
    
    rad_lat = 0._dp
    rad_lon = 0._dp
    q_in = 0.001_dp*9._dp
    q_in = 0.01_dp
    !ch4_in =  0.001_dp*8._dp
    !q_in = 0.01_dp
    h2_in =  1._dp - q_in! - ch4_in!1._dp - q_in
    albedo_in = 0._dp
    t_surf_in = Te(ne)
    
    !write(*,*) '-----------------------------------------------------------------------------'
    call run_socrates(rad_lat, rad_lon, Tf, q_in, h2_in, ch4_in, t_surf_in, pf, pe, pf, pe, albedo_in, &
         temp_tend, net_surf_sw_down, surf_lw_down, net_F)
    !do i=1,size(net_F)
    !   write(*,*) net_F(i)
    !end do
    
#elif defined SHORT_CHAR
    
    !call Kitzmann_TS_noscatt(nf, ne, Te, pe, &
    !     net_F, mu_s, Finc, Fint, olr, q, fup, fdn)

    if ((moist_rad .and. surface)) then
       call toon_driver(Te, pe, net_F, q, Ts)
    else if (moist_rad) then
       call toon_driver(Te, pe, net_F, q)
    else if (surface) then
       call toon_driver(Te, pe, net_F, Ts=Ts)
    else
       call toon_driver(Te, pe, net_F)
    endif
       

#elif defined PICKET
    met =0.0_dp
    Tint = (Fint/sb)**0.25_dp
    Tirr = (Finc/sb)**0.25_dp
    call gam_func_Parmentier(Tint, Tirr, 2, 0._dp, 0._dp, gam_V, Beta_V, Beta, gam_1, gam_2, A_Bond)
    do i=1,nf
       call k_func_Freedman_local(Tf(i), pf(i)*10, met, kIR_Ross(1,i))
       kIR_Ross(1,i) = kIR_Ross(1,i)*0.1_dp
       kV_Ross(:,i) = kIR_Ross(1,i)*gam_V

       kIR_Ross(2,i) = kIR_Ross(1,i) * gam_2
       kIR_Ross(1,i) = kIR_Ross(1,i) * gam_1
    end do

    !A_Bond = 0.0_dp
    call short_char_ross_driver(nf,ne,Te,pe,net_F,1.0_dp,Finc, Fint,olr,kV_Ross,kIR_Ross,Beta_V, &
    Beta, A_Bond)!, &
    !kV_Ross, kIR_Ross, Beta_V, Beta, A_Bond)

    !call radiation_interface(pe,pf,Tf,Tf(nf),net_F, &
    !     kV_Ross, kIR_Ross, Beta_V, Beta, A_bond)

#elif defined TWOSTR

    if (invert_grid) then
       allocate(delp(nf+1))
       do i=2,nf
          delp(i) = pf(i) - pf(i-1)
       enddo
       delp(1) = pf(1) - pe(1)
       delp(nf+1) = pe(nf+1) - pf(nf)
    else
       allocate(delp(nf))
       do i=1,nf
          delp(i) = pe(i+1) - pe(i)
       enddo
    endif
    
    call run_twostr(nf, Tf, Te, pf, pe, delp, q, net_F, olr)

    deallocate(delp)
#endif
    
    
  end subroutine get_fluxes
end module flux_mod
