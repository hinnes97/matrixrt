module flux_mod

  use params, only: dp, rad_scheme, sb
!  use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt

#ifdef SOC
  use socrates_interface_mod, only: run_socrates
#elif defined SHORT_CHAR
  use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
#elif defined PICKET
  use k_Rosseland_mod, only: k_func_Freedman_local, gam_func_Parmentier, AB_func_Parmentier
  use short_char_ross, only: short_char_ross_driver
  use radiation_mod, only : radiation_interface
#endif
  
  implicit none
contains
  subroutine get_fluxes(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
       net_F, mu_s, Finc, Fint, olr)

    !! Input variables
    integer, intent(in) :: nf, ne                         ! Number of layers, levels (lev = lay + 1)
    real(dp), dimension(:), intent(in) :: Tf, pf           ! Temperature [K], pressure [pa] at layers
    real(dp), dimension(:), intent(in) :: Te, pe           ! pressure [pa] at levels
    real(dp), dimension(:), intent(in) :: tau_V, tau_IR  ! V and IR band cumulative optical depth at levels
    real(dp), intent(in) :: Finc, mu_s                        ! Incident flux [W m-2] and cosine zenith angle
    real(dp), intent(in) :: Fint                              ! Internal flux [W m-2]
    
    !! Output variables
    real(dp), dimension(ne), intent(out) :: net_F
    real(dp), intent(out) :: olr

    integer :: i
    
#ifdef SOC
    
    real(dp) :: rad_lat, rad_lon, t_surf_in, albedo_in, net_surf_sw_down, surf_lw_down
    real(dp), dimension(size(Tf)) :: q_in, temp_tend, h2_in
    
#elif defined PICKET
    real(dp), dimension(3) :: gam_V, Beta_V, A_Bond
    real(dp), dimension(2) :: beta
    real(dp) :: gam_1, gam_2, Tint, Tirr
    real(dp) :: met
    real(dp), dimension(2,nf) :: kIR_Ross
    real(dp), dimension(3,nf) :: kV_Ross        
#endif

#ifdef SOC
    
    rad_lat = 0._dp
    rad_lon = 0._dp
    q_in = 0.5_dp
    !q_in = 0.01_dp
    h2_in = 0.0_dp!h2_in = 1._dp - q_in
    albedo_in = 0._dp
    t_surf_in = Te(ne)
    
    !write(*,*) '-----------------------------------------------------------------------------'
    call run_socrates(rad_lat, rad_lon, Tf, q_in, h2_in, t_surf_in, pf, pe, pf, pe, albedo_in, &
         temp_tend, net_surf_sw_down, surf_lw_down, net_F)
    !do i=1,size(net_F)
    !   write(*,*) net_F(i)
    !end do
    
#elif defined SHORT_CHAR
    
    call Kitzmann_TS_noscatt(nf, ne, Te, pe, tau_IR, tau_V, &
         net_F, mu_s, Finc, Fint, olr)

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
#endif
    
    
  end subroutine get_fluxes
end module flux_mod
