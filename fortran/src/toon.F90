module toon_mod
  use params, only: ne, nf, sb, Finc, surface, A_s, kappa_sw, kappa_lw, moist_rad, kappa_q, dp
  use tau_mod, only: calc_tau
  
  implicit none
contains

  subroutine toon_driver(Te, pe, net_F, & 
                         q, Ts)
    real(dp), intent(in) :: Te(:), pe(:) !Edge temperature+pressures
    real(dp), intent(in), optional :: q(:) ! Optional moisture variable
    real(dp), intent(in), optional :: Ts ! Surface temperature

    real(dp), intent(out) :: net_F(:) ! Net flux
    
    real(dp) :: tau_lw(ne), tau_sw(ne) !Optical depths
    real(dp) :: mu_av = 0.5_dp
    real(dp) :: q_temp(nf)

    real(dp) :: F_u(ne), F_d(ne) ! Up/down IR fluxes
    real(dp) :: S_u(ne), S_d(ne) ! Up/down SW fluxes

    real(dp) :: F_up_bound ! Boundary condition
    
    tau_lw = 0._dp
    tau_sw = 0._dp
    q_temp = 1.
    
    if(moist_rad) then
       call calc_tau(kappa_q, pe, q, mu_av, tau_lw)
       call calc_tau(kappa_lw, pe,(1.-q), mu_av, tau_lw)
       call calc_tau(kappa_sw, pe, q_temp, mu_av, tau_sw)
    else
       call calc_tau(kappa_lw, pe, q_temp, mu_av, tau_lw)
       call calc_tau(kappa_sw, pe, q_temp, mu_av, tau_sw)
    endif

    ! Perform LW updown
    call lw_down(Te, tau_lw, mu_av, F_d)

    if(surface) then
       F_up_bound = sb*Ts**4
    else
       F_up_bound = sb*Te(ne)**4
    endif

    call lw_up(Te, tau_lw, mu_av, F_up_bound, F_u)

    call sw_down(tau_sw, S_d)
    
    S_u = 0.0_dp
    if (surface) call sw_up(tau_sw, S_d(ne), S_u)

    net_F = F_u + S_u - F_d - S_d
       
  end subroutine toon_driver

  subroutine lw_down(Te, tau_lw, mu_av, F_d)
    real(dp), intent(in)  :: Te(ne) ! Edge temperature
    real(dp), intent(in) :: tau_lw(ne) ! Optical depth
    real(dp), intent(in) :: mu_av   ! Average cos(zen)

    real(dp), intent(out) :: F_d(ne) ! Downwards fluxes
    
    real(dp) :: B_dash, B(ne), B_1m, B_2m, trans, alpha

    integer :: k

    B = sb*Te**4
    alpha = 2*mu_av

    F_d = 0.0_dp

    ! Perform algorithm (Heng 2014, eqns. 41, 47)
    do k=1,nf
       B_dash = sb*(B(k+1) - B(k))/mu_av/(tau_lw(k+1) - tau_lw(k))

       B_1m = B(k) - B_dash/2.
       B_2m = B(k) - B_dash*mu_av*(tau_lw(k+1) - tau_lw(k))
       trans = exp(-alpha*(tau_lw(k+1) - tau_lw(k)))

       F_d(k+1) = F_d(k)*trans + (B_2m - B_1m*trans)
    enddo
  end subroutine lw_down

  subroutine lw_up(Te, tau_lw, mu_av, BC, F_u)
    real(dp), intent(in)  :: Te(ne) ! Edge temperature
    real(dp), intent(in) :: tau_lw(ne) ! Optical depth
    real(dp), intent(in) :: mu_av   ! Average cos(zen)
    real(dp), intent(in) :: BC  ! Boundary condition at bottom of atm.
    real(dp), intent(out) :: F_u(ne) ! Upwards fluxes

    real(dp) :: B_dash, B(ne), B_1p, B_2p, trans, alpha

    integer :: k
    
    B = sb*Te**4
    alpha = 2*mu_av

    F_u(ne) = BC

    do k=nf,1,-1
       B_dash = sb*(B(k+1) - B(k))/mu_av/(tau_lw(k+1) - tau_lw(k))
       
       B_1p = B(k+1) + B_dash*mu_av*(tau_lw(k) - tau_lw(k+1)) + B_dash/2.
       B_2p = B(k+1) + B_dash/2.

       F_u(k) = trans*F_u(k+1) + (B_1p - B_2p*trans)
    enddo
  end subroutine lw_up
  
  subroutine sw_down(tau_sw, S_d)
    real(dp), intent(in) :: tau_sw(:) ! SW optical depth

    real(dp), intent(out) :: S_d(:) ! Downwelling SW radiation
      
    integer :: k

    S_d(1) = Finc
      
    do k=1,nf
       S_d(k+1) = S_d(k)*exp(-(tau_sw(k+1) - tau_sw(k)))
    enddo
      
  end subroutine sw_down

    subroutine sw_up(tau_sw, S_d_bound, S_u)
      ! Right now, only for use when surface is present
      real(dp), intent(in) :: tau_sw(:)
      real(dp), intent(in) :: S_d_bound ! Shortwave hitting surface

      real(dp), intent(out) :: S_u(:)
      
      integer :: k
      
      S_u(ne) = S_d_bound * A_s

      do k=nf,1,-1
         S_u(k) = S_u(k+1)*exp(-(tau_sw(k+1) - tau_sw(k)))
      enddo
      
    end subroutine sw_up
end module toon_mod
