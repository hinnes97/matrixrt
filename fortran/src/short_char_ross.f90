module short_char_ross
  use params, only: sb, pi, twopi, dp
  implicit none

  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  
contains
    subroutine short_char_ross_driver(nlay, nlev, Te, pe, &
          &  net_F, mu_s, Finc, Fint, olr, kV_R, kIR_R, Beta_V, Beta_IR, A_Bond)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev                         ! Number of layers, levels (lev = lay + 1)
    !real(dp), dimension(nlay), intent(in) :: Tl, pl           ! Temperature [K], pressure [pa] at layers
    real(dp), dimension(nlev), intent(in) :: Te, pe           ! pressure [pa] at levels

    real(dp), intent(in) :: Finc, mu_s                        ! Incident flux [W m-2] and cosine zenith angle
    real(dp), intent(in) :: Fint                              ! Internal flux [W m-2]
    real(dp), dimension(:,:), intent(in) :: kV_R, kIR_R
    real(dp), dimension(:), intent(in) :: Beta_V, Beta_IR, A_Bond

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F
    real(dp), intent(out) :: olr

    !! Work variables
    integer :: i, nV_b, nIR_b
    real(dp), dimension(nlev) :: be
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: sw_down_b, sw_up_b, lw_down_b, lw_up_b
    real(dp), dimension(nlev) :: lw_net, sw_net
    real(dp), dimension(nlev) :: tau_IR, tau_V

    real(dp) :: Finc_b
    
    nV_b = size(Beta_V)
    nIR_b = size(Beta_IR)
    ! Find temperature at layer edges through linear interpolation and extrapolation
    !do i = 2, nlay
    !  call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
      !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
    !end do
    ! Extrapolate to find Te at uppermost and lowest levels
    !Te(1) = Tl(1) + (pe(1) - pe(2))/(pl(1) - pe(2)) * (Tl(1) - Te(2))
    !Te(nlev) = Tl(nlay) + (pe(nlev) - pe(nlay))/(pl(nlay) - pe(nlay)) * (Tl(nlay) - Te(nlay))

    !! Shortwave fluxes
    sw_down(:) = 0.0_dp
    sw_up(:) = 0.0_dp

    do i=1,nV_b
       sw_down_b(:) = 0.0_dp
       call tau_struct(nlev,pe,kV_R(n,:), tau_V)
       if (mu_s > 0.0_dp) then
          Finc_b = Finc* Beta_V(n) * (1.0_dp - A_Bond(n))
          call sw_grey_down(nlev, Finc_b, tau_V, sw_down_b, mu_s)
       end if
       sw_down = sw_down + sw_down_b       
    end do
    
    !! Long wave two-stream fluxes
    ! Blackbody fluxes (note divide by pi for correct units)
    lw_down = 0.0_dp
    lw_up = 0.0_dp

    do n=1, nIR_b
       be(:) = sb * Te(:)**4/pi * Beta_IR(n)
       
       lw_down_b = 0.0_dp
       lw_up_b = 0.0_dp

       call tau_struct(nlev, pe, kIR_R(n,:), tau_IR)
       call lw_grey_updown_linear(nf, ne, be, tau_IR, lw_up_b, lw_down_b)

       lw_down = lw_down + lw_down_b
       lw_up   = lw_up   + lw_up_b
    end do
    
    net_F = lw_up + sw_up - lw_down - sw_down

    olr = lw_up(1)

  end subroutine short_char_ross_driver

  subroutine lw_grey_updown_linear(nlay, nlev, be, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, edel
    real(dp) :: del, e0i, e1i, e1i_del
    real(dp), dimension(nlay) :: Am, Bm
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau in each layer
    do k = 1, nlay
      dtau(k) = tau_IRe(k+1) - tau_IRe(k)
    end do

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Prepare loop
      do k = 1, nlay
        ! Olson & Kunasz (1987) linear interpolant parameters
        del = dtau(k)/uarr(m)
        edel(k) = exp(-del)
        e0i = 1.0_dp - edel(k)
        e1i = del - e0i

        e1i_del = e1i/del

        Am(k) = e0i - e1i_del ! Am(k) = Gp(k), just indexed differently
        Bm(k) = e1i_del ! Bm(k) = Bp(k), just indexed differently
      end do

      !! Begin two-stream loops
      !! Peform downward loop first
      ! Top boundary condition - 0 intensity downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*edel(k) + Am(k)*be(k) + Bm(k)*be(k+1) ! TS intensity
      end do

      !! Peform upward loop
      ! Lower boundary condition - planck function intensity upward from lowest level
      lw_up_g(nlev) = be(nlev)
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*edel(k) + Bm(k)*be(k) + Am(k)*be(k+1) ! TS intensity
      end do

      !! Sum up flux arrays with Gauss weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)
      
    end do

    ! Convert to flux by * 2pi
    lw_down(:) = twopi * lw_down(:)
    lw_up(:) = twopi * lw_up(:)


  end subroutine lw_grey_updown_linear

  subroutine sw_grey_down(nlev, solar, solar_tau, sw_down, mu)
    implicit none

    integer, intent(in) :: nlev
    real(dp), intent(in) :: solar, mu
    real(dp), dimension(nlev), intent(in) :: solar_tau
    real(dp), dimension(nlev), intent(out) :: sw_down

    sw_down(:) = solar * mu * exp(-solar_tau(:)/mu)
    
  end subroutine sw_grey_down

  ! Perform linear interpolation in log10 space
  pure subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp

end module short_char_ross
