module tau_mod

  use params, only: grav, nf, dp
  
  implicit none
  
contains

  subroutine calc_tau(kappa, pe, q, mu_av, tau)
    ! Subroutine to calculate the optical depth of an atmospheric layer

    real(dp), intent(in) :: kappa ! Opacity
    real(dp), intent(in) :: pe(:) ! Edge pressure
    real(dp), intent(in) :: q(:)  ! Mass concentration of species
    real(dp), intent(in) :: mu_av ! Average value of cos(zenith)

    real(dp), intent(inout) :: tau(:)

    integer :: k

    ! Find additional optical depth of layer due to component with concentration q
    do k=1,nf
       tau(k+1) = tau(k+1) + tau(k) + (pe(k+1) - pe(k))*kappa/grav * q(k) / mu_av
    enddo
    
  end subroutine calc_tau
  
  
end module tau_mod
