module timestep
  use params, only: dp, Nt, nf, ne, const, Finc, Fint
  use flux_mod, only: get_fluxes
  use utils, only: linear_log_interp
  implicit none

contains
  subroutine step(Tf, pf, pe, tau_IR, tau_v, net_F, dT, olr)
    real(dp), dimension(:), intent(inout) :: Tf
    real(dp), dimension(:), intent(out) :: dT
    real(dp), dimension(:), intent(in) :: pf
    real(dp), dimension(:), intent(in) :: pe, tau_IR, tau_v
    real(dp), dimension(:), intent(out) :: net_F
    real(dp), intent(out) :: olr

    integer :: i,j
    real(dp), dimension(size(pe)) :: Te
    real(dp), dimension(size(Tf)) :: Tf_half
    
    do j =1,Nt
       write(*,*) ' -------------------------------------------------'
       write(*,*) 'Timestep ', j, ': Max(abs(res))', maxval(abs(net_F - Fint))
     if (mod(j, 100000) .eq. 0) then
        write(*,*) Tf
     end if

     do i = 2, nf
       call linear_log_interp(pe(i), pf(i-1), pf(i), Tf(i-1), Tf(i), Te(i))
    end do
    
    ! Extrapolate to find Te at uppermost and lowest levels
    Te(1) = Tf(1) + (pe(1) - pe(2))/(pf(1) - pe(2)) * (Tf(1) - Te(2))
    Te(ne) = Tf(nf) + (pe(ne) - pe(nf))/(pf(nf) - pe(nf)) * (Tf(nf) - Te(nf))

     call get_fluxes(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
          net_F, 1._dp, Finc, Fint, olr)

     do i=1,nf
        dT(i) = const*(net_F(i+1) - net_F(i))
        if (dT(i)>5.0_dp) then
           dT(i) = 5.0_dp
        endif
        if (dT(i)<-5.0_dp) then
           dT(i) = -5.0_dp
        endif
     end do

     Tf_half = Tf + 0.5_dp*dT
     
     do i = 2, nf
       call linear_log_interp(pe(i), pf(i-1), pf(i), Tf_half(i-1), Tf_half(i), Te(i))
    end do
    
    ! Extrapolate to find Te at uppermost and lowest levels
    Te(1) = Tf_half(1) + (pe(1) - pe(2))/(pf(1) - pe(2)) * ( Tf_half(1) - Te(2))
    Te(ne) = Tf_half(nf) + (pe(ne) - pe(nf))/(pf(nf) - pe(nf)) * (Tf_half(nf) - Te(nf))

     call get_fluxes(nf, ne, Tf_half, pf, Te, pe, tau_IR, tau_V, &
          net_F, 1._dp, Finc, Fint, olr)

     ! Internal flux at lowest level     
     net_F(ne) = Fint
     !do i=1,ne        
     !   write(*,*) net_F(i)
     !end do
     
     do i=1,nf
        dT(i) = const*(net_F(i+1) - net_F(i))
        if (dT(i)>5.0_dp) then
           dT(i) = 5.0_dp
        endif
        if (dT(i)<-5.0_dp) then
           dT(i) = -5.0_dp
        endif
        Tf(i) = Tf(i) + dT(i)
        if (Tf(i) < 100._dp) then
           Tf(i) = 100._dp
        end if        
     end do
     if (maxval(abs(net_F-Fint)) .lt. 1_dp) then
        write(*,*) 'Exiting'        
        exit
     end if
     
  end do
  end subroutine step
  
end module timestep
