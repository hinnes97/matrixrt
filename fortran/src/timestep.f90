module timestep
  use params, only: dp, Nt, nf, ne, const, Finc, Fint, q0
  use flux_mod, only: get_fluxes
  use utils, only: linear_log_interp
  use convection, only: dry_adjust
  use condense, only : rain_out, cold_trap
  use io, only: dump_data
  implicit none

contains
  subroutine step(Tf, pf, pe, tau_IR, tau_v, net_F, dT, olr, ncid, q, fup, fdn)
    integer, intent(in) :: ncid
    real(dp), dimension(:), intent(inout) :: Tf
    real(dp), dimension(:), intent(out) :: dT
    real(dp), dimension(:), intent(in) :: pf
    real(dp), dimension(:), intent(in) :: pe, tau_IR, tau_v
    real(dp), dimension(:), intent(out) :: net_F
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: q(:), fup(:), fdn(:)

    integer :: i,j
    real(dp), dimension(size(pe)) :: Te, q_half, q_sat
    real(dp), dimension(size(Tf)) :: Tf_half, temp, dT_old
    logical, dimension(size(Tf)) :: dry_mask
    real(dp), dimension(size(Tf)) :: factor
    integer, dimension(size(Tf)) :: oscillate
    
    
    factor = const
    temp = 0
    oscillate = 0
    dT = 1._dp
    dry_mask = .false.

    q(1:nf) = q0
    do j =1,Nt
       write(*,*) ' -------------------------------------------------'
       write(*,*) 'Timestep ', j, ': Max(abs(res))', maxval(abs(net_F - Fint)), ', Max(abs(dT)): ', maxval(abs(Tf - temp))
       write(*,*) 'OLR - Fint - Finc', net_F(1) - Fint
!     if (mod(j, 100) .eq. 0) then
!        write(*,*) Tf
!        call dump_data(ncid, nf, ne, Tf, pf, pe, olr, tau_IR(ne), tau_V(ne), Finc, Fint, Te)
!     end if
     
     temp = Tf
     
     do i = 2, nf
       call linear_log_interp(pe(i), pf(i-1), pf(i), Tf(i-1), Tf(i), Te(i))
    end do
    
    ! Extrapolate to find Te at uppermost and lowest levels
    !Te(1) = Tf(1) + (pe(1) - pe(2))/(pf(1) - pe(2)) * (Tf(1) - Te(2))
    !Te(ne) = Tf(nf) + (pe(ne) - pe(nf))/(pf(nf) - pe(nf)) * (Tf(nf) - Te(nf))
    call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

!    do i=1,ne
!       call rain_out(pe(i), Te(i), q(i), q_sat(i))
!    enddo
!    call cold_trap(q)
    
     call get_fluxes(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
          net_F, 1._dp, Finc, Fint, olr, q, fup, fdn)
!!$     do i=1,nf+1
!!$        write(*,*) 'NET FLUX'
!!$        write(*,*) net_F(i)
!!$     enddo
     
     do i=1,nf
        dT_old(i) = dT(i)
        dT(i) = const*(net_F(i+1) - net_F(i))/(pe(i+1) - pe(i))*pe(nf+1)
        !dT(i) = factor(i)*(net_F(i+1) - net_F(i))/(abs(net_F(i+1) - net_F(i))**0.9_dp)
        !*(pe(i+1) - pe(i))/pe(nf+1)
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
    !Te(1) = Tf_half(1) + (pe(1) - pe(2))/(pf(1) - pe(2)) * ( Tf_half(1) - Te(2))
    !Te(ne) = Tf_half(nf) + (pe(ne) - pe(nf))/(pf(nf) - pe(nf)) * (Tf_half(nf) - Te(nf))
    call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

!    do i=1,ne
!       call rain_out(pe(i), Te(i), q(i), q_sat(i))
!    enddo
!    call cold_trap(q)
    
    call get_fluxes(nf, ne, Tf_half, pf, Te, pe, tau_IR, tau_V, &
          net_F, 1._dp, Finc, Fint, olr, q, fup, fdn)

     ! Internal flux at lowest level     
     net_F(ne) = Fint
     !do i=1,ne        
     !   write(*,*) net_F(i)
     !end do
     
     do i=1,nf
        dT(i) = const*(net_F(i+1) - net_F(i))*pe(nf+1)/(pe(i+1) - pe(i))
        !dT(i) = factor(i)*(net_F(i+1) - net_F(i))/(abs(net_F(i+1) - net_F(i))**0.9_dp)
        if (dT(i)>5.0_dp) then
           dT(i) = 5.0_dp
        endif
        if (dT(i)<-5.0_dp) then
           dT(i) = -5.0_dp
        endif
        Tf(i) = Tf(i) + dT(i)

        !write(*,*) dT(i)
        !writE(*,*) Tf(i)
        if (Tf(i) < 100._dp) then
           Tf(i) = 100._dp
        else if (Tf(i) .gt. 1000. ) then
           Tf(i) = 1000._dp
        end if
        ! Check for oscillations
!!$
!!$        if (dry_mask(i) .eqv. .false.) then
!!$        
!!$          if (dT(i)*dT_old(i) .lt. 0._dp) then
!!$             oscillate(i) = oscillate(i) + 1
!!$          endif
!!$        
!!$          if (oscillate(i) .gt. 5) then
!!$             factor(i) = 0.66_dp*factor(i)
!!$             oscillate(i) = oscillate(i) - 1
!!$          else if ((oscillate(i)) .lt. 6 .and. (factor(i) .lt. 0.01)) then
!!$             factor(i) = 1.1_dp* factor(i)
!!$          endif
!!$                
!!$!          write(*,*) factor(i)
!!$       endif
     
     end do

     
     !do i = 2, nf
      ! call linear_log_interp(pe(i), pf(i-1), pf(i), Tf(i-1), Tf(i), Te(i))
    !end do
    
    ! Extrapolate to find Te at uppermost and lowest levels
    !call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
    !call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )
    !Te(1) = Tf(1) + (pe(1) - pe(2))/(pf(1) - pe(2)) * ( Tf(1) - Te(2))
    !Te(ne) = Tf(nf) + (pe(ne) - pe(nf))/(pf(nf) - pe(nf)) * (Tf(nf) - Te(nf))
    !call dry_adjust(Tf, pf, dry_mask)
    
    !temp = Tf
    !do i = 1, nf
    !   call linear_log_interp(pf(i), pe(i), pe(i+1), Te(i), Te(i+1), Tf(i))
    !end do

    !do i=1,nf
    !   write(*,*) Tf(i) - temp(i)
    !enddo
    
     !if (maxval(abs(net_F-Fint)) .lt. 1_dp) then
     !   write(*,*) 'Exiting'        
     !   exit
     !end if
     
  end do

  do i=1,nf
     write(*,*) Tf(i), Te(i)
  enddo
  
     
! do i =1,nf-1
!    writE(*,*) dry_mask(i)
!end do    
 
  end subroutine step
  
end module timestep
