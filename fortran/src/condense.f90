module condense

  use params, only : rdgas, q0, dp

  implicit none

  real(dp) :: rvgas = 461.52
  real(dp) :: L_sub = 2.84e6 ! Defined at triple point
  real(dp) :: L_vap = 2.49e6 ! As above
  real(dp) :: T_tp  = 2.7315e2 ! Triple point temp
  real(dp) :: p_tp  = 611.
  
contains

  subroutine rain_out(p,T,q, qsat)
    real(dp), intent(in) :: p, T
    real(dp), intent(inout) :: q
    real(dp), intent(out) :: qsat
    real(dp) :: eps, sat_p

    eps = rdgas/rvgas

    ! Find saturation vapour pressure
    call sat_vp(p, T, sat_p)
    
    qsat = (eps*sat_p/p)/(1. + (eps - 1.)* sat_p/p )
    
    !where (q .gt. qsat)
    !   q = qsat
    !elsewhere
    !   q = q0
    !endwhere
    if (q .gt. qsat) then
       q = qsat
    else
       q = q0
    endif
    
  end subroutine rain_out

  subroutine sat_vp(p, T, sat_p)
    real(dp), intent(in) :: p
    real(dp), intent(in) :: T
    real(dp), intent(out) :: sat_p

    real(dp) :: L

    if (T .gt. T_tp) then
       L = L_vap
    else
       L = L_sub
    endif

    sat_p = p_tp*exp(-L/rvgas * (1./T - 1./T_tp) )
    
  end subroutine sat_vp

  subroutine cold_trap(q)
    real(dp), intent(inout) :: q(:)

    integer :: i(1), j

    i = minloc(q)
    
    do j=1, i(1)
       q(j) = q(i(1))
    enddo
    
  end subroutine cold_trap
end module condense
