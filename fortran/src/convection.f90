module convection

  use params, only : Rcp, n_iter, dp
  
  implicit none

contains

  subroutine dry_adjust(T, p, dry_mask)
    real(dp), dimension(:), intent(inout) :: T
    real(dp), dimension(:), intent(in) :: p
    logical, dimension(:), intent(inout) :: dry_mask

    ! Work variables
    integer :: n,k
    integer :: nlay 

    real(dp) :: pfact, Tbar
    real(dp), dimension(size(T)) :: dp

    nlay = size(T) - 1

    
    do k=2,nlay+1
       dp(k) = p(k) - p(k-1)
    end do
    dp(1) = dp(2)
    
    ! Number of passes
    do n=1,n_iter
       !Downwards pass
       do k=1,nlay
          pfact = (p(k)/p(k+1))**Rcp

          if ( T(k+1)/T(k) * pfact .gt. 1. ) then
             Tbar = (dp(k)*T(k) + dp(k+1)*T(k+1))/(dp(k) + dp(k+1))
             T(k+1) = (dp(k) + dp(k+1))*Tbar/(dp(k+1) + dp(k) * pfact)
             T(k) = T(k+1)*pfact
             dry_mask(k) = .True.
             dry_mask(k+1) = .True.

          else
             dry_mask(k) = .False.
             dry_mask(k+1) = .False.
          endif
       end do

       ! Upwards pass
       do k=nlay-1, 1, -1
                    pfact = (p(k)/p(k+1))**Rcp

          if ( T(k+1)/T(k) * pfact .gt. 1. ) then
             Tbar = (dp(k)*T(k) + dp(k+1)*T(k+1))/(dp(k) + dp(k+1))
             T(k+1) = (dp(k) + dp(k+1))*Tbar/(dp(k+1) + dp(k) * pfact)
             T(k) = T(k+1)*pfact
             dry_mask(k) = .True.
             dry_mask(k+1) = .True.
          else
             dry_mask(k) = .False.
             dry_mask(k+1) = .False.
          endif
       end do       
    end do
        
  end subroutine dry_adjust
  
  
end module convection
