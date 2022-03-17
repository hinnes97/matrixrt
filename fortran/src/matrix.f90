module matrix
  use params, only: dp, mat_iters, q0, alpha, error_frac, Rcp, passes
  use flux_mod, only: get_fluxes
  use utils, only: linear_log_interp
  use condense, only: rain_out, cold_trap
!  use lapack, only: dgetrf
  
  implicit none

  real(dp) :: pert_T = 0.1_dp
contains
  
  subroutine do_matrix(nf, ne, Tf, pf, Te, pe, &
        mu_s, Finc, Fint, olr, q, Ts)

    integer, intent(in) :: nf, ne
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(in), dimension(:) :: pf, pe
    real(dp), intent(inout) :: Te(:), Tf(:)
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: q(:)
    real(dp), intent(inout) :: Ts
    
    ! Work variables
    real(dp), dimension(ne) :: residual, del_T, delp, T_old
    real(dp), dimension(ne,ne) :: mat
    real(dp), dimension(nf) :: q_sat
    logical :: dry_mask(ne), do_conv
    real(dp) grad, pfact
    integer :: n,i, m,k

    q(1:nf) = q0
    ! Rain out excess and then
!!$       open(10,file='hi.out')
!       do i=1,ne
!          call rain_out(pe(i), Te(i), q(i), q_sat(i))
          !write(*,*) q_sat(i)
!!$          write(10,*) q_sat(i)
 !      enddo
!!$       close(10)
 !      call cold_trap(q)

    dry_mask = .false.
    do_conv = .false.
    do m=1,2
    do n=1,mat_iters
       write(*,*) n
!!$       
!!$       if (mod(n,5) .eq. 0) then
!!$          do i=1,nf
!!$             call rain_out(pf(i), Tf(i), q(i), q_sat(i))
!!$          enddo
!!$          call cold_trap(q)
!!$       endif
!!$       write(*,*) q
       call invert_this_matrix(del_T, nf, ne, Tf, pf, pe,  mu_s, Finc, Fint, olr, &
            residual, Te, q, alpha, dry_mask, do_conv, Ts)
!       do i=1,ne
!          write(*,*) del_T(i)
!       enddo
       
!       call calc_matrix(mat, nf, ne, Tf, pf, pe, tau_IR, tau_V, mu_s, Finc, Fint, olr, residual, Te, q)
!       call solve_matrix(mat, del_T, residual, ne)

       !call invert_this_matrix(del_T, nf, ne, Tf, pf, pe, tau_IR, tau_V, mu_s, Finc, Fint, olr, residual, Te, q, alpha)
       
       Te = Te + del_T
       do i=1,ne
          write(*,*) residual(i), del_T(i), Te(i), dry_mask(i)
          if (Te(i) .lt.  100.0_dp) Te(i) = 100._dp
          if (Te(i) .gt. 2000.0_dp) Te(i) = 2000._dp
       end do
    end do

    ! Now do adiabatic adjustment

!!$    do k=2,nf
!!$       delp(k) = pf(k) - pf(k-1)
!!$    enddo
!!$    delp(1) = delp(2)
!!$    delp(ne) = delp(nf)
!!$
!!$    T_old = Te
!!$
!!$    do i=1,1000
!!$    do k=2,ne
!!$          grad = log(Te(k)/Te(k-1))/log(pe(k)/pe(k-1))
!!$          !pfact = exp(Rcp * log(p(k+1)/p(k)))
!!$          !if ( T(k)/T(k+1) * pfact .gt. 1. ) then
!!$
!!$          if ( grad .gt. Rcp) then
!!$             !write(*,*) 'TRUE', k
!!$             pfact = (pe(k)/pe(k-1))**(Rcp + 0.015)
!!$
!!$             Te(k-1) = ( Te(k-1)*delp(k-1) + Te(k)*delp(k))/( delp(k-1) + delp(k) * pfact)
!!$
!!$             Te(k) = Te(k-1) * pfact
!!$             
!!$             dry_mask(k) = .True.
!!$
!!$             !Te(k-1), Te(k)
!!$             !dry_mask(k-1) = .True.
!!$          else
!!$             dry_mask(k) = .False.
!!$             !Te(k) = 1.
!!$             !dry_mask(k-1) = .False.
!!$          endif
!!$
!!$          
!!$       end do
!!$
!!$       ! Upwards pass
!!$       do k=ne, 2, -1
!!$
!!$          grad = log(Te(k)/Te(k-1))/log(pe(k)/pe(k-1))
!!$          !pfact = expe(Rcp * log(pe(k+1)/pe(k)))
!!$          !if ( Te(k)/Te(k+1) * pfact .gt. 1. ) then
!!$         if ( grad .gt. Rcp) then
!!$             !write(*,*) 'TRUE', k
!!$             pfact = (pe(k)/pe(k-1))**(Rcp + 0.015)
!!$
!!$             Te(k-1) = ( Te(k-1)*delp(k-1) + Te(k)*delp(k))/( delp(k-1) + delp(k) * pfact)
!!$
!!$             Te(k) = Te(k-1) * pfact
!!$
!!$             dry_mask(k) = .True.
!!$            !dry_mask(k-1) = .True.
!!$
!!$          else
!!$             dry_mask(k) = .False.
!!$             !dry_mask(k-) = .False.
!!$          endif
!!$
!!$       end do       
!!$       enddo

    !call matrix_adiabat(ne, pe, delp, Te, dry_mask)
    do i=1,nf
       call linear_log_interp(pf(i), pe(i), pe(i+1), Te(i), Te(i+1), Tf(i))
    end do

!    do_conv = .true.
    
 enddo
   end subroutine do_matrix

   subroutine invert_this_matrix(delT, nf, ne, Tf, pf, pe,  mu_s, Finc, Fint, olr, residual, &
        Te, q, alpha, dry_mask, do_conv, Ts)
    integer, intent(in) :: nf, ne
    real(dp), intent(in) ::  pe(:), pf(:), alpha
    real(dp), intent(out) :: delT(:)
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: residual(:), Te(:), Tf(:)
    real(dp), intent(in) :: q(:)
    real(dp), intent(inout) :: Ts
    logical, intent(in) :: do_conv
    logical, intent(out) :: dry_mask(:)

    real(dp) :: H(ne,ne), S(ne,ne), E(ne,ne), SH_T(ne,ne), HSH_T(ne,ne), mat(ne,ne), F(ne,ne)
    real(dp) :: A(ne,ne), B(ne,ne), X(ne), Y(ne), Z(ne)
    
    integer :: i,j,k, ierr
    
    ! Calculate the flux matrix
    call calc_matrix(H, nf, ne, Tf, pf, pe,  mu_s, Finc, Fint, olr, residual, Te, q, dry_mask, do_conv, Ts)

    ! Calculate prior covariance matrix
    call s_matrix(ne, pe, S)
    ! Calculate error matrix
    call e_matrix(ne, Fint, E)

    ! Calculate SH^T
    SH_T = 0.0_dp
    do i=1,ne
       do j=1,ne
          do k=1,ne
             SH_T(i,j) = SH_T(i,j) + S(i,k)*H(j,k)
          enddo
       enddo
    enddo
    
    ! Calculate HSH^T
    HSH_T = 0.0_dp
    do i=1,ne
       do j=1,ne
          do k=1,ne
             HSH_T(i,j) = HSH_T(i,j) + H(i,k)*SH_T(k,j)
          enddo
       enddo
    enddo

    ! Calculate matrix to invert
    do i=1,ne
       do j=1,ne
          mat(i,j) = alpha*HSH_T(i,j) + E(i,j)
       enddo
    enddo

    ! Now invert this matrix using LAPACK

    ! LU decomposition
    call dgetrf(ne, ne, mat, ne, X, ierr)
    if (ierr .ne. 0) write(*,*) 'WARNING, LU DECOMPOSITION ERROR: ', ierr

    ! Matrix inversion
    call dgetri(ne, mat, ne, X, Y, ne, ierr)
    if (ierr .ne. 0) write(*,*) 'WARNING, INVERSION ERROR: ', ierr
    
    ! Calculate final matrix
    F = 0.0_dp
    do i=1,ne
       do j=1,ne
          do k=1,ne
             F(i,j) = F(i,j) + alpha * SH_T(i,k) * mat(k,j)
          enddo
       enddo
    enddo

    ! Calculate delT
    delT = 0._dp
    do i=1,ne
       do k=1,ne
          delT(i) =  delT(i) - F(i,k) * residual(k)
       enddo
    enddo
    
    
  end subroutine invert_this_matrix
  
  subroutine calc_matrix(mat, nf, ne, Tf, pf, pe,  mu_s, Finc, Fint, olr, residual,Te,q, dry_mask, do_conv, Ts)
    integer, intent(in) :: nf, ne
    real(dp), intent(in) ::  pe(:), pf(:)
    real(dp), intent(out) :: mat(:,:)
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: residual(:), Te(:), Tf(:)
    real(dp), intent(in) :: q(:)
    real(dp), intent(inout) :: Ts
    logical, intent(out) :: dry_mask(:)
    logical, intent(in) :: do_conv
    
    
    real(dp), dimension(ne) :: Tpert, flux, flux_pert, qpert, qsattemp, fup,fdn
    integer :: i,j

    do i=1,nf
       call linear_log_interp(pf(i), pe(i), pe(i+1), Te(i), Te(i+1), Tf(i))
    end do

    call get_fluxes(nf, ne, Tf, pf, Te, pe,  &
         flux, mu_s, Finc, Fint, olr,q, Ts, fup,fdn)

    residual = flux - Fint
    
    do i=1,ne
       Tpert = Te
       !qpert = q
       
       Tpert(i) = Tpert(i) + pert_T
!       call rain_out(pe(i), Tpert(i), qpert(i), qsattemp(i))
       
       do j=1,nf
          call linear_log_interp(pf(j), pe(j), pe(j+1), Tpert(j), Tpert(j+1), Tf(j))
       end do
       
!       call get_fluxes(nf, ne, Tf, pf, Tpert, pe, tau_IR, tau_V, &
!            flux_pert, mu_s, Finc, Fint, olr, qpert,fup,fdn)
       call get_fluxes(nf, ne, Tf, pf, Tpert, pe, &
            flux_pert, mu_s, Finc, Fint, olr, q, Ts, fup,fdn)
       
       mat(:,i) = (flux_pert - flux)/pert_T
    end do

    ! Convective part
    if (do_conv) call adiabat_flux(ne, dry_mask, residual, mat, Fint, pe, Te)
        
    open(8,file='mat.out')
    do i=1,ne
       write(8,'(ES13.4)') (mat(i,j),j=1,ne)
    end do
    close(8)
  end subroutine calc_matrix

  subroutine solve_matrix(mat, del_T, residual, ne)
    integer, intent(in) :: ne
    real(dp), intent(in) :: mat(ne,ne), residual(ne)
    real(dp), intent(out) :: del_T(ne)

    ! Work variables
    real(dp), dimension(ne,ne) :: factors,newmat
    integer :: ipiv(ne), info, lwork,nrhs
    real(dp), dimension(:), allocatable :: work
    real(dp),dimension(ne) ::  R, C
    character(len=1) :: equed
    real(dp) :: rcond
    real(dp), dimension(1) :: ferr, berr
    integer, dimension(ne) :: iwork
    real(dp), dimension(ne,1) :: input
    real(dp), dimension(ne,1) :: output
    nrhs = 1
    
    lwork = ne*4
    allocate(work(lwork))

    newmat = mat
    input(:,nrhs) = -residual
    Call dgesvx('E', 'N', ne, 1, newmat, ne, factors, ne, &
         ipiv, equed, r, c, input, &
         ne, output, ne, rcond, ferr, berr, work, iwork, &
         info)


    del_T = output(:,nrhs)

  end subroutine solve_matrix

  subroutine s_matrix(ne, pe, mat_s)
    integer, intent(in) :: ne
    real(dp), intent(in) :: pe(ne)
    real(dp), intent(out) :: mat_s(ne,ne)

    integer :: i,j

    do i=1,ne
       do j=1,ne
          mat_s(i,j) = exp(-(log(pe(i)/pe(j)))**2._dp/2._dp/0.25_dp)
       enddo
    enddo
    
  end subroutine s_matrix

  subroutine e_matrix(ne, Fint, mat_e)
    integer, intent(in) :: ne
    real(dp), intent(in) :: Fint
    real(dp), intent(out) :: mat_e(ne,ne)

    integer :: i

    mat_e = 0.0_dp

    do i=1,ne
       mat_e(i,i) = Fint*error_frac
    enddo
    
  end subroutine e_matrix

  subroutine matrix_adiabat(nlay, p, delp, T, dry_mask)

    integer, intent(in) :: nlay
    real(dp), intent(in) :: p(nlay), delp(nlay)
    real(dp), intent(inout) :: T(nlay)
    logical, intent(out) :: dry_mask(nlay)
    
    integer :: k,n
    real(dp) :: pfact, grad

    do n=1,1
       !Downwards pass
       do k=2,nlay
          grad = log(T(k)/T(k-1))/log(p(k)/p(k-1))
          !pfact = exp(Rcp * log(p(k+1)/p(k)))
          !if ( T(k)/T(k+1) * pfact .gt. 1. ) then
          if ( grad .gt. Rcp) then
             !write(*,*) 'TRUE', k
             pfact = (p(k)/p(k-1))**(Rcp + 0.015)

             T(k-1) = ( T(k-1)*delp(k-1) + T(k)*delp(k))/( delp(k-1) + delp(k) * pfact)

             T(k) = T(k-1) * pfact
             
             dry_mask(k) = .True.
             write(*,*) 'HERE', (log(T(k)) - log(T(k-1))) / (log(p(k)) - log(p(k-1))), log(T(k)/T(k-1))/log(p(k)/p(k-1))
             !dry_mask(k-1) = .True.

          else
             dry_mask(k) = .False.
             !dry_mask(k-1) = .False.
          endif
          write(*,*) T(k)
       end do

       

       ! Upwards pass
       do k=nlay, 2, -1

          grad = log(T(k)/T(k-1))/log(p(k)/p(k-1))
          !pfact = exp(Rcp * log(p(k+1)/p(k)))
          !if ( T(k)/T(k+1) * pfact .gt. 1. ) then
         if ( grad .gt. Rcp) then
             !write(*,*) 'TRUE', k
             pfact = (p(k)/p(k-1))**(Rcp + 0.015)

             T(k-1) = ( T(k-1)*delp(k-1) + T(k)*delp(k))/( delp(k-1) + delp(k) * pfact)

             T(k) = T(k-1) * pfact

             dry_mask(k) = .True.
            !dry_mask(k-1) = .True.

          else
             dry_mask(k) = .False.
             !dry_mask(k-) = .False.
          endif

       end do       
    end do
         
  end subroutine matrix_adiabat

  subroutine adiabat_flux(nlay, dry_mask, flux, H, Fint, p, T)
    integer, intent(in) :: nlay
    logical, intent(out) :: dry_mask(nlay)
    real(dp), intent(in) :: p(nlay), T(nlay), Fint
    
    real(dp), intent(inout) :: flux(nlay), H(nlay, nlay)

    integer :: k
    real(dp) :: grad

    do k=nlay,2,-1
       grad = (log(T(k)/T(k-1)) / log(p(k)/p(k-1)) )
       if (grad .gt. Rcp) then
          dry_mask(k) = .true.
          !grad = (log(T(k)/T(k-1)) / log(p(k)/p(k-1)) )
       !   write(*,*) grad, grad/Rcp - 1._dp
          flux(k) = flux(k) +  Fint * 1e3 * (grad/Rcp - 1._dp)**2!exp(200 * (grad/Rcp  - 1._dp))

          H(k,k) = H(k,k) + 2e3 * Fint * (grad/Rcp - 1)/ (T(k) * Rcp * log(p(k)/p(k-1)))!exp(200* (grad/Rcp - 1._dp)) / (T(k) * Rcp * log(p(k)/p(k-1)))
          H(k,k-1) = H(k,k) - 2e3 * Fint * (grad/Rcp - 1)/ (T(k-1) * Rcp * log(p(k)/p(k-1)))
          !H(k,k-1) = H(k,k-1) - 0.2 * Fint * exp(200* (grad/Rcp - 1._dp)) / (T(k-1) * Rcp * log(p(k)/p(k-1)))
       endif
    enddo
       
  end subroutine adiabat_flux
end module matrix
