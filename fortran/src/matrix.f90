module matrix
  use params, only: dp, N_max
  use flux_mod, only: get_fluxes
  use utils, only: linear_log_interp
!  use lapack, only: dgetrf
  
  implicit none

  real(dp) :: pert_T = 0.1_dp
contains
  
  subroutine do_matrix(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
        mu_s, Finc, Fint, olr)

    integer, intent(in) :: nf, ne
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(in), dimension(:) :: pf, pe, tau_IR, tau_V
    real(dp), intent(inout) :: Te(:), Tf(:)
    real(dp), intent(out) :: olr

    ! Work variables
    real(dp), dimension(ne) :: residual, del_T
    real(dp), dimension(ne,ne) :: mat
    integer :: n,i
    
    do n=1,N_max
       write(*,*) n
       call calc_matrix(mat, nf, ne, Tf, pf, pe, tau_IR, tau_V, mu_s, Finc, Fint, olr, residual, Te)
       call solve_matrix(mat, del_T, residual, ne)

       Te = Te + del_T
       do i=1,ne
          write(*,*) residual(i)
          if (Te(i) .lt.  100.0_dp) Te(i) = 100.0_dp
          if (Te(i) .gt. 5000.0_dp) Te(i) = 5000._dp
       end do
       write(*,*) '------------------------------------------------------------------------------------'
       do i=1,ne
          write(*,*) Te(i) 
       end do
    end do

  end subroutine do_matrix

  subroutine calc_matrix(mat, nf, ne, Tf, pf, pe, tau_IR, tau_V, mu_s, Finc, Fint, olr, residual,Te)
    integer, intent(in) :: nf, ne
    real(dp), intent(in) ::  pe(:), pf(:), tau_IR(:), tau_V(:)
    real(dp), intent(out) :: mat(:,:)
    real(dp), intent(in) :: mu_s, Finc, Fint
    real(dp), intent(out) :: olr
    real(dp), intent(inout) :: residual(:), Te(:), Tf(:)

    real(dp), dimension(ne) :: Tpert, flux, flux_pert
    integer :: i,j

    do i=1,nf
       call linear_log_interp(pf(i), pe(i), pe(i+1), Te(i), Te(i+1), Tf(i))
    end do
    
    call get_fluxes(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, &
         flux, mu_s, Finc, Fint, olr)

    residual = flux - Fint
    
    do i=1,ne
       Tpert = Te
       Tpert(i) = Tpert(i) + pert_T
       
       do j=1,nf
          call linear_log_interp(pf(j), pe(j), pe(j+1), Tpert(j), Tpert(j+1), Tf(j))
       end do
       
       call get_fluxes(nf, ne, Tf, pf, Tpert, pe, tau_IR, tau_V, &
            flux_pert, mu_s, Finc, Fint, olr)
       
       mat(:,i) = (flux_pert - flux)/pert_T
    end do

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
end module matrix
