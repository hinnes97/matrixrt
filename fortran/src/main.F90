program main
  !use radiation_Kitzmann_noscatt, only: Kitzmann_TS_noscatt
  !use, intrinsic  :: iso_fortran_env
  use io
  use utils
  use params
  use timestep
  use matrix

#ifdef SOC
  use socrates_interface_mod, only : socrates_init
#elif defined PICKET
  use radiation_mod, only : radiation_init
#elif defined TWOSTR
  use band_grey_mod, only : bg_init
#endif


  implicit none

  real(dp), dimension(:), allocatable :: Tf, pf ! Temp and pressure arrays
  real(dp), dimension(:), allocatable :: Te, pe     ! Edge pressure array
  real(dp), dimension(:), allocatable :: tau_IR, tau_V

  real(dp), dimension(:), allocatable :: net_F
  real(dp), dimension(:), allocatable :: dT
  real(dp), dimension(:), allocatable :: q, fdn, fup
  real(dp) :: olr

  integer :: ncid

  integer :: i

  ! Initialise parameters and allocate arrays
  call read_constants()
  call allocate_arrays(Tf, pf, pe, tau_IR, tau_V, net_F, dT, Te, q, fup, fdn)

#ifdef SOC
  call socrates_init()
#elif defined PICKET
  call radiation_init(nf)
#elif defined TWOSTR
  call bg_init(nf)
#endif

  !Initialise output file
  call file_setup(output_file, nf, ne, ncid)

  if (init_file .eqv. .true.) then
     call read_initial_data(input_file, Tf, Te)
     call logspace(log_top_p, log_bot_p, pe)
     ! Initialise pf array from pe
     do i=1,nf
        pf(i) = (pe(i+1) - pe(i)) / (log(pe(i+1)) - log(pe(i)))
     end do

  else
     ! Initialise pressure and temperature arrays as log and lin spaced respectively
     call logspace(log_top_p, log_bot_p, pe)
     ! Initialise pf array from pe
     do i=1,nf
        pf(i) = (pe(i+1) - pe(i)) / (log(pe(i+1)) - log(pe(i)))
     end do

     call linspace(top_t, bot_t, Tf)
     call linspace(top_t, bot_t, Te)
     write(*,*) (Tf(i), i=1,nf)
     write(*,*) (Te(i), i=1,nf+1)
     do i=1,nf
        if (Tf(i) .lt. 150) then
           Tf(i) = 150.
        endif
     enddo

     do i=1,nf+1
        if (Te(i) .lt. 150.) then
           Te(i) = 150.
        endif
     enddo
     
!!$
!!$  do i=1,ne
!!$     Te(i) = 600._dp*(pe(i)/pe(ne))**(2._dp/7._dp)
!!$  end do
!!$  
!!$  do i=1, nf
!!$     Te(i) = 600._dp*(pf(i)/pf(nf))**(2._dp/7._dp)
!!$  end do
!!$
!!$  do i=1,ne
!!$     if (Te(i) .lt. 150._dp) then
!!$        Te(i) = 150._dp
!!$     end if
!!$  end do
  
!!$     do i=1,ne
!!$        write(*,*) Te(i)
!!$     end do
  endif
  
  
  ! Initialise tau arrays
  do i=1, ne
     tau_V(i) = tau_V_inf*(pe(i)/pe(ne))
     tau_IR(i) = tau_IR_inf*(pe(i)/pe(ne))
     !tau_IR(i)= tau_IR_inf*(  0.8*(pe(i)/pe(ne)) +(1-0.8) * (pe(i)/pe(ne))**2 )
  end do

  if (matrix_rt) then
     ! Do matrix method
     call do_matrix(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, 1.0_dp, Finc, Fint, olr,q)
  else
     ! Do timestepping
     call step(Tf, pf, pe, tau_IR, tau_v, net_F, dT, olr, ncid,q, fup, fdn)

     ! Interpolate to Te
     do i=2,nf
        call linear_log_interp(pe(i), pf(i-1), pf(i), Tf(i-1), Tf(i), Te(i))
     enddo

     call linear_log_interp(pe(1), pf(1), pf(2), Tf(1), Tf(2), Te(1))
     call linear_log_interp(pe(ne), pf(nf-1), pf(nf), Tf(nf-1), Tf(nf), Te(ne) )

  endif
  
  call dump_data(ncid, nf, ne, Tf, pf, pe, olr, tau_IR_inf, tau_V_inf, Finc, Fint,Te, q, fup, fdn)
  call close_file(ncid)
  call deallocate_arrays(Tf, pf, pe, tau_IR, tau_V, net_F, dT,Te)
  
end program main
