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
#endif

  implicit none

  real(dp), dimension(:), allocatable :: Tf, pf ! Temp and pressure arrays
  real(dp), dimension(:), allocatable :: Te, pe     ! Edge pressure array
  real(dp), dimension(:), allocatable :: tau_IR, tau_V

  real(dp), dimension(:), allocatable :: net_F
  real(dp), dimension(:), allocatable :: dT
  real(dp) :: olr

  integer :: ncid

  integer :: i

  ! Initialise parameters and allocate arrays
  call read_constants()
  call allocate_arrays(Tf, pf, pe, tau_IR, tau_V, net_F, dT, Te)

#ifdef SOC
  call socrates_init()
#elif defined PICKET
  call radiation_init(nf)
#endif
  
  !Initialise output file
  call file_setup("output.nc", nf, ne, ncid)
  
  ! Initialise pressure and temperature arrays as log and lin spaced respectively
  call logspace(log_top_p, log_bot_p, pe)
  call linspace(top_t, bot_t, Tf)
  call linspace(top_t, bot_t, Te)
  
  ! Initialise pf array from pe
  do i=1,nf
     pf(i) = (pe(i+1) - pe(i)) / (log(pe(i+1)) - log(pe(i)))
  end do

!  do i=1,ne
!     Te(i) = 500._dp*(pe(i)/pe(ne))**(2._dp/7._dp)
!  end do
  
!  do i=1, nf
!     Te(i) = 500._dp*(pf(i)/pf(nf))**(2._dp/7._dp)
!  end do

!  do i=1,ne
!     if (Te(i) .lt. 200._dp) then
!        Te(i) = 200._dp
!     end if
!  end do
  
  do i=1,ne
     write(*,*) Te(i)
  end do
  
  ! Initialise tau arrays
  do i=1, ne
     tau_V(i) = tau_V_inf*(pe(i)/pe(ne))
     tau_IR(i)= tau_IR_inf*(pe(i)/pe(ne))
  end do

  ! Do timestepping
  !call step(Tf, pf, pe, tau_IR, tau_v, net_F, dT, olr)
  ! Do matrix method
  call do_matrix(nf, ne, Tf, pf, Te, pe, tau_IR, tau_V, 1.0_dp, Finc, Fint, olr)
  call dump_data(ncid, nf, ne, Tf, pf, pe, olr, tau_IR_inf, tau_V_inf, Finc, Fint,Te)
  call close_file(ncid)
  call deallocate_arrays(Tf, pf, pe, tau_IR, tau_V, net_F, dT,Te)
  
end program main
