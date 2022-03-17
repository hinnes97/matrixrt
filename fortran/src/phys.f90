module phys_mod
  use params, only : dp
  
  implicit none

  real(dp), parameter :: u_h2 = 2. 
  real(dp), parameter :: u_ch4 = 16.
  real(dp), parameter :: u_h2o = 18.
  real(dp), parameter :: u_nh3 = 17.
  real(dp), parameter :: u_co2 = 44.
  real(dp), parameter :: u_co  = 28.

  ! Give each of the molecules a unique integer

  integer, parameter :: i_h2 = 1
  integer, parameter :: i_ch4 = 2
  integer, parameter :: i_nh3 = 3
  integer, parameter :: i_co2 = 4
  integer, parameter :: i_co = 5
  integer, parameter :: i_h2o = 6
  
  real(dp), parameter, dimension(6) :: u_vec = (/u_h2,  u_ch4, u_nh3, u_co2, u_co, u_h2o /)
  
end module phys_mod
