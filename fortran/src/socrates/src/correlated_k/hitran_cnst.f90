! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Constants for HITRAN molecules and isotopes
!
!- ---------------------------------------------------------------
MODULE hitran_cnst

  USE realtype_rd

  IMPLICIT NONE


  INTEGER, Parameter :: number_molecules = 31
  INTEGER, Parameter :: number_species   = 88
  INTEGER, Parameter :: number_temperatures = 931 ! 70 - 1000K for qcoeff

  CHARACTER (LEN=6), Parameter, DIMENSION(number_molecules) :: &
    molecule_names = &
       (/"   H2O", "   CO2", "    O3", "   N2O", "    CO", "   CH4", "    O2", & !  1-7
         "    NO", "   SO2", "   NO2", "   NH3", "  HNO3", "    OH", "    HF", & !  8-14
         "   HCl", "   HBr", "    HI", "   ClO", "   OCS", "  H2CO", "  HOCl", & ! 15-21
         "    N2", "   HCN", " CH3Cl", "  H2O2", "  C2H2", "  C2H6", "   PH3", & ! 22-28
         "  COF2", "   SF6", "   H2S" /)                                         ! 29-31

  INTEGER, Parameter, DIMENSION(1:number_molecules) :: number_isotopes = &
       (/ 6,10, 5, 5, 6, 4, 3, &
          3, 2, 1, 2, 1, 3, 1, &
          2, 2, 1, 2, 5, 3, 2, &
          1, 3, 2, 1, 3, 2, 1, &
          2, 1, 3 /)

  REAL  (RealK), Parameter, DIMENSION(1:number_species) :: q296 = &
     (/1.7464E+02,  1.7511E+02,  1.0479E+03,  8.5901E+02,  8.7519E+02,  5.2204E+03, & ! H20
       2.8694E+02,  5.7841E+02,  6.0948E+02,  3.5527E+03,  1.2291E+03,  7.1629E+03, & ! CO2
       3.2421E+02,  3.7764E+03,  1.1002E+04,  6.5350E+02, &
       3.4838E+03,  7.4657E+03,  3.6471E+03,  4.3331E+04,  2.1405E+04, &              ! O3
       5.0018E+03,  3.3619E+03,  3.4586E+03,  5.3147E+03,  3.0971E+04, &              ! N2O
       1.0712E+02,  2.2408E+02,  1.1247E+02,  6.5934E+02,  2.3582E+02,  1.3809E+03, & ! CO
       5.9045E+02,  1.1808E+03,  4.7750E+03,  1.5980E+03, &                           ! CH4
       2.1577E+02,  4.5230E+02,  2.6406E+03, &                                        ! O2
       1.1421E+03,  7.8926E+02,  1.2045E+03, &                                        ! NO
       6.3403E+03,  6.3689E+03, &                                                     ! SO2
       1.3578E+04, &                                                                  ! NO2
       1.7252E+03,  1.1527E+03, &                                                     ! NH3
       2.1412E+05, &                                                                  ! HNO3
       8.0362E+01,  8.0882E+01,  2.0931E+02, &                                        ! OH
       4.1466E+01, &                                                                  ! HF
       1.6066E+02,  1.6089E+02, &                                                     ! HCl
       2.0018E+02,  2.0024E+02, &                                                     ! HBr
       3.8900E+02, &                                                                  ! HI
       3.2746E+03,  3.3323E+03, &                                                     ! ClO
       1.2210E+03,  1.2535E+03,  2.4842E+03,  4.9501E+03,  1.3137E+03, &              ! OCS
       2.8467E+03,  5.8376E+03,  2.9864E+03, &                                        ! H2CO
       1.9274E+04,  1.9616E+04, &                                                     ! HOCl
       4.6598E+02, &                                                                  ! N2
       8.9529E+02, 1.8403E+03, 6.2141E+02, &                                          ! HCN
       5.7916E+04, 5.8834E+04, &                                                      ! CH3Cl
       9.8198E+03, &                                                                  ! H2O2
       4.1403E+02, 1.6562E+03, 1.5818E+03, &                                          ! C2H2
       7.0881E+04, 3.6191E+04, &                                                      ! C2H6
       3.2486E+03, &                                                                  ! PH3
       7.0044E+04, 3.7844E+04, &                                                      ! COF2
       1.6233E+06, &                                                                  ! SF6
       5.0307E+02, 5.0435E+02, 2.0149E+03 /)                                          ! H2S

                 
  REAL  (RealK), DIMENSION(1:number_temperatures,1:number_species) :: qcoeff
! Total internal partition sums for 70 <= T <=1000 K range, read from parsum.dat file
! using internal subroutine read_parsum_dat.

  REAL  (RealK), Parameter, DIMENSION(1:number_species) :: iso_mass = &
   (/18.010565, 20.014811, 19.014780, 19.016740, 21.020985, 20.020956, & ! H20
     43.989830, 44.993185, 45.994076, 44.994045, 46.997431, 45.997400, & ! CO2
     47.998322, 46.998291, 45.998262, 49.001675, &
     47.984745, 49.988991, 49.988991, 48.988960, 48.988960, &            ! O3
     44.001062, 44.998096, 44.998096, 46.005308, 45.005278, &            ! N2O
     27.994915, 28.998270, 29.999161, 28.999130, 31.002516, 30.002485, & ! CO
     16.031300, 17.034655, 17.037475, 18.040830, &                       ! CH4
     31.989830, 33.994076, 32.994045, &                                  ! O2
     29.997989, 30.995023, 32.002234, &                                  ! NO
     63.961901, 65.957695, &                                             ! SO2
     45.992904, &                                                        ! NO2
     17.026549, 18.023583, &                                             ! NH3
     62.995644, &                                                        ! HNO3
     17.002740, 19.006986, 18.008915, &                                  ! OH
     20.006229, &                                                        ! HF
     35.976678, 37.973729, &                                             ! HCl
     79.926160, 81.924115, &                                             ! HBr
    127.912297, &                                                        ! HI
     50.963768, 52.960819, &                                             ! ClO
     59.966986, 61.962780, 60.970341, 60.966371, 61.971231, &            ! OCS
     30.010565, 31.013920, 32.014811, &                                  ! H2CO
     51.971593, 53.968644, &                                             ! HOCl
     28.006147, &                                                        ! N2
     27.010899, 28.014254, 28.007933, &                                  ! HCN
     49.992328, 51.989379, &                                             ! CH3Cl
     34.005480, &                                                        ! H2O2
     26.015650, 27.019005, 27.021825, &                                  ! C2H2
     30.046950, 31.050305, &                                             ! C2H6
     33.997238, &                                                        ! PH3
     65.991722, 66.995083, &                                             ! COF2
     145.962492, &                                                       ! SF6
     33.987721, 35.983515, 34.987105 /)                                  ! H2S
     

  REAL  (RealK), PARAMETER :: &  ! Physical constants
       speed_of_light     = 2.997925E8,   & ! m s-1
       boltzmann_constant = 1.380622E-23, & ! J K-1
       planck_constant    = 6.626196E-34, & ! J s
       avogadro_number    = 6.022169E23,  & ! mol-1
       atomic_mass_unit   = 1.660531E-27, & ! kg
       molar_gas_constant = 8.31434E0,    & ! J mol-1 K-1
       c1                 = 2.0*planck_constant*speed_of_light * &
                                speed_of_light, & ! W m2
       c2                 = planck_constant*speed_of_light / &
                            boltzmann_constant, & ! m K
       PI                 = 3.1415926   

CONTAINS

  SUBROUTINE read_parsum_dat

    IMPLICIT NONE

    INTEGER :: ierr, iu_parsum, i, ios
    REAL :: temp
    CHARACTER (LEN=2917) :: line

    CALL get_free_unit(ierr, iu_parsum)
    CALL open_file_in(ierr, iu_parsum, 'Enter location of parsum.dat.')
    READ(iu_parsum, '(a)', iostat=ios) line
    DO i=1, number_temperatures
      READ(iu_parsum, '(f6.1,f20.6,999f27.6)', iostat=ios) temp, qcoeff(i,:)
    END DO
    CLOSE(iu_parsum)

  END SUBROUTINE read_parsum_dat

END MODULE hitran_cnst
