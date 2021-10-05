! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the atmospheric density.
!
! Method:
!	The atmospheric density is calculated from T and q.
!
!- ---------------------------------------------------------------------
      FUNCTION fnc_density(p, t, q)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE rad_ccf, ONLY: r_gas_dry, ratio_molar_weight
!
!
      IMPLICIT NONE
!
!
!
!     Dummy variables
      REAL  (RealK), Intent(IN) ::
     &    p
!           Pressure
     &  , t
!           Temperature
     &  , q
!           Humidity
      REAL  (RealK) ::
     &    fnc_density
!           Calculated density
!
!
      fnc_density=p/(r_gas_dry*t
     &  *(1+(ratio_molar_weight-1.0_RealK)*q))
!
!
!
      RETURN
      END
