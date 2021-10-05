! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to interpolate to the value of a field at a point.
!
! Method:
!	The field to be evaluated is passed to the routine, together
!	with its second derivatives if a splined fit is used. Its
!	value at the pressure supplied is evaluated according to the
!	type of interpolation prescribed. The result is returned.
!
!- ---------------------------------------------------------------------
      SUBROUTINE interpolate_p(ierr, n, p, a, x, y, y2, pp, aa
     &  , i_mode, l_splined)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE interp_mode_pcf
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments.
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      INTEGER, Intent(IN) ::
     &    n
!           Number of levels
     &  , i_mode
!           Mode of interpolation
      LOGICAL, Intent(INOUT) ::
     &     l_splined
!           True if field already splined
      REAL  (RealK), Intent(IN) ::
     &    p(n)
!           Pressure levels
     &  , a(n)
!           Field to be interpolated
     &  , pp
!           Interpolating pressure
      REAL  (RealK), Intent(INOUT) ::
     &    x(n)
!           Converted abscissa
     &  , y(n)
!           Converted ordinate
     &  , y2(n)
!           Second derivative of ordinate
      REAL  (RealK), Intent(OUT) ::
     &    aa
!           Interpolated value
!
!     Local variables.
      INTEGER
     &    i
!           Loop variable
      REAL  (RealK) ::
     &    xx
!           Converted evalaution point
     &  , yy
!           Raw interpolate
!
!     Subroutines called:
      EXTERNAL
     &    spline_fit, spline_evaluate
!
!
!
!     Perform the initial spline fit on the first call.
      IF (.NOT.l_splined) THEN
        IF ( (i_mode == IP_1_lin_lin).OR.
     &       (i_mode == IP_3_lin_lin) ) THEN
!         Linear-linear fit.
          DO i=1, n
            x(i)=p(i)
            y(i)=a(i)
          ENDDO
        ELSE IF ( (i_mode == IP_1_log_lin).OR.
     &            (i_mode == IP_3_log_lin) ) THEN
!         Logarithmic-linear fit.
          DO i=1, n
            x(i)=log(p(i))
            y(i)=a(i)
          ENDDO
        ELSE IF ( (i_mode == IP_1_lin_log).OR.
     &            (i_mode == IP_3_lin_log) ) THEN
!         Linear-logarithmic fit.
          DO i=1, n
            x(i)=p(i)
            y(i)=log(a(i))
          ENDDO
        ELSE IF ( (i_mode == IP_1_log_log).OR.
     &            (i_mode == IP_3_log_log) ) THEN
!         Logarithmic-logarithmic fit.
          DO i=1, n
            x(i)=log(p(i))
            y(i)=log(a(i))
          ENDDO
        ELSE
          WRITE(iu_err, '(/a)')
     &      '*** Error: Unrecognized fitting mode.'
          ierr=i_err_fatal
          RETURN
        ENDIF
!
        IF ( (i_mode == IP_1_lin_lin).OR.
     &       (i_mode == IP_1_log_lin).OR.
     &       (i_mode == IP_1_lin_log).OR.
     &       (i_mode == IP_1_log_log) ) THEN
!         This is a spline fit without second derivatives.
          DO i=1, n
            y2(i)=0.0_RealK
          ENDDO
        ELSE IF ( (i_mode == IP_3_lin_lin).OR.
     &            (i_mode == IP_3_log_lin).OR.
     &            (i_mode == IP_3_lin_log).OR.
     &            (i_mode == IP_3_log_log) ) THEN
!         Cubic splines.
          CALL spline_fit(n, x, y, y2)
        ENDIF
        l_splined=.true.
!
      ENDIF
!
!
!     Calculate the splining point from PP.
      IF ( (i_mode == IP_1_lin_lin).OR.
     &     (i_mode == IP_3_lin_lin) ) THEN
        xx=pp
      ELSE IF ( (i_mode == IP_1_log_lin).OR.
     &          (i_mode == IP_3_log_lin) ) THEN
        xx=log(pp)
      ELSE IF ( (i_mode == IP_1_lin_log).OR.
     &          (i_mode == IP_3_lin_log) ) THEN
        xx=pp
      ELSE IF ( (i_mode == IP_1_log_log).OR.
     &          (i_mode == IP_3_log_log) ) THEN
        xx=log(pp)
      ELSE
        WRITE(iu_err, '(/a)')
     &    '*** Error: Unrecognized fitting mode.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
      CALL spline_evaluate(ierr, n, x, y, y2, xx, yy)
      IF (ierr /= i_normal) THEN
        IF (ierr == i_err_range) THEN
          IF (xx < x(1)) THEN
            yy=y(1)
          ELSE IF (xx > x(n)) THEN
            yy=y(n)
          ENDIF
!         Recover form the error
          ierr=i_normal
        ELSE
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDIF
!
!     Convert back to the actual data field.
      IF ( (i_mode == IP_1_lin_lin).OR.
     &     (i_mode == IP_3_lin_lin) ) THEN
        aa=yy
      ELSE IF ( (i_mode == IP_1_log_lin).OR.
     &          (i_mode == IP_3_log_lin) ) THEN
        aa=yy
      ELSE IF ( (i_mode == IP_1_lin_log).OR.
     &          (i_mode == IP_3_lin_log) ) THEN
        aa=exp(yy)
      ELSE IF ( (i_mode == IP_1_log_log).OR.
     &          (i_mode == IP_3_log_log) ) THEN
        aa=exp(yy)
      ENDIF
!
!
!
      RETURN
      END
