  MODULE twostr_module

  implicit none

  !private
  !public :: twostr_module

  CONTAINS

!     -----------------------------------------------------------
!                       HERE IS CALL_TWOSTR
!     -----------------------------------------------------------
!     =============================================================
!                        CALL_TWOSTR()
!     -------------------------------------------------------------
!     |* Call TWOSTR() from here, with input from RADTRAN() and
!     |* prescribed values.
!     |* The TWOSTR is described in Kylling, A., K. Stamnes and
!     |* S.-C. Tsay (1995). The RT solver can be dwonloaded from
!     |* the websit:
!     |* ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple_Scatt/
!     |* I slightly modify the TWOSTR to include a internal heat
!     |* source in the bottom of domain for giant planets condition,
!     |* which is characterized by Tint. By setting albedo = 1 and
!     |* emissivity=1, and temperature of lower boundary = Tint for
!     |* the lower boundary, the bottom condition is such that the
!     |* net upcoming flux = sigma*Tint^4.
!     |* Xianyu Tan, Aug. 2015
!      * Tan, March 2019: make the twostr subroutine double precision
!     -------------------------------------------------------------
!     |* Meaning of each variables is referred to documentation of
!     |* the TWOSTR code, available online:
!     |* ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple_Scatt/TWOSTR/

      SUBROUTINE CALL_TWOSTR(ii,jj,nlyr,temper,gg,ssalb,dtauc,ntau,utau,  &
          planck,wvnmlo,wvnmhi,Tint,fbeam,umu0,fnetup,flup)


      IMPLICIT NONE
      INTEGER maxcly  
      INTEGER maxulv 
      INTEGER nerr
      parameter (nerr=22)
      PARAMETER ( maxcly = 650, maxulv = 623) 
!     ------------------ INPUT VARIABLES -----------------
      LOGICAL planck
      INTEGER nlyr,ntau
      REAL*8 temper(0:maxcly),utau(maxulv)
      REAL*8 gg(maxcly),ssalb(maxcly),dtauc(maxcly)
      REAL*8 fbeam,umu0,Tint,wvnmlo,wvnmhi
!     ------------------ OUTPUT VARIABLES ----------------
      REAL*8 fnetup(maxulv)

!     ------------------ LOCAL VARIABLES -----------------
      LOGICAL :: quiet = .TRUE.
      LOGICAL :: spher = .false.
      LOGICAL :: deltam = .true.
      LOGICAL :: usrtau = .true.
      LOGICAL prnt(2)
      REAL*8 :: albedo = 0.0  ! 1 if use net upward flux
      REAL*8 :: fisot = 0.d0
      REAL*8 :: temis = 0.0
      REAL*8 :: ttemp = 0.0
      REAL*8 btemp,radius
      REAL*8 dfdt(maxulv), flup(maxulv), rfldir(maxulv), &
             rfldn(maxulv), uavg(maxulv),zd(0:maxcly)
      REAL*8 taucum
      INTEGER ierror(nerr),i,j,k,ierr,ii,jj
      CHARACTER header*127



!     INITIALIZATION FOR FNETUP
      DO i = 1,maxulv
         fnetup(i) = 0.d0
      ENDDO

      taucum = 0.
      do i = 1,nlyr
         taucum = taucum + dtauc(i)
      enddo
      if (utau(ntau).gt.taucum) utau(ntau) = taucum
      btemp = Tint
      prnt(1) = .false.
      prnt(2) = .false.

      CALL twostr(albedo, btemp, deltam, dtauc, fbeam, fisot,  &
           gg, header, ierror, maxcly, maxulv, nlyr, planck,  &
           ntau, prnt, quiet, radius, spher, ssalb, temis,  &
           temper, ttemp, umu0,  usrtau, utau, wvnmlo,  &
           wvnmhi, zd, dfdt, flup, rfldir, rfldn, uavg )
      
      do i = 1,ntau
         fnetup(i) = flup(i) - rfldir(i) - rfldn(i)
      enddo
      !IF ((ii == ni).AND.(jj == nj)) print*,'fnetup',fnetup(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'flup',flup(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'rfldir',rfldir(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'rfldn',rfldn(1:ntau)


         DO ierr = 1, nerr
            IF ( ierror(ierr) .NE. 0 ) THEN
               WRITE(*,'(/,A,I4,/)')  "TWOSTR REPORTS FATAL ERROR: ", &
                   ierr
            ENDIF
         ENDDO

      RETURN
      END SUBROUTINE
!     *****  END CALL_TWOSTR  *****


!     =============================================================
!     -------------------------------------------------------------
!         BELOW IS THE DISORT-TWOSTR RADIATIVE TRANSFER SOLVER
!     -------------------------------------------------------------
!     =============================================================

      SUBROUTINE twostr(albedo, btemp, deltam, dtauc, fbeam, fisot, &
           gg, header, ierror, maxcly, maxulv, nlyr, planck, &
           ntau, prnt, quiet, radius, spher, ssalb, temis, & 
           temper, ttemp, umu0,  usrtau, utau, wvnmlo, &
           wvnmhi, zd, dfdt, flup, rfldir, rfldn, uavg )
!***********************************************************************
! Copyright (C) 1993, 1994, 1995 Arve Kylling
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! To obtain a copy of the GNU General Public License write to the
! Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
! USA.
!
!***********************************************************************
!     twostr solves the radiative transfer equation in the two-stream
!     approximation. Based on the general-purpose algorithm
!     DISORT, but both simplified and extended as follows:
!
!     1) The only quantities calculated are mean intensities
!        and fluxes.
!         
!     2) The medium may be taken to be pseudo-spherical (spher=.TRUE.).
!
!     3) Only Lambertian reflection at the bottom boundary is allowed
!
!    ( See twostr.doc for complete documentation )
!
!***********************************************************************
!
!                 I/O variable specifications
!
      implicit none
      CHARACTER  header*127
      LOGICAL deltam, planck, prnt(2), quiet, spher, usrtau
      INTEGER ierror(22), maxcly, maxulv, nlyr, ntau
      REAL*8  albedo, btemp, fbeam, fisot, dtauc(maxcly), gg(maxcly),&
              radius, ssalb(maxcly), temper(0:maxcly ), temis, ttemp, &
              wvnmlo, wvnmhi, umu0, utau( maxulv ), zd( 0:maxcly )
!
      REAL*8    rfldir(maxulv), rfldn(maxulv), flup(maxulv), &
              uavg( maxulv ), dfdt( maxulv )
      !
!+---------------------------------------------------------------------+
!      Routines called (in order):  tzeroal, tchekin, settwo, tprtinp,
!                                   hopsol, setmtx, solbcd, tfluxes,
!+---------------------------------------------------------------------+
!
!  Index conventions (for all do-loops and all variable descriptions):
!
!     iq     :  For quadrature angles
!
!     lu     :  For user levels
!
!     lc     :  For computational layers (each having a different
!               single-scatter albedo and/or phase function)
!
!     lev    :  For computational levels
!
!     ls     :  Runs from 0 to 2*mxcly+1, ls = 1,2,3 refers to top,
!               center and bottom of layer 1, ls = 3,4,5 refers to
!               top, center and bottom of layer 2, etc.....
!
!+---------------------------------------------------------------------+
!               I n t e r n a l    v a r i a b l e s
!
!   b()               Right-hand side vector of Eqs. KST(38-41), set 
!                     in *solbcd*
!
!   bplank            Intensity emitted from bottom boundary
!
!   cband()           Matrix of left-hand side of the linear system
!                     Eqs. KST(38-41);  in tridiagonal form
!
!   ch(lc)            The Chapman-factor to correct for pseudo-
!                     spherical geometry in the direct beam.
!
!   chtau(lc)         The optical depth in spherical geometry.
!
!   cmu               Computational polar angle, single or double
!                     Gaussian quadrature rule used, see routine
!                     -settwo-
!
!   expbea(lc)        Transmission of direct beam in delta-m optical
!                     depth coordinates
!
!   fldn(lu)          Diffuse down flux (delta-m scaled)
!
!   fldir(lu)         Direct beam flux (delta-m scaled)
!
!   flyr(lc)          Truncated fraction in delta-m method
!
!   kk(lc)            Eigenvalues in Eq. KST(20)
!
!   layru(lu)         Computational layer in which user output level
!                     -utau(lu)- is located
!
!   ll(iq,lc)         Constants of integration C-tilde in Eqs.
!   KST(42-43)
!                     obtained by solving Eqs. KST(38-41)
!
!   lyrcut            True, radiation is assumed zero below layer
!                     -ncut- because of almost complete absorption
!
!   ncut              Computational layer number in which absorption
!                     optical depth first exceeds -abscut-
!
!   oprim(lc)         Single scattering albedo after delta-m scaling
!
!   pass1             True on first entry, false thereafter
!
!   pkag(0:lc)        Integrated Planck function for internal emission
!                     at layer boundaries
!
!   pkagc(lc)         Integrated Planck function for internal emission
!                     at layer center
!
!   rr(lc)            Eigenvectors at polar quadrature angles.
!
!   tauc(0:lc)        Cumulative optical depth (un-delta-m-scaled)
!
!   taucpr(0:lc)      Cumulative optical depth (delta-m-scaled if
!                     deltam = true, otherwise equal to -tauc-)
!
!   tplank            Intensity emitted from top boundary
!
!   u0c(iq,lu)        Azimuthally-averaged intensity
!
!   utaupr(lu)        Optical depths of user output levels in delta-m
!                     coordinates;  equal to  -utau(lu)- if no delta-m
!
!   xb_0d(lc)         x-sub-zero-sup-minus in expansion of pseudo 
!                     spherical beam source, Eq. KST(22)
!
!   xb_0u(lc)         x-sub-zero-sup-plus in expansion of pseudo 
!                     spherical beam source, Eq. KST(22)
!
!   xb_1d(lc)         x-sub-one-sup-minus in expansion of pseudo 
!                     spherical beam source, Eq. KST(22)
!
!   xb_1u(lc)         x-sub-one-sup-plus in expansion of pseudo 
!                     spherical beam source, Eq. KST(22)
!
!   xp_0(lc)          x-sub-zero in expansion of thermal source func-
!                     tion; see Eq. KST(22) (has no (mu) dependence)
!
!   xp_1(lc)          x-sub-one in expansion of thermal source func-
!                     tion; see Eq. KST(22) (has no (mu) dependence)
!
!   yb_0d(lc)         y-sub-zero-sup-minus in Eq. KST(23), solution for
!                     pseudo spherical beam source.
!
!   yb_0u(lc)         y-sub-zero-sup-plus in Eq. KST(23), solution for
!                     pseudo spherical beam source.
!
!   yb_1d(lc)         y-sub-one-sup-minus in Eq. KST(23), solution for
!                     pseudo spherical beam source.
!
!   yb_1u(lc)         y-sub-one-sup-plus in Eq. KST(23), solution for
!                     pseudo spherical beam source.
!
!   yp_0d(lc)         y-sub-zero-sup-minus in Eq. KST(23), solution for
!                     thermal source.
!
!   yp_0u(lc)         y-sub-zero-sup-plus in Eq. KST(23), solution for
!                     thermal source.
!
!   yp_1d(lc)         y-sub-one-sup-minus in Eq. KST(23), solution for
!                     thermal source.
!
!   yp_1u(lc)         y-sub-one-sup-plus in Eq. KST(23), solution for
!                     thermal source.
!
!   zb_a(lc)          Alfa coefficient in Eq.  KST(22) for pseudo-
!                     spherical beam source.
!
!   zp_a(lc)          Alfa coefficient in Eq. KST(22) for thermal
!                     source.
!
!+---------------------------------------------------------------------+
!   Local symbolic dimensions:
!
!       mxcly  = max no. of computational layers
!       mxulv  = max no. of output levels
!       mxcmu  = max no. of computation polar angles (=2 in two-stream)
!+---------------------------------------------------------------------+
      integer mxcly,mxulv,mxcmu,mi,mi9m2,nnlyri
      PARAMETER ( mxcly = 650, mxulv = 651, mxcmu = 2, mi = mxcmu/2, &
                   mi9m2 = 4, nnlyri = mxcly*mxcmu )
!
      LOGICAL lyrcut, pass1
      INTEGER ipvt(nnlyri), iret, layru(mxulv), nerr
      REAL*8 b(nnlyri),cband(mi9m2,nnlyri),ch(mxcly),chtau(2*mxcly+1), &
           cmu, dtaucpr(mxcly), expbea(0:mxcly), &
           flyr(mxcly), fldn(mxulv), fldir(mxulv), ggprim(mxcly), &
           kk(mxcly),  ll(mxcmu,mxcly), oprim( mxcly ), pkag(0:mxcly), &
           pkagc(mxcly), rr(mxcly), tauc(0:mxcly), & 
           taucpr(0:mxcly), u0c(mxcmu,mxulv), utaupr(mxulv), &
           xb_0d(mxcly), xb_0u(mxcly), xb_1d(mxcly), xb_1u(mxcly), &
           xp_0(mxcly), xp_1(mxcly), yb_0d(mxcly), yb_0u(mxcly), &
           yb_1d(mxcly), yb_1u(mxcly), yp_0d(mxcly), yp_0u(mxcly), &
           yp_1d(mxcly), yp_1u(mxcly), zb_a(mxcly), &
           zp_a(mxcly), subd(2*mxcly), &
           diag(2*mxcly), superd(2*mxcly), sol(2*mxcly) 
!
      real*8 biggest,smallest, bplank,epsil,pi,tplank
      integer i,ierr,lc,ncol,ncut,nn,nstr

      SAVE  pass1, pi, epsil
      DATA  pass1 / .TRUE. /, nerr / 22 /
      !

      IF ( pass1 )  THEN
         pi = 2. * ASIN(1.d0)
         epsil =  0.00001             ! Typical value for 32 bits machine
         pass1 = .FALSE.
      END IF
!
      IF ( prnt(1) )  WRITE( *,1010 )  header
!
! Zero some arrays (not strictly necessary, but otherwise unused 
! parts of arrays collect garbage)
!
      DO i = 1, nerr
         ierror(i) = 0
      ENDDO
      CALL tzeroal( ch, kk, ll, mxcmu, mxcly, &
           nnlyri, rr, xb_0d, xb_0u, xb_1d, xb_1u, &
           xp_0, xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, &
           yp_0u, yp_1d, yp_1u, zb_a, zp_a )
!
! Calculate cumulative optical depth and dither single-scatter
! albedo to improve numerical behavior of eigenvalue/vector
! computation
!
      tauc( 0 ) = 0.
      CALL  zeroit( tauc(0), mxcly+1 )
      DO 20  lc = 1, nlyr
         IF( ssalb(lc) .GT. 1.0 - epsil )  ssalb(lc) = 1.0 - epsil
         tauc(lc) = tauc(lc-1) + dtauc(lc)
 20   CONTINUE
!
! Check input dimensions and variables
!
      CALL tchekin( albedo, btemp, dtauc, fbeam, fisot, gg, &
           ierror, maxcly, maxulv, mxcly, mxulv, nlyr, ntau, planck, &
           quiet, spher, ssalb, tauc, temis, temper, ttemp, &
           umu0, usrtau, utau, wvnmlo, wvnmhi, zd )
!
      iret = 0
      DO ierr = 1, nerr
         IF ( ierror(ierr) .NE. 0 ) THEN
            iret = 1
            IF ( .NOT. quiet )  THEN
               WRITE(*,'(/,A,I4,/)')  "TWOSTR REPORTS FATAL ERROR: ", &
                   ierr
            ENDIF
         ENDIF
      ENDDO
      IF ( iret .EQ. 1 ) RETURN
!     
! Perform various setup operations
!
      CALL settwo( albedo, biggest, bplank, btemp, ch, chtau, cmu, &
           deltam,  &
           dtauc, dtaucpr, expbea, fbeam, flyr, gg, ggprim, &
           layru, lyrcut, mxcly, ncut, nlyr, nn, nstr, ntau, oprim, &
           pkag, pkagc, planck, radius, smallest, spher, ssalb, tauc, & 
           taucpr, temis, temper, tplank, ttemp, umu0, usrtau, &
           utau, utaupr, wvnmlo, wvnmhi, zd )
!     
! Print input information
!
      IF ( prnt(1) ) &
           CALL tprtinp( albedo, btemp, deltam, dtauc, fbeam, fisot, &  
           flyr, gg, lyrcut, nlyr, planck, nstr, ntau, oprim, &
           spher, ssalb, tauc, taucpr, temis, temper, ttemp, umu0, &
           utau, wvnmlo, wvnmhi )
!
! Calculate the homogenous and particular solutions
!
      CALL hopsol(biggest, ch, chtau, cmu, fbeam, ggprim, kk, ncut, &
           nlyr,oprim, pi, pkag, pkagc, planck, radius, rr, smallest, &
           spher,taucpr, umu0, xb_0d, xb_0u, xb_1d, xb_1u, xp_0,xp_1, &
           yb_0d,yb_0u,yb_1d,yb_1u, yp_0d,yp_0u, yp_1d,yp_1u,zb_a, &
           zp_a )
!     
! Solve for constants of integration in homogeneous solution
! (general boundary conditions)
!
      CALL solbcd( albedo, b, bplank, cband, cmu, diag, expbea, &
           fbeam, fisot, ipvt, kk, ll, lyrcut, mi, mi9m2, mxcmu, &
           ncol, ncut, nlyr, nn, nnlyri, nstr, pi, rr, subd, &
           superd, sol, tplank, taucpr, umu0, yb_0d, yb_0u, yb_1d, &
           yb_1u, yp_0d, yp_0u, yp_1d, yp_1u, zb_a,zp_a )
!
! Compute upward and downward fluxes, mean intensities and
! flux divergences.
!
      CALL tfluxes(ch, cmu, dfdt, fbeam, flup, fldn, fldir, kk, &
           layru, ll, lyrcut, maxulv, mxcmu, mxulv, ncut, ntau, &
           pi, planck, prnt, oprim, rfldir, rfldn, rr, spher, ssalb, &
           taucpr, u0c, uavg, umu0, utau, utaupr, xp_0, xp_1, &
           yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, yp_1d, yp_1u, &
           zb_a, zp_a )
!     
      RETURN
!
1010  FORMAT ( ////, 1X, 120('*'), /, 25X, &
        'Two stream method radiative transfer program, version 1.13', &
        /, 1X, A, /, 1X, 120('*') )
!
      END SUBROUTINE
!
      REAL*8 FUNCTION chapmn( lc, taup, tauc, nlyr, zd, dtauc, &
           zenang, r )
!
! Calculates the Chapman-factor
!
! I n p u t       v a r i a b l e s:
!
!      lc        : Computational layer
!      nlyr      : Number of layers in atmospheric model
!      zd(lc)    : lc = 0, nlyr. zd(lc) is distance from bottom 
!                  surface to top of layer lc. zd(nlyr) = 0.0 km
!      dtauc     : Optical thickness of layer lc (un-delta-m-scaled)
!      zenang    : Solar zenith angle as seen from bottom surface
!      r         : Radius of earth. NOTE: Use the same dimension as zd,
!                  for instance both in km.
!
! O u t p u t      v a r i a b l e s:
!
!      ch        : Chapman-factor. in a pseudo-spherical atmosphere 
!                  replace EXP( -tau/umu0 ) by EXP( -ch(lc) ) in the
!                  beam source in 
!
! I n t e r n a l     v a r i a b l e s:
!
!      dhj       : delta-h-sub-j in Eq. B2 (DS)
!      dsj       : delta-s-sub-j in Eq. B2 (DS)
!      fact      : =1 for first sum in Eq. B2 (DS)
!                  =2 for second sum in Eq. B2 (DS)
!      rj        : r-sub-j in Eq. B1 (DS)
!      rjp1      : r-sub-j+1 in Eq. B1 (DS)
!      xpsinz    : The length of the line OG in Fig. 1, (DS)
!
!
      implicit none
      integer lc,nlyr,id,j
      real*8 pi,rj,rjp1,sum,xp,xpsinz,zenrad
      real*8 taup, zenang, r, dhj, dsj, fact, fact2
      REAL*8 dtauc(*) ,tauc(0:*), zd(0:*)
!
      pi     = 2.0 * ASIN( 1.d0 )
      zenrad = zenang * (2. * ASIN(1.d0)) / 180.0
      xp     = r +  zd(lc) + (zd(lc-1) - zd(lc) ) * &
          ( tauc(lc) - taup ) / dtauc(lc)
      xpsinz = xp * SIN( zenrad )
!
      IF( (zenang.GT.90.0) .AND. (xpsinz.LT.r) ) THEN
        chapmn = 1.0E+20
        RETURN
      END IF
!
! Find index of layer in which the screening height lies
!
      id = lc
      IF( zenang.GT.90.0 ) THEN
         DO 100 j = lc, nlyr
            IF( (xpsinz.LT.( zd(j-1) + r ) ) .AND. &
                              (xpsinz.GE.( zd(j) + r )) ) id = j
 100     CONTINUE
      END IF
!
      sum = 0.0
      DO 200 j = 1, id
        fact = 1.0
        fact2= 1.0
!
! Include factor of 2 for zenang .gt. 90, second sum in Eq. B2 (DS)
!
        IF( j.GT.lc ) fact = 2.0
        IF(j.EQ.id .AND. id.EQ.lc .AND. zenang.GT.90.0) fact2 = -1.0
!
        rj = r + zd(j-1)
        rjp1 = r + zd(j)
        IF(j.EQ.lc .AND. id.EQ.lc) rjp1 = xp
!
        dhj = zd(j-1) -zd(j)
        IF(id.GT.lc .AND. j.EQ.id) THEN
           dsj = SQRT(rj*rj - xpsinz*xpsinz )
        ELSE
           dsj = SQRT( rj*rj - xpsinz*xpsinz ) - &
                fact2 * SQRT( rjp1*rjp1 - xpsinz*xpsinz )
        END IF
!
        sum = sum + dtauc(j)*fact* dsj / dhj
!
 200  CONTINUE
!
! Add third term in Eq. B2 (DS)
!
      IF( id.GT.lc ) THEN
        dhj = zd(lc-1) -zd(lc)
        dsj = SQRT( xp*xp - xpsinz*xpsinz ) - &
             SQRT( (zd(lc)+r)*(zd(lc)+r) - xpsinz*xpsinz )
        sum = sum + dtauc(lc) * dsj / dhj
      END IF
!
      chapmn = sum
      RETURN
      END FUNCTION
!
      SUBROUTINE  tchekin(albedo, btemp, dtauc, fbeam, fisot, gg, &
           ierror, maxcly, maxulv, mxcly, mxulv, nlyr, ntau, planck, &
           quiet, spher, ssalb, tauc, temis, temper, ttemp, umu0, &
           usrtau, utau, wvnmlo, wvnmhi, zd )
!
! Checks the input dimensions and variables
!
      implicit none
      integer lc,lu
      real*8 umumin
      LOGICAL inperr, planck, quiet, spher, usrtau!, wrtbad, wrtdim
      INTEGER  maxcly, maxulv, mxcly, mxulv, nlyr, ntau, ierror(*)
      REAL*8 albedo, btemp, dtauc(*), fbeam, fisot, gg(*), &
           ssalb(*), tauc(0:*), temis, temper(0:*), ttemp, &
           umu0, utau(*), wvnmlo, wvnmhi, zd(0:*)
!
      inperr = .FALSE.
      IF ( nlyr.LT.1 ) THEN
         inperr = wrtbad(quiet, 'nlyr' )
         ierror(1) = 1
      ENDIF
      IF ( nlyr.GT.maxcly ) THEN
         inperr = wrtbad(quiet, 'maxcly' )
         ierror(2) = 1
      ENDIF
!
      DO 10  lc = 1, nlyr
         IF ( dtauc(lc).LT.0.0 )  THEN
            inperr = wrtbad( quiet, 'dtauc' )
            ierror(3) = ierror(3) + 1
            print*,'dtauc,lc',dtauc(lc),lc
            !print*,'dtauc',dtauc(1:nlyr),
            print*, 'temperature', temper(lc-1), temper(lc), temper(lc+1)
            
         ENDIF
         IF ( ssalb(lc).LT.0.0 .OR. ssalb(lc).GT.1.0 ) THEN
            inperr = wrtbad( quiet, 'ssalb' )
            ierror(4) = ierror(4) + 1
         ENDIF
         IF ( planck )  THEN
            IF( lc.EQ.1 .AND. temper(0).LT.0.0 ) THEN
               inperr = wrtbad( quiet, 'temper' )
               ierror(5) = ierror(5) + 1
            ENDIF
            IF( temper(lc).LT.0.0 ) THEN 
               inperr = wrtbad( quiet, 'temper' )
               ierror(5) = ierror(5) + 1
            ENDIF
         ENDIF
         IF( gg(lc).LT.-1.0 .OR. gg(lc).GT.1.0 ) THEN
            inperr = wrtbad( quiet, 'gg' )
            ierror(6) = ierror(6) + 1
         ENDIF
 10   CONTINUE
      IF ( spher ) THEN
         DO 11 lc = 1, nlyr
            IF ( zd(lc) .GT. zd(lc-1) ) THEN
               inperr = wrtbad( quiet, 'zd' )
               ierror(7) = ierror(7) + 1
            ENDIF
 11      CONTINUE          
      ENDIF
!
      IF ( usrtau )  THEN
         IF ( ntau.LT.1 ) THEN
            inperr = wrtbad( quiet, 'ntau' )
            ierror(8) = 1
         ENDIF
         IF ( maxulv.LT.ntau ) THEN
            inperr = wrtbad( quiet, 'maxulv' )
            ierror(9) = 1
         ENDIF
         DO 20  lu = 1, ntau
            IF( ABS(utau(lu)-tauc(nlyr)).LE.1.E-4) utau(lu) = tauc(nlyr)
            IF( utau(lu).LT.0.0 .OR. utau(lu).GT.tauc(nlyr) ) THEN
               inperr = wrtbad( quiet, 'utau' )
               ierror(10) = ierror(10) + 1
               
               print*, utau(1:lu)
               print*, '-------------------------------'
               print*, tauc(1:nlyr)
            ENDIF
 20      CONTINUE
      ELSE
         IF ( maxulv.LT.nlyr+1 ) THEN
            inperr = wrtbad( quiet, 'maxulv' )
            ierror(11) = 1
         ENDIF
      END IF
!
      IF ( fbeam.LT.0.0 ) THEN
         inperr = wrtbad( quiet, 'fbeam' )
         ierror(12) = 1
      ENDIF
      umumin = 0.0
      IF ( spher ) umumin = - 1.0
      IF ( fbeam.GT.0.0 .AND. ( umu0.LE.umumin .OR. umu0.GT.1.0 ) ) &
           THEN
         inperr = wrtbad( quiet, 'umu0' )
         ierror(13) = 1
      ENDIF
      IF ( fisot.LT.0.0 ) THEN
         inperr = wrtbad( quiet, 'fisot' )
         ierror(14) = 1
      ENDIF
      IF ( albedo.LT.0.0 .OR. albedo.GT.1.0 ) THEN
         inperr = wrtbad( quiet, 'albedo' )
         ierror(15) = 1
      ENDIF
!
      IF ( planck )  THEN
         IF ( wvnmlo.LT.0.0 .OR. wvnmhi.LT.wvnmlo ) THEN
            inperr = wrtbad( quiet, 'wvnmlo,hi' )
            ierror(16) = 1
         ENDIF
         IF ( temis.LT.0.0 .OR. temis.GT.1.0 ) THEN
            inperr = wrtbad( quiet, 'temis' )
            ierror(17) = 1
         ENDIF
         IF ( btemp.LT.0.0 ) THEN
            inperr = wrtbad( quiet, 'btemp' )
            ierror(18) = 1
         ENDIF
         IF ( ttemp.LT.0.0 ) THEN
            inperr = wrtbad( quiet, 'ttemp' )
            ierror(19) = 1
         ENDIF
      END IF
!
      IF ( mxcly.LT.nlyr ) THEN
         inperr = wrtdim( quiet, 'mxcly', nlyr )
         ierror(20) = 1
      ENDIF
      IF ( usrtau .AND. mxulv.LT.ntau ) THEN
         inperr = wrtdim( quiet, 'mxulv', ntau )
         ierror(21) = 1
      ENDIF
      IF ( .NOT.usrtau .AND. mxulv.LT.nlyr+1 ) THEN
         inperr = wrtdim( quiet, 'mxulv', nlyr+1 )
         ierror(22) = 1
      ENDIF
!
      IF ( inperr .AND. .NOT. quiet ) &
         CALL errmsg( 'twostr--input and/or dimension errors', .TRUE. )
!
      DO 100  lc = 1, nlyr
         IF ( (planck .AND. ABS(temper(lc)-temper(lc-1)) .GT. 50.0) &
                .AND. .NOT. quiet ) &
                CALL errmsg( 'chekin--vertical temperature step may' &
                // ' be too large for good accuracy', .FALSE. )
100   CONTINUE
!
      RETURN
      END SUBROUTINE
!
      SUBROUTINE tfluxes(ch, cmu, dfdt, fbeam, flup, fldn, fldir, kk, &
           layru, ll, lyrcut, maxulv, mxcmu, mxulv, ncut, ntau, &
           pi, planck, prnt, oprim, rfldir, rfldn, rr, spher, ssalb, &
           taucpr, u0c, uavg, umu0, utau, utaupr, xp_0, xp_1, &
           yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, yp_1d, yp_1u, &
           zb_a, zp_a )
!
! Calculates the radiative fluxes, mean intensity, and flux
! derivative with respect to optical depth from the 
! azimuthally-averaged intensity
!
! I n p u t     v a r i a b l e s:
!
!       ch       :  Chapman factor
!       cmu      :  Abscissae for gauss quadrature over angle cosine
!       kk       :  Eigenvalues 
!       layru    :  Layer number of user level -utau-
!       ll       :  Constants of integration in Eqs. KST(42-43),
!       obtained
!                   by solving Eqs. KST(38-41)
!       lyrcut   :  Logical flag for truncation of comput. layer
!       ncut     :  Number of computational layer where absorption
!                     optical depth exceeds -abscut-
!       oprim    :  Delta-m scaled single scattering albedo
!       rr       :  Eigenvectors at polar quadrature angles 
!       taucpr   :  Cumulative optical depth (delta-m-scaled)
!       utaupr   :  Optical depths of user output levels in delta-m
!                     coordinates;  equal to  -utau- if no delta-m
!       xp_0,    :  Thermal source coefficients x-sub-zero and
!        xp_1         x-sub-one in Eq. KST(22)
!       yb_0d,u, :  Thermal source vectors, Eq. KST(23)
!        yb_1d,u     
!       yp_0d,u, :  Beam source vectors, Eq. KST(23)
!        yp_1d,u     
!       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
!       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
!       (remainder are 'twostr' input variables)
!
! O u t p u t     v a r i a b l e s:
!
!       u0c      :  Azimuthally averaged intensities at polar
!                     quadrature angle cmu
!       (rfldir, rfldn, flup, dfdt, uavg are 'twostr' output variables)
!
! I n t e r n a l       v a r i a b l e s:
!
!       dirint   :  direct intensity attenuated
!       fdntot   :  total downward flux (direct + diffuse)
!       fldir    :  direct-beam flux (delta-m scaled)
!       fldn     :  diffuse down-flux (delta-m scaled)
!       fnet     :  net flux (total-down - diffuse-up)
!       fact     :  EXP( - utaupr / ch ), where ch is the Chapman factor
!       plsorc   :  Planck source function (thermal)
!+---------------------------------------------------------------------+


      implicit none
      INTEGER layru(*)
      integer mxcmu,maxulv,mxulv,ncut,ntau,iq,lu,lyu
      LOGICAL lyrcut, planck, prnt(*), spher
      REAL*8 ch(*), cmu, dfdt(*), flup(*), fldir(*), fldn(*), kk( * ), &
           ll( mxcmu,* ), oprim(*), rfldir(*), rfldn(* ), rr(*), &
           ssalb(*), taucpr( 0:* ), u0c( mxcmu,mxulv ), uavg(*), &
           utau(*), utaupr(*), xp_0(*), xp_1(*), yb_0d(*), yb_0u(*), &
           yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*), yp_1u(*), &
           zb_a(*), zp_a(*)
!
      real*8 fbeam,umu0,pi,dirint,fact,fdntot,fnet,plsorc
      IF ( prnt(2) )  WRITE( *,1010 )
!
! Zero twostr output arrays
!
      CALL  zeroit( u0c, mxulv*mxcmu )
      CALL  zeroit( rfldir, maxulv )
      CALL  zeroit( fldir,  mxulv )
      CALL  zeroit( rfldn,  maxulv )
      CALL  zeroit( fldn,   mxulv )
      CALL  zeroit( flup,   maxulv )
      CALL  zeroit( uavg,   maxulv )
      CALL  zeroit( dfdt,   maxulv )
!
! Loop over user levels
!

      !IF ((ii == ni).AND.(jj == nj)) print*,'planck,fbeam',planck,fbeam
      !IF ((ii == ni).AND.(jj == nj)) print*,'utaupr',utaupr(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'taucpr',taucpr(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'-zp_a',-zp_a(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'yp_0u',yp_0d(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'yp_1d',yp_1d(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'yp_0u',yp_0u(1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'yp_1u',yp_1u(1:ntau)
      IF ( planck ) THEN
         DO 100  lu = 1, ntau
            lyu = layru(lu)
            iq = 1
            u0c( iq,lu ) = u0c( iq,lu ) + &
                DEXP(-zp_a(lyu)*utaupr(lu))* &
                (yp_0d(lyu)+yp_1d(lyu)*utaupr(lu))
            iq = 2
            u0c( iq,lu ) = u0c( iq,lu ) + &
                DEXP(-zp_a(lyu)*utaupr(lu))* &
                (yp_0u(lyu)+yp_1u(lyu)*utaupr(lu))
 100     CONTINUE
      ENDIF
!
      DO 102  lu = 1, ntau
!
         lyu = layru(lu)
         IF ( lyrcut .AND. lyu.GT.ncut ) THEN
!
! No radiation reaches this level
!
            fdntot = 0.0
            fnet   = 0.0
            plsorc = 0.0
            GO TO 90                    ! ** Done with this layer
         END IF
!
         IF ( fbeam.GT.0.0 ) THEN
            iq = 1
            u0c( iq,lu ) = u0c(iq,lu) + &
                 DEXP(-zb_a(lyu)*utaupr(lu)) * &
                 (yb_0d(lyu)+yb_1d(lyu)*utaupr(lu) )
            iq = 2
            u0c( iq,lu ) = u0c(iq,lu) + &
                DEXP(-zb_a(lyu)*utaupr(lu))* &
                ( yb_0u(lyu)+yb_1u(lyu)*utaupr(lu) )
            IF ( umu0.GT.0.0 .OR. spher  ) THEN
               fact         = DEXP( - utaupr(lu) / ch(lyu) )
               dirint       = fbeam * fact
               fldir(  lu ) = ABS(umu0) * ( fbeam * fact )
               rfldir( lu ) = ABS(umu0)*fbeam*DEXP( -utau(lu)/ch(lyu) )
            ELSE
               dirint       = 0.0
               fldir(  lu ) = 0.0
               rfldir( lu ) = 0.0
            ENDIF
         ELSE
            dirint       = 0.0
            fldir(  lu ) = 0.0
            rfldir( lu ) = 0.0
         END IF
!
         iq = 1
         u0c( iq,lu ) =  u0c( iq,lu ) + rr(lyu) * ll(1,lyu) * &
                      DEXP( - kk(lyu) * (taucpr(lyu) - utaupr(lu)) ) &
                    + ll(2,lyu) * &
                      DEXP( - kk(lyu) * (utaupr(lu) - taucpr(lyu-1)) )
         iq = 2
         u0c(iq,lu ) = u0c( iq,lu ) + ll(1,lyu) * &
                      DEXP( - kk(lyu) * (taucpr(lyu) - utaupr(lu) ) ) &
                    + rr(lyu) * ll(2,lyu) * &
                      DEXP( - kk(lyu) * (utaupr(lu) - taucpr(lyu-1)) )
!
! Calculate fluxes and mean intensities
!
! Downward and upward fluxes from Eq. KST(9)
!
         fldn( lu )  = 2.0 * 2. * ASIN(1.d0) * cmu * u0c(1,lu)
         flup( lu )  = 2.0 * 2. * ASIN(1.d0) *  cmu * u0c( 2,lu )
         fdntot      = fldn( lu ) + fldir( lu )
         fnet        = fdntot - flup( lu )
         rfldn( lu ) = fdntot - rfldir( lu )
!
! Mean intensity from Eq. KST(10)
!
         uavg(lu)   = u0c( 1,lu ) + u0c( 2, lu )
         uavg( lu ) = (2.0*2.*ASIN(1.d0)*uavg(lu)+dirint) / (4.*2.*ASIN(1.d0))
!
! Flux divergence from Eqs. KST(11-12)
!
         plsorc = (1./(1.-oprim(lyu)))*DEXP(-zp_a(lyu)*utaupr(lu))* &
                                (xp_0(lyu) + xp_1(lyu)* utaupr(lu))
         dfdt( lu ) =  (1.0-ssalb(lyu)) * 4.*2.*ASIN(1.d0)* (  uavg(lu) - plsorc )
!
 90      IF( prnt(2) )  WRITE( *,1020 ) utau(lu), lyu, rfldir(lu), &
                                       rfldn(lu), fdntot, flup(lu), &
                                       fnet, uavg(lu), plsorc, dfdt(lu)
 102  CONTINUE
!
 1010 FORMAT( //, 21X, &
       '<----------------------- Fluxes ----------------------->', /, &
       '   optical  compu    downward    downward    downward     ', &
       ' upward                    mean      Planck   d(net flux)', /, &
       '     depth  layer      direct     diffuse       total     ', &
       'diffuse         net   intensity      source   / d(op dep)', / )
 1020 FORMAT( F10.4, I7, 1P,7E12.3, E14.3 )
!
      !IF ((ii == ni).AND.(jj == nj)) print*,'u0c 1',u0c(1,1:ntau)
      !IF ((ii == ni).AND.(jj == nj)) print*,'u0c 2',u0c(2,1:ntau)
      RETURN
      END SUBROUTINE
!
      SUBROUTINE hopsol(biggest, ch, chtau, cmu, fbeam, ggprim, kk, &
           ncut,nlyr, oprim, pi, pkag, pkagc, planck, radius, rr, &
           smallest,spher,taucpr, umu0, xb_0d, xb_0u, xb_1d, xb_1u, &
           xp_0,xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, &
           yp_1d, yp_1u, zb_a, zp_a )
!
! Calculates the homogenous and particular solutions to the
! radiative transfer equation in the two-stream approximation,
! for each layer in the medium.
!
!    I n p u t     v a r i a b l e s:
!
!       ch       :  Chapman correction factor
!       cmu      :  Abscissae for gauss quadrature over angle cosine
!       ncut     :  Number of computational layer where absorption
!                     optical depth exceeds -abscut-
!       oprim    :  Delta-m scaled single scattering albedo
!       pkag,c   :  Planck function in each layer
!       spher    :  spher = true => spherical geometry invoked
!       taucpr   :  Cumulative optical depth (delta-m-scaled)
!       (remainder are 'twostr' input variables)
!
!   O u t p u t     v a r i a b l e s:
!
!       kk       :  Eigenvalues 
!       rr       :  Eigenvectors at polar quadrature angles 
!       xp_0,    :  Thermal source coefficients x-sub-zero and
!        xp_1         x-sub-one in Eq.  KST(22)
!       yb_0d,u, :  Thermal source vectors, Eqs. KST(24-25)
!        yb_1d,u     
!       yp_0d,u, :  Beam source vectors, Eq. KST(24-25)
!        yp_1d,u     
!       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
!       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
!
!


      implicit none
      LOGICAL planck, spher
      REAL*8 ch(*), chtau(*), ggprim(*), kk(*), large, oprim(*), &
           pkag(0:*), pkagc(*), rr(*), taucpr(0:*), xb_0d(*), &
           xb_0u(*), xb_1d(*), xb_1u(*), xp_0(*), xp_1(*), yb_0d(*), &
           yb_0u(*), yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*), &
           yp_1u(*), zb_a(*), zp_a(*)
!
      real*8 cmu,fbeam,radius,umu0,pi,fact1,fact2,fact3,q0,q0a
      real*8 q2,q2a,q_1,q_2,qq,sgn,small,z0m,z0p,q1,q1a
      real*8 biggest,smallest,arg,beta,big,deltat,denomb,denomp
      integer ncut,nlyr,LC,ls
! The calculation of the particular solutions require some care, small
! and big have been set so that no problems should occurr on 32-bits
! machine running single precision
!
      !IF ((ii == ni).AND.(jj == nj)) print*,'biggest,smallest',biggest,smallest
      !big   = sqrt(biggest) / 1.E+10
      biggest = 3.4E+38
      big   = sqrt(3.4E+38) / 1.E+10
      !small = 1.E+30*smallest
      smallest = 1.2E-38
      small = 1.E+30*1.2E-38
      !IF ((ii == ni).AND.(jj == nj)) print*,'big,small',big,small

!
! ===================  Begin loop on computational layers  =============
!
      DO 100  LC = 1, NCUT
         ls = 2*(lc-1) + 1
!
! Calculate eigenvalues -kk- and eigenvector -rr-, Eqs. KST(20-21) 
!
         beta   = 0.5 * ( 1.-3.*ggprim(lc)*cmu*cmu )
         fact1  = 1. - oprim(lc)
         fact2  = 1. - oprim(lc) + 2.*oprim(lc)*beta 
         kk(lc) = (1./cmu) * sqrt( fact1 * fact2 )
         rr(lc) = ( sqrt(fact2)-sqrt(fact1) ) / &
                                          ( sqrt(fact2)+sqrt(fact1) )
!
         IF ( fbeam.GT.0.0 ) THEN
!
! Set coefficients in KST(22) for beam source
!
           q_1 = (fbeam/(4.*2. * ASIN(1.d0)))*oprim(lc)*(1.-3.*ggprim(lc)*cmu*umu0)
           q_2 = (fbeam/(4.*2. * ASIN(1.d0)))*oprim(lc)*(1.+3.*ggprim(lc)*cmu*umu0) 
!
           IF ( umu0 .GE. 0.0 ) THEN
              qq = q_2
           ELSE
              qq = q_1
           ENDIF
!
           IF ( spher ) THEN
             q0a = DEXP(-chtau(ls-1) )
             q0 = q0a*qq
             IF ( q0 .LE. small) THEN
                q1a = 0.0
                q2a = 0.0
             ELSE
                q1a = DEXP(-chtau(ls) )
                q2a = DEXP(-chtau(ls+1) )
             ENDIF
           ELSE IF ( .NOT. spher ) THEN
             q0a = DEXP(-taucpr(lc-1)/umu0)
             q0 = q0a*qq
             IF ( q0 .LE. small) THEN
                q1a = 0.0
                q2a = 0.0
             ELSE
                q1a = DEXP(- ( taucpr(lc-1)+taucpr(lc) )/ (2.*umu0) )
                q2a = DEXP(-taucpr(lc)/umu0)
             ENDIF
           ENDIF
           q1 = q1a*qq
           q2 = q2a*qq
!
! Calculate alpha coefficient 
!
           deltat = taucpr(lc) - taucpr(lc-1)
           zb_a(lc) = 1./ch(lc)
           large = LOG(biggest)-20. 
           IF( ABS(zb_a(lc)*taucpr(lc-1)) .GT. large .OR. &
                ABS(zb_a(lc)*taucpr(lc)) .GT. large ) zb_a(lc) = 0.0
!
! Dither alpha if it is close to an eigenvalue
!
           denomb =  fact1 * fact2 - (zb_a(lc)*cmu)*(zb_a(lc)*cmu)
           IF ( denomb .LT. 1.E-03 ) THEN
              zb_a(lc) = 1.02*zb_a(lc)
           ENDIF
           q0 = q0a * q_1
           q2 = q2a * q_1
!
! Set constants in Eq. KST(22)
!
           IF ( deltat .LT. 1.E-07 ) THEN
              xb_1d(lc) = 0.0
           ELSE
              xb_1d(lc) = (1./deltat)*(q2*DEXP(zb_a(lc)*taucpr(lc)) &
                   -q0*DEXP(zb_a(lc)*taucpr(lc-1)))
           ENDIF
           xb_0d(lc) = q0 * DEXP(zb_a(lc)*taucpr(lc-1)) - &
                xb_1d(lc)*taucpr(lc-1)
           q0 = q0a * q_2
           q2 = q2a * q_2
           IF ( deltat .LT. 1.E-07 ) THEN
              xb_1u(lc) = 0.0
           ELSE
              xb_1u(lc) = (1./deltat)*(q2*DEXP(zb_a(lc)*taucpr(lc)) &
                   -q0*DEXP(zb_a(lc)*taucpr(lc-1)))
           ENDIF
           xb_0u(lc) = q0 * DEXP(zb_a(lc)*taucpr(lc-1)) - &
                xb_1u(lc)*taucpr(lc-1)
!     
! Calculate particular solutions for incident beam source in 
! pseudo-spherical geometry, Eqs. KST(24-25)
!
           denomb =  fact1 * fact2 - (zb_a(lc)*cmu)*(zb_a(lc)*cmu)
           yb_1d(lc) = (  oprim(lc)*beta*xb_1d(lc) + &
                (1.-oprim(lc)+oprim(lc)*beta+zb_a(lc)*cmu)*xb_1u(lc))/ &
               denomb
           yb_1u(lc) = (  oprim(lc)*beta*xb_1u(lc) + &
                (1.-oprim(lc)+oprim(lc)*beta-zb_a(lc)*cmu)*xb_1d(lc))/ &
                denomb 
           z0p = xb_0u(lc) - cmu*yb_1d(lc)
           z0m = xb_0d(lc) + cmu*yb_1u(lc)
           yb_0d(lc) = (  oprim(lc)*beta*z0m + &
                (1.-oprim(lc)+oprim(lc)*beta+zb_a(lc)*cmu)*z0p)/ &
                denomb
           yb_0u(lc) = (  oprim(lc)*beta*z0p + &
                (1.-oprim(lc)+oprim(lc)*beta-zb_a(lc)*cmu)*z0m)/ &
                denomb
!     
         ENDIF
!
         IF ( planck  ) THEN
!
! Set coefficients in KST(22) for thermal source
!
! Calculate alpha coefficient 
!
            small = 1.E+20*smallest
            q0 = (1.-oprim(lc)) * pkag(lc-1)
            q1 = (1.-oprim(lc)) * pkagc(lc)
            q2 = (1.-oprim(lc)) * pkag(lc)
            deltat = taucpr(lc) - taucpr(lc-1)
!
! Case 1: source small at bottom layer
!
            !IF ((ii == ni).AND.(jj == nj)) print*,'q0,q1,q2',q0,q1,q2 
            IF ( (q2.LT.(q0*1.E-02) .OR. q2.LE.small ) &
                 .AND. q1.GT.small .AND. q0.GT.small ) THEN
!
! alpha Eq. KS(50)
!
               !IF ((ii == ni).AND.(jj == nj)) print*,'CASE 1'
               zp_a(lc) = (2./deltat) * LOG( q0/q1 )
               IF ( zp_a(lc) .GT. big )    zp_a(lc) = big
               IF ( zp_a(lc)*taucpr(lc-1) .GE. LOG(big) ) THEN
                  xp_0(lc) =  big
               ELSE
                  xp_0(lc) = q0
               ENDIF
               xp_1(lc) = 0.0
!     
! Case 2: Source small at center and bottom of layer
!
            ELSE IF ( (q2.LE.(q1*1.E-02) .OR. q2.LE.small ) .AND. &
                    ((q1.LE.(q0*1.E-02)) .OR. (q1.LE.small)) &
                    .AND. (q0.GT.small) ) THEN
!     
               !IF ((ii == ni).AND.(jj == nj)) print*,'CASE 1'
               zp_a(lc)  =   big / taucpr(ncut)
               xp_0(lc) = q0
               xp_1(lc) = 0.0
!     
!     Case 3:All sources zero
!
            ELSE IF ( q2.LE.small .AND. q1.LE.small &
                    .AND. q0.LE.small) THEN
               zp_a(lc)  = 0.0
               xp_0(lc) = 0.0
               xp_1(lc) = 0.0
               !IF ((ii == ni).AND.(jj == nj)) print*,'CASE 3'
!     
!     Case 4: Sources same at center, bottom and top of layer
!     or layer optically very thin
!
            ELSE IF ( (ABS((q2-q0)/q2).LT.1.E-04) .AND. &
                    (ABS((q2-q1)/q2).LT.1.E-04) &
                    .OR. deltat.LT. 1.E-04           ) THEN
!     
               zp_a(lc)  = 0.0
               xp_0(lc) = q0
               xp_1(lc) = 0.0
               !IF ((ii == ni).AND.(jj == nj)) print*,'CASE 4'
!     **  Case 5: Normal case
            ELSE
               arg = (q1/q2)**2. - q0/q2
               IF ( arg .LT. 0.0 ) arg = 0.0
               !IF ((ii == ni).AND.(jj == nj)) print*,'CASE 5'
!     
! alpha Eq. (44). For source that has its maximum at the top of the
! layer, use negative solution
!
               sgn = 1.0
               IF ( pkag(lc-1) .GT. pkag(lc) ) sgn = -1.
               fact3 = LOG(q1/q2 + sgn*SQRT(arg) )
               IF ( ABS(fact3) .LE. 0.005 ) THEN ! Be careful with log of
                  q1 = 0.99 * q1 ! numbers close to one
                  fact3 = LOG(q1/q2 + sgn*SQRT(arg) )
               ENDIF 
               zp_a(lc) = (2./deltat) * fact3
               IF(ABS(zp_a(lc)*taucpr(lc)) .GT. &
                    (LOG(biggest)-LOG(q0*100.) ) )    zp_a(lc) = 0.0 
!     
! Dither alpha if it is close to an eigenvalue
!
               denomp =  fact1 * fact2 - (zp_a(lc)*cmu)*(zp_a(lc)*cmu)
               IF ( denomp .LT. 1.E-03 ) THEN
                  zp_a(lc) = 1.01*zp_a(lc)
               ENDIF
!
! Set constants in Eqs. KST(22)
!
               IF ( deltat .LT. 1.E-07 ) THEN
                  xp_1(lc) = 0.0
               ELSE
                  xp_1(lc) = (1./deltat)*(q2*DEXP(zp_a(lc)*taucpr(lc)) &
                       -q0*DEXP(zp_a(lc)*taucpr(lc-1)))
               ENDIF
               xp_0(lc) = q0 * DEXP(zp_a(lc)*taucpr(lc-1)) - &
                    xp_1(lc)*taucpr(lc-1)
            ENDIF
!     
! Calculate particular solutions Eqs. KST(24-25) for internal thermal
! source
!
            denomp =  fact1 * fact2 - (zp_a(lc)*cmu)*(zp_a(lc)*cmu)
            yp_1d(lc) = (  oprim(lc)*beta*xp_1(lc) + &
                 (1.-oprim(lc)+oprim(lc)*beta+zp_a(lc)*cmu)*xp_1(lc))/ &
                 denomp
            yp_1u(lc) = (  oprim(lc)*beta*xp_1(lc) + &
                 (1.-oprim(lc)+oprim(lc)*beta-zp_a(lc)*cmu)*xp_1(lc))/ &
                 denomp
            z0p = xp_0(lc) - cmu*yp_1d(lc)
            z0m = xp_0(lc) + cmu*yp_1u(lc)
            yp_0d(lc) = (  oprim(lc)*beta*z0m + &
                 (1.-oprim(lc)+oprim(lc)*beta+zp_a(lc)*cmu)*z0p)/ &
                 denomp
            yp_0u(lc) = (  oprim(lc)*beta*z0p + &
                 (1.-oprim(lc)+oprim(lc)*beta-zp_a(lc)*cmu)*z0m)/ &
                 denomp
!     
         END IF
!     
 100  CONTINUE
!
! ===================  End loop on computational layers  ===============
!
      RETURN
      END SUBROUTINE
!
      SUBROUTINE tprtinp( albedo, btemp, deltam, dtauc, fbeam, fisot, &
           flyr, gg, lyrcut, nlyr, planck, nstr, ntau, oprim, &
           spher, ssalb, tauc, taucpr, temis, temper, ttemp, umu0, &
           utau, wvnmlo, wvnmhi )
!
! Print values of input variables
!
      implicit none
      integer nlyr,nstr,ntau,lc,lu
      real*8 yessct
      LOGICAL  deltam, lyrcut, planck, spher
      REAL*8     dtauc(*), flyr(*), gg(*), oprim(*), ssalb(*), &
               tauc( 0:* ), taucpr( 0:* ), temper( 0:* ), utau(*)
!
      real*8 albedo,fbeam,btemp,fisot,temis,ttemp,umu0
      real*8 wvnmlo, wvnmhi

      WRITE( *,1010 )  nstr, nlyr
      WRITE( *,1030 )  ntau, (utau(lu), lu = 1, ntau)
      IF ( spher ) WRITE(*,1090)
      IF ( .NOT. planck  )  WRITE( *,1100 )
      WRITE( *,1060 ) fbeam, umu0,  fisot
      WRITE( *,1080 ) albedo
      IF ( planck )  WRITE( *,1110 ) wvnmlo, wvnmhi, btemp, &
                                          ttemp, temis
      IF ( deltam )       WRITE( *,1120 )
      IF ( .NOT.deltam )  WRITE( *,1130 )
      IF ( lyrcut )       WRITE( *,1170 )
      IF ( planck )       WRITE ( *,1190 )
      IF ( .NOT. planck ) WRITE ( *,1191 )
      yessct = 0.0
      DO 10 lc = 1, nlyr
         yessct = yessct + ssalb(lc)
         IF( planck ) &
             WRITE( *,1200 )  lc, dtauc(lc), tauc(lc), ssalb(lc), &
              flyr(lc), taucpr(lc)-taucpr(lc-1), taucpr(lc), &
              oprim(lc), gg(lc), temper(lc-1)
         IF( .NOT.  planck ) &
             WRITE( *,1200 )  lc, dtauc(lc), tauc(lc), ssalb(lc), &
              flyr(lc), taucpr(lc)-taucpr(lc-1), taucpr(lc), &
              oprim(lc), gg(lc)
 10   CONTINUE
      IF( planck )  WRITE( *,1210 ) temper(nlyr)
!
      RETURN
!
1010  FORMAT ( /, ' No. streams =', I4, &
        '     No. computational layers =', I4 )
1030  FORMAT( I4,' User optical depths :',10F10.4, /, (26X,10F10.4) )
1060  FORMAT( '    Incident beam with intensity =', 1P,E11.3, ' and', & 
        ' polar angle cosine = ', 0P,F8.5, & 
        /,'    plus isotropic incident intensity =', 1P,E11.3 )
1080  FORMAT( '    Bottom albedo (lambertian) =', 0P,F8.4 )
1090  FORMAT( ' Pseudo-spherical geometry invoked' )
1100  FORMAT( ' No thermal emission' )
1110  FORMAT( '    Thermal emission in wavenumber interval :', 2F14.4,/, & 
          '    bottom temperature =', F10.2, '     top temperature =', & 
          F10.2,'    top emissivity =', F8.4 )
1120  FORMAT( ' Uses delta-m method' )
1130  FORMAT( ' Does not use delta-m method' )
1150  FORMAT( ' Calculate fluxes and intensities' )
1170  FORMAT( ' Sets radiation = 0 below absorption optical depth 10' )
1190  FORMAT( /, 37X, '<------------- delta-m --------------->', /, & 
       '                   total    single                           ', & 
                      'total    single', /, & 
       '       optical   optical   scatter   truncated   ', & 
          'optical   optical   scatter    asymm', /, & 
       '         depth     depth    albedo    fraction     ', & 
            'depth     depth    albedo   factor   temperature' )
1191  FORMAT( /, 37X, '<------------- delta-m --------------->', /, & 
       '                   total    single                           ', & 
                      'total    single', /, & 
       '       optical   optical   scatter   truncated   ', & 
          'optical   optical   scatter    asymm', /, & 
       '         depth     depth    albedo    fraction     ', & 
            'depth     depth    albedo   factor' )
1200  FORMAT( I4, 2F10.4, F10.5, F12.5, 2F10.4, F10.5, F9.4,F14.3 )
1210  FORMAT( 85X, F14.3 )
1300  FORMAT( I6, 10F11.6, /, (6X,10F11.6) )
!
      END SUBROUTINE
!
      SUBROUTINE settwo( albedo, biggest, bplank, btemp, ch, chtau, & 
            cmu, deltam,  & 
            dtauc, dtaucpr, expbea, fbeam, flyr, gg, ggprim, & 
            layru, lyrcut, mxcly, ncut, nlyr, nn, nstr, ntau, oprim, & 
            pkag, pkagc, planck, radius, smallest, spher, ssalb, tauc,  & 
            taucpr, temis, temper, tplank, ttemp, umu0, usrtau, & 
            utau, utaupr, wvnmlo, wvnmhi, zd)
!
! Perform miscellaneous setting-up operations
!
! Routines called:  errmsg, zeroit
!
! Input :  All are 'twostr' input variables (see doc file)
!
! Output:  ntau,utau   If usrtau = false
!          bplank      Intensity emitted from bottom boundary
!          ch          The Chapman factor
!          cmu         Computational polar angle
!          expbea      Transmission of direct beam
!          flyr        Truncated fraction in delta-m method
!          layru       Computational layer in which -utau- falls
!          lyrcut      Flag as to whether radiation will be zeroed
!                      below layer -ncut-
!          ncut        Computational layer where absorption
!                      optical depth first exceeds -abscut-
!          nn          nstr / 2  =  1
!          nstr        No. of streams (=2)
!          oprim       Delta-m-scaled single-scatter albedo
!          pkag,c      Planck function in each layer
!          taucpr      Delta-m-scaled optical depth
!          tplank      Intensity emitted from top boundary
!          utaupr      Delta-m-scaled version of -utau-
!
! Internal variables
!          abscut      Absorption optical depth, medium is cut off below
!                      this depth
!          tempc       Temperature at center of layer, assumed
!                      to be average of layer boundary temperatures
!
      implicit none
      real*8 biggest,bplank,smallest,tplank,abscut,abstau,pi,f
      integer mxcly,ncut,nlyr,nn,nstr,ntau,lc,lev,ls,lu
      LOGICAL  deltam, lyrcut, planck, spher, usrtau, first
      INTEGER  layru(*)
      REAL*8 ch(*), chtau(0:*), cmu, dtauc(*), dtaucpr(*), &
           expbea(0:*),flyr(*),gg(*),ggprim(*),oprim(*),pkag(0:*), &
           pkagc(*), ssalb(*), tauc(0:*), taucpr(0:*), &
           temper(0:*), utau(*), utaupr(*), zd(0:*)
!
      real*8 umu0,zenang,radius,taup,wvnmlo, wvnmhi, ttemp,btemp
      real*8 tempc,albedo,fbeam,temis
      DATA  abscut / 10. /, first / .TRUE. /
      !real*8 chapmn,tplkavg,d1mach
!
      IF ( first ) THEN
         first    = .FALSE.
         !smallest = d1mach(1)
         smallest = 1.2E-38
         !biggest  = d1mach(2)
         biggest  = 3.4E+38
         pi       = 2. * ASIN( 1.0 )
         nstr     = 2
         nn       = nstr / 2
      ENDIF
!
      IF ( .NOT.usrtau ) THEN
!
! Set output levels at computational layer boundaries
!
         ntau = nlyr + 1
         DO 30  lc = 0, ntau-1
            utau(lc+1) = tauc(lc)
 30      CONTINUE
      END IF
!
! Apply delta-m scaling and move description of computational
! layers to local variables
!
      expbea( 0 ) = 1.0
      !print*,'umu0,pi',umu0,pi
      zenang      = DACOS(umu0) * 180. / (2. * ASIN(1.d0))
      IF( spher .AND. umu0 .LT. 0.0 ) expbea(0) = &
                DEXP(-chapmn(1,0.d0,tauc,nlyr, zd,dtauc,zenang,radius) )
      CALL  zeroit( taucpr(0), mxcly+1 )
      CALL  zeroit( expbea(1), mxcly )
      CALL  zeroit( flyr, mxcly )
      CALL  zeroit( oprim, mxcly )
      abstau = 0.0
      DO  60  lc = 1, nlyr
         IF ( abstau.LT.abscut )  ncut = lc
         abstau = abstau + ( 1. - ssalb(lc) ) * dtauc(lc)
!
         IF ( .NOT.deltam )  THEN
            oprim(lc)  = ssalb(lc)
            taucpr(lc) = tauc(lc)
            f          = 0.0
            ggprim(lc) = gg(lc)
            dtaucpr(lc)= dtauc(lc)
         ELSE
!
! Do delta-m transformation Eqs. WW(20a,20b,14)
!
            f = gg(lc) * gg(lc)
            taucpr(lc) = taucpr(lc-1) + ( 1. - f*ssalb(lc) ) * dtauc(lc)
            oprim(lc)  = ssalb(lc) * ( 1. - f ) / ( 1. - f * ssalb(lc) )
            ggprim(lc) =  (gg(lc)-f) / (1.-f)
            dtaucpr(lc)= taucpr(lc) - taucpr(lc-1)
         ENDIF
!
         flyr(lc)   = f
!
 60   CONTINUE
!
! If no thermal emission, cut off medium below absorption optical
! depth = abscut ( note that delta-m transformation leaves absorption
! optical depth invariant ).  Not worth the trouble for one-layer
! problems, though.
!
      lyrcut = .FALSE.
      IF ( abstau.GE.abscut .AND. .NOT. planck &
                                 .AND. nlyr.GT.1 )  lyrcut =.TRUE.
      IF ( .NOT.lyrcut )  ncut = nlyr
!
! Calculate chapman function is spherical geometry, set expbea and ch
! for beam source.
!
      IF ( fbeam.GT.0.0 )  THEN
         chtau(0) = 0.0
         DO lc = 1, ncut
            expbea(lc) = 0.0
            IF ( spher ) THEN
               ls = 2*(lc-1) + 1
               taup   = taucpr(lc-1) + dtaucpr(lc)/2.0
               chtau(ls) = chapmn(lc,taup,taucpr,nlyr, &
                    zd,dtaucpr,zenang,radius)
               chtau(ls+1) = chapmn(lc,taucpr(lc),taucpr,nlyr, &
                    zd,dtaucpr,zenang,radius)
               ch(lc) = taup/chtau(ls)
               expbea(lc) = DEXP(-chtau(ls+1) )
            ELSE IF ( .NOT. spher ) THEN
               ch(lc)     = umu0
               expbea(lC) = DEXP( - taucpr(lc) / umu0 )
            ENDIF
         ENDDO
      ENDIF
!
! Set arrays defining location of user output levels within 
! delta-m-scaled computational mesh
! 
      DO 90  lu = 1, ntau
         DO 70 lc = 1, nlyr
            IF ( utau(lu).GE.tauc(lc-1) .AND. utau(lu).LE.tauc(lc) ) &
                 GO TO 80
 70      CONTINUE
         lc = nlyr
!
 80      utaupr(lu) = utau(lu)
         IF(deltam) utaupr(lu) = taucpr(lc-1)+(1.-ssalb(lc)*flyr(lc)) &
                                             * (utau(lu) - tauc(lc-1))
         layru(lu) = lc
 90   CONTINUE
!
! Set computational polar angle cosine for double gaussian 
! quadrature; cmu = 0.5, or  single gaussian quadrature; cmu =
! 1./sqrt(3)
! See KST for discussion of which is better for your specific
! application
!
      IF ( planck .AND. fbeam .EQ. 0.0 ) THEN
         cmu =  0.5
      ELSE
         cmu = 1./SQRT(3.0)
      ENDIF
!
! Calculate planck functions
!
      IF ( .NOT. planck )  THEN
         bplank = 0.0
         tplank = 0.0
         CALL  zeroit( pkag, mxcly+1 )
         CALL  zeroit( pkagc, mxcly )
      ELSE
         tplank = temis * tplkavg( wvnmlo, wvnmhi, ttemp )
         bplank =         tplkavg( wvnmlo, wvnmhi, btemp )
         DO 180  lev = 0, nlyr            
            pkag( lev ) = tplkavg( wvnmlo, wvnmhi, temper(lev) )
180         CONTINUE
         DO 190 lc = 1, nlyr
            tempc = 0.5 * ( temper(lc-1) + temper(lc) )
            pkagc( lc ) = tplkavg( wvnmlo, wvnmhi, tempc )
 190     CONTINUE
      END IF
      RETURN
      END SUBROUTINE
!
      SUBROUTINE solbcd( albedo, b, bplank, cband, cmu, diag, expbea, &
           fbeam, fisot, ipvt, kk, ll, lyrcut, mi, mi9m2, mxcmu, &
           ncol, ncut, nlyr, nn, nnlyri, nstr, pi, rr, subd, &
           superd, sol,tplank, taucpr, umu0, yb_0d, yb_0u, yb_1d, &
           yb_1u, yp_0d, yp_0u, yp_1d, yp_1u, zb_a, zp_a )
!
! Construct right-hand side vector -b- for general boundary conditions
! and solve system of equations obtained from the boundary conditions
! and the continuity-of-intensity-at-layer-interface equations.
!
! Routines called: sgbfa, sgbsl, zeroit
!
! I n p u t      v a r i a b l e s:
!
!       bplank   :  Bottom boundary thermal emission
!       cband    :  Left-hand side matrix of linear system Eqs.
!       KST(38-41)
!                   in banded form required by linpack solution routines
!       cmu      :  Abscissae for gauss quadrature over angle cosine
!       expbea   :  Transmission of incident beam, EXP(-taucpr/ch)
!       lyrcut   :  Logical flag for truncation of comput. layer
!       ncol     :  Counts of columns in -cband-
!       nn       :  Order of double-gauss quadrature (nstr/2)
!       ncut     :  Total number of computational layers considered
!       nstr     :  No. of streams (=2)
!       tplank   :  Top boundary thermal emission
!       taucpr   :  Cumulative optical depth (delta-m-scaled)
!       yb_0d,u, :  Thermal source vectors, Eq. KST(24-25)
!        yb_1d,u     
!       yp_0d,u, :  Beam source vectors, Eq. KST(24-25)
!        yp_1d,u     
!       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
!       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
!       (Remainder are 'twostr' input variables)
!
! O u t p u t     v a r i a b l e s:
!
!       b        :  Right-hand side vector of Eqs. KST(38-41) going into
!                   *sgbsl*; returns as solution vector of Eqs.
!                   KST(38-41)
!                   constants of integration 
!      ll        :  Permanent storage for -b-, but re-ordered
!
! I n t e r n a l    v a r i a b l e s:
!
!       it       :  Pointer for position in  -b-
!       ncd      :  Number of diagonals below or above main diagonal
!       rcond    :  Indicator of singularity for -cband-
!       z        :  Scratch array required by *sgbco*
!+---------------------------------------------------------------------+
!
      implicit none
      real*8 bplank,pi,tplank,refflx,rp_m,rp_p,rpp1_p,sum,wk,wk0,wk1
      real*8 rpp1_m
      integer mxcmu,mi9m2,mi,ncol,ncut,nlyr,nn,nnlyri,nstr,lc,nloop
      integer info,irow,job,nrow
      INTEGER  ipvt(*)
      LOGICAL  lyrcut
      REAL*8 b(*), cband( mi9m2,nnlyri ), cmu, diag(*), expbea(0:*), &
           kk(*), ll( mxcmu,* ), rr(*), subd(*), superd(*), sol(*), &
           taucpr( 0:* ), yb_0d(*), yb_0u(*), yb_1d(*), yb_1u(*), &
           yp_0d(*), yp_0u(*), yp_1d(*), yp_1u(*), zb_a(*), zp_a(*)
!
      real*8 albedo,fbeam,fisot,umu0
! First top row, top boundary condition
!
      irow        = 1
      lc          = 1
!     subd(irow)  = undefined
      diag(irow)  = rr(lc) * DEXP(-kk(lc) * taucpr(lc))
      superd(irow)= 1.0
!
! next from layer no. 2 to nlyr -1
!
      nloop = ncut - 1 
      DO lc = 1, nloop
         irow         = irow + 1
         wk0          = DEXP(-kk(lc) * (taucpr(lc) - taucpr(lc-1)))
         wk1          = DEXP(-kk(lc+1) * (taucpr(lc+1) - taucpr(lc)))
         subd(irow)   = 1.0 - rr(lc) * rr(lc+1)
         diag(irow)   = ( rr(lc) -  rr(lc+1 ) ) * wk0
         superd(irow) = - ( 1. - rr(lc+1) * rr(lc+1 ) ) * wk1
         irow         = irow + 1
         subd(irow)   = ( 1.0 - rr(lc) * rr(lc) ) * wk0 
         diag(irow)   = ( rr(lc) -  rr(lc+1 ) ) * wk1
         superd(irow) = - ( 1. - rr(lc+1) * rr(lc ) ) 
      ENDDO
!
! bottom layer
!
      irow         = irow + 1
      lc           = ncut
!     superd(irow) = undefined
      wk           = DEXP( -kk(lc) * (taucpr(lc) - taucpr(lc-1)) )
      IF ( lyrcut ) THEN
         subd(irow) = 1.0
         diag(irow) = rr(lc) * wk
      ELSE
         subd(irow) = 1.0 - 2.0*albedo*cmu*rr(lc)
         diag(irow) = ( rr(lc) - 2.0*albedo*cmu ) * wk         
      ENDIF
!
      CALL  zeroit( b, nnlyri )
!
! Construct -b-,  for parallel beam + bottom reflection +
! thermal emission at top and/or bottom
!
! Top boundary, right-hand-side of Eq. KST(28)
!
      lc      = 1
      irow    = 1
      b(irow) = - yb_0d(lc) - yp_0d(lc) +fisot +tplank
!
! Continuity condition for layer interfaces,
! right-hand-side of Eq. KST(29)
!
      DO   lc = 1, nloop
         rpp1_m = DEXP(-zb_a(lc+1)*taucpr(lc))* & 
               (yb_0d(lc+1)+yb_1d(lc+1)*taucpr(lc))  & 
               + DEXP(-zp_a(lc+1)*taucpr(lc))* & 
               (yp_0d(lc+1)+yp_1d(lc+1)*taucpr(lc))
         rp_m   =  DEXP(-zb_a(lc)*taucpr(lc))* & 
               (yb_0d(lc)+yb_1d(lc)*taucpr(lc)) & 
               +  DEXP(-zp_a(lc)*taucpr(lc))* & 
               (yp_0d(lc)+yp_1d(lc)*taucpr(lc))
         rpp1_p = DEXP(-zb_a(lc+1)*taucpr(lc))* & 
               (yb_0u(lc+1)+yb_1u(lc+1)*taucpr(lc)) & 
               +  DEXP(-zp_a(lc+1)*taucpr(lc))* & 
               (yp_0u(lc+1)+yp_1u(lc+1)*taucpr(lc))
         rp_p   = DEXP(-zb_a(lc)*taucpr(lc))* & 
               (yb_0u(lc)+yb_1u(lc)*taucpr(lc)) & 
               +  DEXP(-zp_a(lc)*taucpr(lc))* & 
               (yp_0u(lc)+yp_1u(lc)*taucpr(lc))
         irow    = irow + 1
         b(irow) = rpp1_p - rp_p - rr(lc+1) * ( rpp1_m - rp_m )
         irow    = irow + 1
         b(irow) = rpp1_m - rp_m - rr(lc) * ( rpp1_p - rp_p )
      ENDDO
!
! Bottom boundary
!
      irow = irow + 1
      lc   = ncut
      IF ( lyrcut ) THEN
!
! Right-hand-side of Eq. KST(30)
!
         b(irow) = - DEXP(-zb_a(ncut)*taucpr(ncut))* & 
               (yb_0u(ncut)+yb_1u(ncut)*taucpr(ncut)) & 
               - DEXP(-zp_a(ncut)*taucpr(ncut))* & 
               (yp_0u(ncut)+yp_1u(ncut)*taucpr(ncut))
      ELSE
         sum = cmu * albedo*(DEXP(-zb_a(ncut)*taucpr(ncut))* & 
               (yb_0d(ncut)+yb_1d(ncut)*taucpr(ncut)) & 
               +  DEXP(-zp_a(ncut)*taucpr(ncut))* & 
               (yp_0d(ncut)+yp_1d(ncut)*taucpr(ncut)))
         refflx = 1.
         IF ( umu0 .LE. 0.0 ) refflx = 0.0
         b(irow) = 2.*sum +  & 
               ( albedo * umu0*fbeam/2.*ASIN(1.d0)*refflx ) * expbea(ncut)  & 
               + (1.-albedo) * bplank &
               !+  bplank       &  ! Modify by Xianyu Tan to
!                                         account for a self-luminous
!                                         boundary & 
               -  DEXP(-zb_a(ncut)*taucpr(ncut))* & 
               (yb_0u(ncut)+yb_1u(ncut)*taucpr(ncut)) & 
               -  DEXP(-zp_a(ncut)*taucpr(ncut))* & 
               (yp_0u(ncut)+yp_1u(ncut)*taucpr(ncut))

!         write(*,*) sum, ( albedo * umu0*fbeam/2.*ASIN(1.d0)*refflx ) * expbea(ncut), bplank  
      END IF
!
! solve for constants of integration by inverting matrix KST(38-41)
!
      nrow = irow
!
          CALL zeroit( cband, mi9m2*nnlyri)
          DO irow = 1, nrow
             cband(1,irow)   = 0.0
             cband(3,irow)   = diag(irow)
          ENDDO
          DO irow = 1, nrow-1
             cband(2,irow+1) = superd(irow)
          ENDDO
          DO irow = 2, nrow
             cband(4,irow-1) = subd(irow)
          ENDDO
!
          CALL sgbfa(cband, mi9m2, nrow, 1, 1, ipvt, info )
          job = 0
          CALL sgbsl(cband, mi9m2, nrow, 1, 1, ipvt, b, job )
!
! unpack
!
          irow = 1
          DO lc = 1, ncut
             ll(1,lc) = b(irow)      ! downward direction
             irow     = irow + 1
             ll(2,lc) = b(irow)      ! upward direction
             irow     = irow + 1
          ENDDO
!
      RETURN
      END SUBROUTINE
!
      SUBROUTINE tzeroal( ch, kk, ll, mxcmu, mxcly, &
           nnlyri, rr, xb_0d, xb_0u, xb_1d, xb_1u, &
           xp_0, xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, &
           yp_0u, yp_1d, yp_1u, zb_a, zp_a )
!
! Zero arrays
!
      implicit none
      integer mxcmu,mxcly,nnlyri
      REAL*8  ch(*), kk(*), ll(mxcmu,*), rr(*), &
            xb_0d(*), xb_0u(*), xb_1d(*), &
            xb_1u(*), xp_0(*), xp_1(*), yb_0d(*), yb_0u(*), &
            yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*), & 
            yp_1u(*), zb_a(*), zp_a(*)
!
      CALL  zeroit( ch    , mxcly )
      CALL  zeroit( kk,     mxcly )
      CALL  zeroit( ll,     mxcmu*mxcly )
      CALL  zeroit( rr,     mxcly )
      CALL  zeroit( xb_0d, mxcly )
      CALL  zeroit( xb_0u, mxcly )
      CALL  zeroit( xb_1d, mxcly )
      CALL  zeroit( xb_1u, mxcly )
      CALL  zeroit( xp_0, mxcly )
      CALL  zeroit( xp_1, mxcly )
      CALL  zeroit( yb_0d, mxcly )
      CALL  zeroit( yb_0u, mxcly )
      CALL  zeroit( yb_1d, mxcly )
      CALL  zeroit( yb_1u, mxcly )
      CALL  zeroit( yp_0d, mxcly )
      CALL  zeroit( yp_0u, mxcly )
      CALL  zeroit( yp_1d, mxcly )
      CALL  zeroit( yp_1u, mxcly )
      CALL  zeroit( zb_a, mxcly )
      CALL  zeroit( zp_a, mxcly )
!
      RETURN
      END SUBROUTINE
!
      SUBROUTINE zeroit( A, LENGTH )
!
!         ZEROS A REAL ARRAY -A- HAVING -LENGTH- ELEMENTS
!
      implicit none
      integer length,L
      REAL*8  A(*)
!
      DO 10  L = 1, LENGTH
         A( L ) = 0.0
10    CONTINUE
!
      RETURN
      END SUBROUTINE

      real*8 function F(X)

      implicit none
      real*8 X
      
      F = X**3 / ( DEXP(X) - 1 )  
      return
      end function


      REAL*8   FUNCTION  TPLKAVG ( WNUMLO, WNUMHI, T )
!
!        COMPUTES PLANCK FUNCTION INTEGRATED BETWEEN TWO WAVENUMBERS,
!        except if wnmulo .EQ. wnmuhi, then the Planck function at 
!        wnumlo is returned
!
!  NOTE ** CHANGE 'R1MACH' TO 'D1MACH' TO RUN IN DOUBLE PRECISION
!
!  I N P U T :  WNUMLO : LOWER WAVENUMBER ( INV CM ) OF SPECTRAL
!                           INTERVAL
!               WNUMHI : UPPER WAVENUMBER
!               T      : TEMPERATURE (K)
!
!  O U T P U T :  PLKAVG : INTEGRATED PLANCK FUNCTION ( WATTS/SQ M )
!                           = INTEGRAL (WNUMLO TO WNUMHI) OF
!                              2H C**2  NU**3 / ( EXP(HC NU/KT) - 1)
!                              (WHERE H=PLANCKS CONSTANT, C=SPEED OF
!                              LIGHT, NU=WAVENUMBER, T=TEMPERATURE,
!                              AND K = BOLTZMANN CONSTANT)
!
!  REFERENCE : SPECIFICATIONS OF THE PHYSICAL WORLD: NEW VALUE
!                 OF THE FUNDAMENTAL CONSTANTS, DIMENSIONS/N.B.S.,
!                 JAN. 1974
!
!  METHOD :  FOR  -WNUMLO-  CLOSE TO  -WNUMHI-, A SIMPSON-RULE
!            QUADRATURE IS DONE TO AVOID ILL-CONDITIONING; OTHERWISE
!
!            (1)  FOR WAVENUMBER (WNUMLO OR WNUMHI) SMALL,
!                 INTEGRAL(0 TO WNUM) IS CALCULATED BY EXPANDING
!                 THE INTEGRAND IN A POWER SERIES AND INTEGRATING
!                 TERM BY TERM;
!
!            (2)  OTHERWISE, INTEGRAL(WNUMLO/HI TO INFINITY) IS
!                 CALCULATED BY EXPANDING THE DENOMINATOR OF THE
!                 INTEGRAND IN POWERS OF THE EXPONENTIAL AND
!                 INTEGRATING TERM BY TERM.
!
!  ACCURACY :  AT LEAST 6 SIGNIFICANT DIGITS, ASSUMING THE
!              PHYSICAL CONSTANTS ARE INFINITELY ACCURATE
!
!  ERRORS WHICH ARE NOT TRAPPED:
!
!      * POWER OR EXPONENTIAL SERIES MAY UNDERFLOW, GIVING NO
!        SIGNIFICANT DIGITS.  THIS MAY OR MAY NOT BE OF CONCERN,
!        DEPENDING ON THE APPLICATION.
!
!      * SIMPSON-RULE SPECIAL CASE IS SKIPPED WHEN DENOMINATOR OF
!        INTEGRAND WILL CAUSE OVERFLOW.  IN THAT CASE THE NORMAL
!        PROCEDURE IS USED, WHICH MAY BE INACCURATE IF THE
!        WAVENUMBER LIMITS (WNUMLO, WNUMHI) ARE CLOSE TOGETHER.
! ----------------------------------------------------------------------
!                                   *** ARGUMENTS
      implicit none
      REAL*8     T, WNUMLO, WNUMHI
!                                   *** LOCAL VARIABLES
!
!        A1,2,... :  POWER SERIES COEFFICIENTS
!        C2       :  H * C / K, IN UNITS CM*K (H = PLANCKS CONSTANT,
!                      C = SPEED OF LIGHT, K = BOLTZMANN CONSTANT)
!        D(I)     :  EXPONENTIAL SERIES EXPANSION OF INTEGRAL OF
!                       PLANCK FUNCTION FROM WNUMLO (I=1) OR WNUMHI
!                       (I=2) TO INFINITY
!        EPSIL    :  SMALLEST NUMBER SUCH THAT 1+EPSIL .GT. 1 ON
!                       COMPUTER
!        EX       :  EXP( - V(I) )
!        EXM      :  EX**M
!        MMAX     :  NO. OF TERMS TO TAKE IN EXPONENTIAL SERIES
!        MV       :  MULTIPLES OF 'V(I)'
!        P(I)     :  POWER SERIES EXPANSION OF INTEGRAL OF
!                       PLANCK FUNCTION FROM ZERO TO WNUMLO (I=1) OR
!                       WNUMHI (I=2)
!        PI       :  3.14159...
!        SIGMA    :  STEFAN-BOLTZMANN CONSTANT (W/M**2/K**4)
!        SIGDPI   :  SIGMA / PI
!        SMALLV   :  NUMBER OF TIMES THE POWER SERIES IS USED (0,1,2)
!        V(I)     :  C2 * (WNUMLO(I=1) OR WNUMHI(I=2)) / TEMPERATURE
!        VCUT     :  POWER-SERIES CUTOFF POINT
!        VCP      :  EXPONENTIAL SERIES CUTOFF POINTS
!        VMAX     :  LARGEST ALLOWABLE ARGUMENT OF 'EXP' FUNCTION
!
      REAL*8 A1,A2,A3,A4,A5,A6
      real*8 EXM,X,arg,PI,HH,DEL,c1,OLDVAL,plkavg,VAL,VAL0
      real*8 wvn, VMAX
      integer I,N,K,M,MMAX
      PARAMETER  (A1 = 1./3.,A2 = -1./8.,A3 = 1./60.,A4 = -1./5040., &
                   A5 = 1./272160., A6 = -1./13305600. )
      INTEGER  SMALLV
      REAL*8     C2, CONC, D(2), EPSIL, EX, MV, P(2), SIGMA, SIGDPI, &
              V(2), VCUT, VCP(7), VSQ
      SAVE     CONC, VMAX, EPSIL, SIGDPI
      DATA     C2 / 1.438786 /,  SIGMA / 5.67032E-8 /, &
               VCUT / 1.5 /, VCP / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      DATA     PI / 0.0 /
!      real*8 D1MACH, F
!      F(X) = X**3 / ( DEXP(X) - 1 )
!
!
      IF ( PI.EQ.0.0 )  THEN
         PI = 2. * ASIN( 1.0 )
         VMAX = LOG( D1MACH(2) )
         EPSIL = D1MACH(4)
         SIGDPI = SIGMA / PI
         CONC = 15. / PI**4
      END IF
!
      IF( T.LT.0.0 .OR. WNUMHI.LT.WNUMLO .OR. WNUMLO.LT.0. ) &
         CALL ERRMSG( 'PLKAVG--TEMPERATURE OR WAVENUMS. WRONG', .TRUE.)
!
      IF ( T.LT.1.E-4 )  THEN
         TPLKAVG = 0.0
         RETURN
      ENDIF
!
!
      IF ( wnumhi .eq. wnumlo ) THEN
         wvn  =  wnumhi
         arg  = DEXP( - C2 * wvn / T )
         plkavg = c1 * (wvn**3.) * arg / ( 1. - arg )
         RETURN
      ENDIF
!     
      V(1) = C2 * WNUMLO / T
      V(2) = C2 * WNUMHI / T
      IF ( V(1).GT.EPSIL .AND. V(2).LT.VMAX .AND. &
           (WNUMHI-WNUMLO)/WNUMHI .LT. 1.E-2 )  THEN
!
!                          ** WAVENUMBERS ARE VERY CLOSE.  GET INTEGRAL
!                          ** BY ITERATING SIMPSON RULE TO CONVERGENCE.
         HH = V(2) - V(1)
         OLDVAL = 0.0
         VAL0 = F( V(1) ) + F( V(2) )
!
         DO  2  N = 1, 10
            DEL = HH / (2*N)
            VAL = VAL0
            DO  1  K = 1, 2*N-1
               VAL = VAL + 2*(1+MOD(K,2)) * F( V(1) + K*DEL )
    1       CONTINUE
            VAL = DEL/3. * VAL
            IF ( ABS( (VAL-OLDVAL)/VAL ) .LE. 1.E-6 )  GO TO 3
            OLDVAL = VAL
    2    CONTINUE
         CALL ERRMSG( 'PLKAVG--SIMPSON RULE DIDNT CONVERGE', .FALSE. )
!
    3    TPLKAVG = SIGDPI * T**4 * CONC * VAL
         RETURN
      END IF
!
      SMALLV = 0
      DO  50  I = 1, 2
!
         IF( V(I).LT.VCUT )  THEN
!                                   ** USE POWER SERIES
            SMALLV = SMALLV + 1
            VSQ = V(I)**2
            P(I) =  CONC * VSQ * V(I) * ( A1 + V(I) * ( A2 + V(I) * &
                      ( A3 + VSQ * ( A4 + VSQ * ( A5 + VSQ*A6 ) ) ) ) )
         ELSE
!                    ** USE EXPONENTIAL SERIES
            MMAX = 0
!                                ** FIND UPPER LIMIT OF SERIES
   20       MMAX = MMAX + 1
               IF ( V(I).LT.VCP( MMAX ) )  GO TO 20
!
            EX = DEXP( - V(I) )
            EXM = 1.0
            D(I) = 0.0
!
            DO  30  M = 1, MMAX
               MV = M * V(I)
               EXM = EX * EXM
               D(I) = D(I) + &
                      EXM * ( 6. + MV*( 6. + MV*( 3. + MV ) ) ) / M**4
   30       CONTINUE
!
            D(I) = CONC * D(I)
         END IF
!
   50 CONTINUE
!
      IF ( SMALLV .EQ. 2 ) THEN
!                                    ** WNUMLO AND WNUMHI BOTH SMALL
         TPLKAVG = P(2) - P(1)
!
      ELSE IF ( SMALLV .EQ. 1 ) THEN
!                                    ** WNUMLO SMALL, WNUMHI LARGE
         TPLKAVG = 1. - P(1) - D(2)
!
      ELSE
!                                    ** WNUMLO AND WNUMHI BOTH LARGE
         TPLKAVG = D(1) - D(2)
!
      END IF
!
      TPLKAVG = SIGDPI * T**4 * TPLKAVG
      IF( TPLKAVG.EQ.0.0 ) &
          CALL ERRMSG( 'PLKAVG--RETURNS ZERO; POSSIBLE UNDERFLOW', &
                       .FALSE. )
!
      RETURN
      END FUNCTION





      SUBROUTINE  SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )
!
!         FACTORS A REAL BAND MATRIX BY ELIMINATION.
!
!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
!
!     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!
!     INPUT:  SAME AS 'SGBCO'
!
!     ON RETURN:
!
!        ABD,IPVT    SAME AS 'SGBCO'
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
!                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.
!
!     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)
!
!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
!                       FROM FORTRAN: MAX0, MIN0
!
      implicit none
      INTEGER  LDA, N, ML, MU, IPVT(*), INFO
      REAL*8     ABD(LDA,*)
!
      REAL*8     T
      INTEGER  I,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
!
!
      M = ML + MU + 1
      INFO = 0
!                        ** ZERO INITIAL FILL-IN COLUMNS
      J0 = MU + 2
      J1 = MIN0(N, M) - 1
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
      JZ = J1
      JU = 0
!
!                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
      NM1 = N - 1
      DO 120 K = 1, NM1
         KP1 = K + 1
!                                  ** ZERO NEXT FILL-IN COLUMN
         JZ = JZ + 1
         IF (JZ .LE. N) THEN
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0E0
   40       CONTINUE
         ENDIF
!                                  ** FIND L = PIVOT INDEX
         LM = MIN0(ML, N-K)
         L = ISAMAX(LM+1, ABD(M,K), 1) + M - 1
         IPVT(K) = L + K - M
!
         IF (ABD(L,K) .EQ. 0.0E0) THEN
!                                      ** ZERO PIVOT IMPLIES THIS COLUMN
!                                      ** ALREADY TRIANGULARIZED
            INFO = K
         ELSE
!                                ** INTERCHANGE IF NECESSARY
            IF (L .NE. M) THEN
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
            ENDIF
!                                   ** COMPUTE MULTIPLIERS
            T = -1.0E0 / ABD(M,K)
            CALL SSCAL(LM, T, ABD(M+1,K), 1)
!
!                               ** ROW ELIMINATION WITH COLUMN INDEXING
!
            JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
            MM = M
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .NE. MM) THEN
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
               ENDIF
               CALL SAXPY(LM, T, ABD(M+1,K), 1, ABD(MM+1,J), 1)
   80       CONTINUE
!
         ENDIF
!
  120 CONTINUE
!
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0E0) INFO = N
      RETURN
      END SUBROUTINE
      SUBROUTINE  SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )
!
!         SOLVES THE REAL BAND SYSTEM
!            A * X = B  OR  TRANSPOSE(A) * X = B
!         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.
!
!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
!
!     INPUT:
!
!        ABD     REAL(LDA, N)
!                THE OUTPUT FROM SBGCO OR SGBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.
!
!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!
!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SBGCO OR SGBFA.
!
!        B       REAL(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
!                            TRANS(A)  IS THE TRANSPOSE.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
!        OR SGBFA HAS SET INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
!                       FROM FORTRAN: MIN0
!
      implicit none
      INTEGER  LDA, N, ML, MU, IPVT(*), JOB
      REAL*8     ABD(LDA,*), B(*)
!
      REAL*8     T!,SDOT
      INTEGER  K,KB,L,LA,LB,LM,M,NM1
!
!
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .EQ. 0) THEN
!                               ** JOB = 0 , SOLVE  A * X = B
!                               ** FIRST SOLVE L*Y = B
         IF (ML .NE. 0) THEN
            DO 20 K = 1, NM1
               LM = MIN0(ML, N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .NE. K) THEN
                  B(L) = B(K)
                  B(K) = T
               ENDIF
               CALL SAXPY( LM, T, ABD(M+1,K), 1, B(K+1), 1 )
   20       CONTINUE
         ENDIF
!                           ** NOW SOLVE  U*X = Y
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K) / ABD(M,K)
            LM = MIN0(K, M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL SAXPY(LM, T, ABD(LA,K), 1, B(LB), 1)
   40    CONTINUE
!
      ELSE
!                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
!                                  ** FIRST SOLVE  TRANS(U)*Y = B
         DO 60 K = 1, N
            LM = MIN0(K, M) - 1
            LA = M - LM
            LB = K - LM
            T = SDOT(LM, ABD(LA,K), 1, B(LB), 1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
!                                  ** NOW SOLVE TRANS(L)*X = Y
         IF (ML .NE. 0) THEN
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML, N-K)
               B(K) = B(K) + SDOT(LM, ABD(M+1,K), 1, B(K+1), 1)
               L = IPVT(K)
               IF (L .NE. K) THEN
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
               ENDIF
   80       CONTINUE
         ENDIF
!
      ENDIF
!
      RETURN
      END SUBROUTINE
      REAL*8 FUNCTION  SASUM( N, SX, INCX )
!
!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
!
! --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
!
      implicit none
      REAL*8 SX(*)
      integer N,incx,I,M
!
!
      SASUM = 0.0
      IF( N.LE.0 )  RETURN
      IF( INCX.NE.1 ) THEN
!                                          ** NON-UNIT INCREMENTS
          DO 10 I = 1, 1+(N-1)*INCX, INCX
             SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      ELSE
!                                          ** UNIT INCREMENTS
         M = MOD(N,6)
         IF( M.NE.0 ) THEN
!                             ** CLEAN-UP LOOP SO REMAINING VECTOR
!                             ** LENGTH IS A MULTIPLE OF 6.
            DO 30  I = 1, M
              SASUM = SASUM + ABS(SX(I))
   30       CONTINUE
         ENDIF
!                              ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 6
           SASUM = SASUM + ABS(SX(I))   + ABS(SX(I+1)) + ABS(SX(I+2)) &
                         + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50    CONTINUE
      ENDIF
!
      RETURN
      END FUNCTION
      SUBROUTINE     SAXPY( N, SA, SX, INCX, SY, INCY )
!
!          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)
!
!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
!
! --OUTPUT--
!       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH
!                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
!            WHERE LX = 1          IF INCX .GE. 0,
!                     = (-INCX)*N  IF INCX .LT. 0
!            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
!
      implicit none
      REAL*8 SX(*), SY(*), SA
      integer N,incx,incy,M,I,IX,IY
!
!
      IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN
!
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
!
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SY(I) = SY(I) + SA * SX(I)
   10     CONTINUE
!
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
!
!                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,4)
         IF( M .NE. 0 ) THEN
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
              SY(I) = SY(I) + SA * SX(I)
   20       CONTINUE
         ENDIF
!                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 4
            SY(I)   = SY(I)   + SA * SX(I)
            SY(I+1) = SY(I+1) + SA * SX(I+1)
            SY(I+2) = SY(I+2) + SA * SX(I+2)
            SY(I+3) = SY(I+3) + SA * SX(I+3)
   30    CONTINUE
!
      ELSE
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            SY(IY) = SY(IY) + SA*SX(IX)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
!
      ENDIF
!
      RETURN
      END SUBROUTINE
      REAL*8 FUNCTION  SDOT( N, SX, INCX, SY, INCY )
!
!          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'
!
!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
!
! --OUTPUT--
!     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
!            WHERE  LX = 1          IF INCX .GE. 0,
!                      = (-INCX)*N  IF INCX .LT. 0,
!            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
!
      implicit none
      integer N,incx,incy,i,ix,iy,m
      REAL*8 SX(*), SY(*)
!
!
      SDOT = 0.0
      IF( N.LE.0 )  RETURN
!
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
!
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SDOT = SDOT + SX(I) * SY(I)
   10     CONTINUE
!
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
!
!                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,5)
         IF( M .NE. 0 ) THEN
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 4.
            DO 20  I = 1, M
               SDOT = SDOT + SX(I) * SY(I)
   20       CONTINUE
         ENDIF
!                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 5
            SDOT = SDOT + SX(I)*SY(I)     + SX(I+1)*SY(I+1) &
                        + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) &
                        + SX(I+4)*SY(I+4)
   30    CONTINUE
!
      ELSE
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            SDOT = SDOT + SX(IX) * SY(IY)
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
!
      ENDIF
!
      RETURN
      END FUNCTION
      SUBROUTINE     SSCAL( N, SA, SX, INCX )
!
!         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)
!
!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
!            SA  SINGLE PRECISION SCALE FACTOR
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
!
! --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX)
!                FOR I = 0 TO N-1
!
      implicit none
      integer n,incx,i,M
      REAL*8 SA, SX(*)
!
!
      IF( N.LE.0 ) RETURN
!
      IF( INCX.NE.1 ) THEN
!
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             SX(I) = SA * SX(I)
   10     CONTINUE
!
      ELSE
!
         M = MOD(N,5)
         IF( M.NE.0 ) THEN
!                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                           ** IS A MULTIPLE OF 5.
            DO 30  I = 1, M
               SX(I) = SA * SX(I)
   30       CONTINUE
         ENDIF
!                             ** UNROLL LOOP FOR SPEED
         DO 50  I = M+1, N, 5
            SX(I)   = SA * SX(I)
            SX(I+1) = SA * SX(I+1)
            SX(I+2) = SA * SX(I+2)
            SX(I+3) = SA * SX(I+3)
            SX(I+4) = SA * SX(I+4)
   50    CONTINUE
!
      ENDIF
!
      RETURN
      END SUBROUTINE
      SUBROUTINE     SSWAP( N, SX, INCX, SY, INCY )
!
!          INTERCHANGE S.P VECTORS  X  AND  Y
!
!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
!
! --OUTPUT--
!       SX  INPUT VECTOR SY (UNCHANGED IF N .LE. 0)
!       SY  INPUT VECTOR SX (UNCHANGED IF N .LE. 0)
!
!     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
!     WHERE LX = 1          IF INCX .GE. 0,
!              = (-INCX)*N  IF INCX .LT. 0
!     AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
!
      implicit none
      integer n,incx,incy,ix,iy,m,i
      REAL*8 SX(*), SY(*), STEMP1, STEMP2, STEMP3
!
!
      IF( N.LE.0 ) RETURN
!
      IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN
!
          DO 10  I = 1, 1+(N-1)*INCX, INCX
             STEMP1 = SX(I)
             SX(I) = SY(I)
             SY(I) = STEMP1
   10     CONTINUE
!
      ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN
!
!                                        ** EQUAL, UNIT INCREMENTS
         M = MOD(N,3)
         IF( M .NE. 0 ) THEN
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 3.
            DO 20  I = 1, M
               STEMP1 = SX(I)
               SX(I) = SY(I)
               SY(I) = STEMP1
   20       CONTINUE
         ENDIF
!                              ** UNROLL LOOP FOR SPEED
         DO 30  I = M+1, N, 3
            STEMP1  = SX(I)
            STEMP2  = SX(I+1)
            STEMP3  = SX(I+2)
            SX(I)   = SY(I)
            SX(I+1) = SY(I+1)
            SX(I+2) = SY(I+2)
            SY(I)   = STEMP1
            SY(I+1) = STEMP2
            SY(I+2) = STEMP3
   30    CONTINUE
!
      ELSE
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
         IX = 1
         IY = 1
         IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
         IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
         DO 40  I = 1, N
            STEMP1 = SX(IX)
            SX(IX) = SY(IY)
            SY(IY) = STEMP1
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE
!
      ENDIF
!
      RETURN
      END SUBROUTINE
      INTEGER FUNCTION  ISAMAX( N, SX, INCX )
!
!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
!
! --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
!                         ABS(SX(1+(I-1)*INCX))
!
      implicit none
      integer n,incx,i,ii
      REAL*8 SX(*), SMAX, XMAG
!
!
      IF( N.LE.0 ) THEN
         ISAMAX = 0
      ELSE IF( N.EQ.1 ) THEN
         ISAMAX = 1
      ELSE
         SMAX = 0.0
         II = 1
         DO 20  I = 1, 1+(N-1)*INCX, INCX
            XMAG = ABS(SX(I))
            IF( SMAX.LT.XMAG ) THEN
               SMAX = XMAG
               ISAMAX = II
            ENDIF
            II = II + 1
   20    CONTINUE
      ENDIF
!
      RETURN
      END FUNCTION



! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! RCS version control information:
! $Header: D1MACH.f,v 1.2 97/03/18 17:05:25 wiscombe Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      DOUBLE PRECISION FUNCTION D1MACH(I)
      REAL*8 FUNCTION D1MACH(I)  
!  Double-precision machine constants (see R1MACH for documentation).

!  By default, returns values appropriate for a computer with IEEE 
!  arithmetic.  This is an abbreviated version of a routine widely
!  used for 20+ years by numerical analysts.  Most of the values in
!  the original version pertain to computers which went to computer
!  heaven years ago and are of little if any interest.
! 
!  If the values herein do not work for any reason, just look in
!  your Fortran manual for the correct values (usually in the part
!  discussing representations of numbers) and insert them. The exact
!  values are not that important; they can be a factor of 2-3 off
!  without causing any harm.

!  Only I = 1,2,4 is actually used by DISORT. 

!  This routine is superseded in Fortran-90 by the intrinsic numeric 
!  inquiry functions HUGE(1.D0), TINY(1.D0), and EPSILON(1.D0).

!  The original version can be found on NetLib (search by name):
!      http://www.netlib.org/
! ====================================================================

      implicit none
      INTEGER   I
      !EXTERNAL  ERRMSG

      IF( I.EQ.1 )  THEN
         D1MACH = 2.3D-308
!        D1MACH = TINY(1.D0)
      ELSE IF( I.EQ.2 )  THEN  
         D1MACH = 1.7D+308
!        D1MACH = HUGE(1.D0)
      ELSE IF( I.EQ.4 )  THEN  
         D1MACH = 2.3D-16
!        D1MACH = EPSILON(1.D0)
      ELSE
         CALL ERRMSG( 'D1MACH--argument incorrect', .TRUE.)
      END IF

      RETURN
      END FUNCTION







      SUBROUTINE  ERRMSG( MESSAG, FATAL )
!
!        print out a warning or error message;  abort if error
!
      implicit none
      LOGICAL       FATAL, ONCE
      CHARACTER*(*) MESSAG
      INTEGER       MAXMSG, NUMMSG
      SAVE          MAXMSG, NUMMSG, ONCE
      DATA NUMMSG / 0 /,  MAXMSG / 100 /,  ONCE / .FALSE. /

      IF ( FATAL )  THEN
         WRITE ( *, '(/,2A)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NUMMSG = NUMMSG + 1
      IF ( NUMMSG.GT.MAXMSG )  THEN
         IF ( .NOT.ONCE )  WRITE ( *,99 )
         ONCE = .TRUE.
      ELSE
         WRITE ( *, '(/,2A)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ENDIF

      RETURN

   99 FORMAT( ///,' >>>>>>  TOO MANY WARNING MESSAGES --  ', &
         'THEY WILL NO LONGER BE PRINTED  <<<<<<<', /// )
      END SUBROUTINE
      
      LOGICAL FUNCTION  WRTBAD ( quiet, VARNAM )
!
!          write names of erroneous variables and return 'true'
!
!      input :   VARNAM = name of erroneous variable to be written
!                         ( character, any length )
! ----------------------------------------------------------------------
      implicit none
      CHARACTER*(*)  VARNAM
      INTEGER        MAXMSG, NUMMSG
      LOGICAL quiet
      SAVE  NUMMSG, MAXMSG
      DATA  NUMMSG / 0 /,  MAXMSG / 50 /


      WRTBAD = .TRUE.
      NUMMSG = NUMMSG + 1
      IF ( .NOT. quiet ) &
          WRITE ( *, '(3A)' )  ' ****  INPUT VARIABLE  ', VARNAM, &
           '  IN ERROR  ****'
      IF ( NUMMSG.EQ.MAXMSG  .AND. .NOT. quiet ) &
       CALL  ERRMSG( 'TOO MANY INPUT ERRORS.  ABORTING...$', .TRUE. )
      RETURN
      END FUNCTION
      
      LOGICAL FUNCTION  WRTDIM ( quiet, DIMNAM, MINVAL )
!
!          write name of too-small symbolic dimension and
!          the value it should be increased to;  return 'true'
!
!      input :  DIMNAM = name of symbolic dimension which is too small
!                        ( character, any length )
!               MINVAL = value to which that dimension should be
!                        increased (at least)
! ----------------------------------------------------------------------
      implicit none
      CHARACTER*(*)  DIMNAM
      INTEGER        MINVAL
      LOGICAL quiet


      IF ( .NOT. quiet ) &
          WRITE ( *, '(3A,I7)' )  ' ****  SYMBOLIC DIMENSION  ', &
          DIMNAM, '  SHOULD BE INCREASED TO AT LEAST ', MINVAL
      WRTDIM = .TRUE.
      RETURN
      END FUNCTION



! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! RCS version control information:
! $Header: R1MACH.f,v 1.2 97/03/18 17:04:13 wiscombe Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REAL FUNCTION R1MACH(I)

!        Single-precision machine constants

!  Assume floating-point numbers are represented in the t-digit,
!  base-b form

!         sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )

!  where 0.le.x(i).lt.b  for  i = 1,...,t,
!  0.lt.x(1), and  emin.LE.e.LE.emax.  then

!  R1MACH(1) = b**(emin-1), the smallest positive magnitude
!              (use TINY(R) in Fortran 90, where R is a single
!              precision variable)

!  R1MACH(2) = b**emax*(1 - b**(-t)), the largest magnitude
!              (use HUGE(R) in Fortran 90, where R is a single
!              precision variable))

!  R1MACH(3) = b**(-t), the smallest relative spacing.

!  R1MACH(4) = b**(1-t), the largest relative spacing.  i.e.,
!              smallest positive eps such that  1+eps .ne. 1
!              (use EPSILON(R) in Fortran 90, where R is a single
!              precision variable))

!  R1MACH(5) = LOG10(b)


!  Reference: Fox P.A., Hall A.D., Schryer N.L.,'Framework For A
!               Portable Library', ACM Transactions On Mathematical
!               Software, Vol. 4, No. 2, June 1978, pp. 177-188.


!  By default, returns values appropriate for a computer with IEEE 
!  arithmetic.  This is an abbreviated version of a routine widely
!  used for 20+ years by numerical analysts.  Most of the values in
!  the original version pertain to computers which went to computer
!  heaven years ago and are of little if any interest.
! 
!  If the values herein do not work for any reason, just look in
!  your Fortran manual for the correct values (usually in the part
!  discussing representations of numbers) and insert them. The exact
!  values are not that important; they can be a factor of 2-3 off
!  without causing any harm.

!  Only I = 1,2,4 is actually used by DISORT. 

!  This routine is superseded in Fortran-90 by the intrinsic numeric 
!  inquiry functions HUGE(1.0), TINY(1.0), and EPSILON(1.0).

!  The original version can be found on NetLib (search by name):
!      http://www.netlib.org/
! ====================================================================
      implicit none
      INTEGER I
      !EXTERNAL  ERRMSG

      IF( I.EQ.1 )  THEN
         R1MACH = 1.2E-38
!        R1MACH = TINY(1.0)
      ELSE IF( I.EQ.2 )  THEN  
         R1MACH = 3.4E+38
!        R1MACH = HUGE(1.0)
      ELSE IF( I.EQ.4 )  THEN  
         R1MACH = 1.2E-07
!        R1MACH = EPSILON(1.0)
      ELSE
         CALL ERRMSG( 'R1MACH--argument incorrect', .TRUE.)
      END IF

      RETURN
      END FUNCTION

  END MODULE twostr_module

