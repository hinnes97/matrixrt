! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE astro_constants_mod

! Description: Parameters of the Earth's orbit.
!
!     Default values of the orbital elements
!      (currently those for the epoch J2000 which is 1.5d Jan. 2000):
!     The Eccentricity and Longitue of perhelion are recommended by NASA
!      see (http://ssd.jpl.nasa.gov/elem_planets.html)
!     The Obliquity value comes from the Astronomical Almanac for 1984
!      page S26 and is used on several webpages e.g.
!      nedwww.ipac.caltech.edu/help/calc_doc.txt
!      www.stargazing.net/kepler/astrovba2.html
!      http://edhs1.gsfc.nasa.gov/waisdata/docsw/txt/tp4450505.txt
!
!     The data in the series expansions are adapted from Berger 1978.
!      Andre' L. Berger Journ. of Atm. Sci. Volume 35 p. 2362-2367,
!      and is available from the Met Office library.

USE rad_ccf, ONLY: pi

IMPLICIT NONE

! Eccentricity of the orbit
REAL,    PARAMETER :: E_DFLT         = 1.6710222E-02

! Longitude of the perihelion in radians
REAL,    PARAMETER :: LPH_DFLT       = 102.94719*pi/180.0

! Obliquity of the orbit - corresponds to 23.43929111 degrees
REAL,    PARAMETER :: OBLQ_DFLT      = 0.409092804

! Reference year for setting the date of the vernal equinox
INTEGER, PARAMETER :: YEAR_REF_VE    = 2000

! Date of the vernal equinox in days after the start of the year
! This date is for the year 2000.
REAL,    PARAMETER :: DATE_VE_DFLT   = 79.3159

!     The final parameter required is the time of the perihelion
!     passage, TAU0. For a pure Keplerian orbit, with a specified
!     eccentricity and longitude of the perihelion, this can be
!     deduced from the date of the vernal equinox (as is specified
!     in AMIP-2, for example). In practice it is somewhat more
!     complicated to calculate the time of the perihelion.
!     For simplicity, a mean value for the years 1995-2005 is used
!     here: note that the range of TAU0 in this period is from
!     1.0 to 3.75 and that there is no simple relationship with leap
!     years.

! Time of the perihelion passage in days
REAL,    PARAMETER :: TAU0_DFLT      = 2.667

!     The tropical year is defined as the mean interval between two
!     successive passages of the sun through the vernal equinox.
!     This value is 365.2424 and will be used here.
!     The value 365.2422 most often used is the mean length of
!     the year starting at different points on the ellipse.

! Number of days in the tropical year
REAL, PARAMETER    :: TropYearLength = 365.2424

!     ------------------------------------------------------------------

!     The parameters used to calculate secular variations of the orbital
!     elements are taken from A. L. Berger (1978), J. Atm. Sci., vol 35,
!     p. 2362.These have been converted so that:
!     amplitudes are in radians,
!     angular frequencies in radians per year and
!     phases in radians
!     Phases have also been converted to be taken relative to
!     J2000 for consistency with the time for the default values of the
!     orbital parameters.
!     Berger's numbers (with time correction) differ slightly from the
!     default values above.
!
!     The obliquity and longitude of the perihelion have been adjusted
!     to agree with the NASA values for J2000, but adjustment of the
!     eccentricity to NASA values is not so easy and has not been done.


! Reference year     YEAR_REF
INTEGER, PARAMETER :: YEAR_REF       = 2000

!  -----------------------------------------------------------------
!     Obliquity (Table 1) from the first 24 terms:
!     (enough for deviations of less than 0.002 degrees)

! Constant term in the obliquity: from the Astr. Almanac for 1984
! The following value corresponds to 23.320870 degrees at J2000
REAL,    PARAMETER :: OBLQ_CNST      = 0.40702597
!
! Number of terms retained in the series for the obliquity
INTEGER, PARAMETER :: N_TERM_OBQ     = 24
!
! Amplitude
REAL, PARAMETER :: A(N_TERM_OBQ) = (/                                   &
          -1.19372E-02, -4.15640E-03, -3.05103E-03, -2.00849E-03        &
        , -1.51146E-03,  1.49778E-03, -7.88065E-04, -5.62917E-04        &
        ,  4.90244E-04, -3.28170E-04,  1.20767E-04,  1.09471E-04        &
        , -1.02587E-04, -7.58733E-05,  7.46128E-05,  7.11222E-05        &
        , -5.68686E-05,  4.97904E-05,  3.14644E-05,  2.83616E-05        &
        , -2.66163E-05, -2.63254E-05,  2.50164E-05,  2.46285E-05 /)
! Angular frequency
REAL, PARAMETER :: F(N_TERM_OBQ) = (/                                   &
           1.5324946E-04,  1.5814864E-04,  1.1719011E-04                &
        ,  1.5506174E-04,  2.1733392E-04,  1.5016256E-04                &
        ,  2.1170962E-04,  1.5633636E-04,  1.4835028E-04                &
        ,  2.0692488E-04,  2.1252514E-04,  2.2999289E-04                &
        ,  3.0649899E-04,  3.1139817E-04,  4.8991877E-06                &
        ,  3.6059331E-05,  2.7043965E-04,  1.8122966E-06                &
        ,  6.4084427E-05,  3.0341210E-04,  3.0831127E-04                &
        ,  3.7058338E-04,  2.2211866E-04,  4.0958519E-05 /)
! Phase in the series
REAL, PARAMETER :: D(N_TERM_OBQ) = (/                                   &
           4.4041E+00,  4.9093E+00,  2.2451E+00,  5.1167E+00            &
        ,  2.7912E-01,  4.6115E+00,  5.3935E+00,  4.1966E+00            &
        ,  3.8990E+00,  4.7014E+00,  5.5397E+00,  5.5896E+00            &
        ,  2.5251E+00,  3.0303E+00,  5.0517E-01,  2.1589E+00            &
        ,  3.6608E-01,  7.1253E-01,  2.1582E+00,  2.7325E+00            &
        ,  3.2376E+00,  4.6833E+00,  9.7121E-01,  2.6640E+00 /)

!   -----------------------------------------------------------------
!     Eccentricity and longitude of the fixed perihelion (Table 4):

! Number of terms retained in the series for the
! eccentricty and longitude of the perihelion
INTEGER, PARAMETER :: N_TERM_ECN_LPH = 19
!
! Amplitude
REAL, PARAMETER :: M(N_TERM_ECN_LPH)= (/                                &
           1.8607980E-02,  1.6275220E-02, -1.3006600E-02                &
        ,  9.8882900E-03, -3.3670000E-03,  3.3307700E-03                &
        , -2.3540000E-03,  1.4001500E-03,  1.0070000E-03                &
        ,  8.5700000E-04,  6.4990000E-04,  5.9900000E-04                &
        ,  3.7800000E-04, -3.3700000E-04,  2.7600000E-04                &
        ,  1.8200000E-04, -1.7400000E-04, -1.2400000E-04                &
        ,  1.2500000E-05 /)
! Angular frequency
REAL, PARAMETER :: G(N_TERM_ECN_LPH)= (/                                &
           2.0397105E-05,  3.5614854E-05,  8.6574454E-05                &
        ,  8.3487563E-05,  8.1675266E-05,  2.5205846E-05                &
        ,  8.8386751E-05,  1.2710243E-04,  3.0830121E-05                &
        ,  7.8588375E-05,  1.4860417E-05,  8.0400672E-05                &
        ,  8.9661345E-05,  3.0014587E-05,  9.1473642E-05                &
        ,  8.4481533E-05,  2.9990579E-05,  8.9290274E-05                &
        ,  3.2378912E-06 /)
! Phase in the series
REAL, PARAMETER :: B(N_TERM_ECN_LPH) = (/                               &
           5.0053E-01,  3.3839E+00,  5.3852E+00,  5.5925E+00            &
        ,  4.8800E+00,  1.5230E+00,  6.0977E+00,  2.2481E+00            &
        ,  2.6918E+00,  5.0874E+00,  2.0054E+00,  5.8001E+00            &
        ,  5.1778E+00,  2.5455E+00,  5.8903E+00,  2.6587E+00            &
        ,  2.2151E+00,  3.6812E+00,  1.2585E+00 /)
!   ------------------------------------------------------------------
!     General Precession (Table 5):

! Linear rate of precession!
! The value corresponds to 50.439273 seconds per year -Berger 1979
REAL,    PARAMETER :: LIN_RATE_GN_PRCS  = 2.44536496E-04

! Constant offset to general precession (in seconds pre year),
! corrected for 50 years difference in reference time.
REAL,    PARAMETER ::  GN_PRCS_CNST      = 7.14372244E-02

! Number of terms kept in the series for the general precession
INTEGER, PARAMETER :: N_TERM_GN_PRCS = 10

! Amplitude
REAL, PARAMETER :: C(N_TERM_GN_PRCS) = (/                               &
           3.58327E-02,  1.23877E-02,  9.80662E-03, -9.56853E-03        &
        ,  6.01280E-03,  4.62449E-03, -4.51725E-03,  4.22942E-03        &
        ,  2.93967E-03, -2.40482E-03 /)
! Angular frequency
REAL, PARAMETER :: H(N_TERM_GN_PRCS) = (/                               &
           1.5324946E-04,  1.5814864E-04,  1.1719011E-04,  3.0868911E-06&
        ,  1.5506174E-04,  1.5217749E-05,  1.5016256E-04,  2.1733392E-04&
        ,  4.8087409E-06,  1.8122966E-06 /)
! Phase in the series
REAL, PARAMETER :: R(N_TERM_GN_PRCS) = (/                               &
           4.4041E+00,  4.9093E+00,  2.2451E+00,  6.0756E+00            &
        ,  5.1167E+00,  2.8833E+00,  4.6115E+00,  2.7912E-01            &
        ,  1.0225E+00,  7.1253E-01 /)
!

END MODULE astro_constants_mod
