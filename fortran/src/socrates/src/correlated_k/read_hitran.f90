! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read in HITRAN line data
!
SUBROUTINE read_hitran ( &
     iu_lbl,          & ! in
     iu_monitor,      & ! in
     molecule,        & ! in
     nd_lines,        & ! in
     n_isotope,       & ! in
     isotope,         & ! in
     wavenum_min,     & ! in
     wavenum_max,     & ! in
     number_of_lines, & ! out
     hitran_data )      ! out
  
  
! Description:
!
! Reads in HITRAN line data for required molecular species across
! required spectral interval into structure defined above
! 
! Method:
!
! Straightforward
!
! Modules used:
  USE realtype_rd
  USE def_hitran_record
  USE hitran_cnst
  
  IMPLICIT NONE
  
  ! Subroutine arguments
  
  ! Scalar arguments with intent(in):
  INTEGER, INTENT(IN) :: nd_lines     ! Size allocated for number of lines
  INTEGER, INTENT(IN) :: molecule     ! molecule number
  INTEGER, INTENT(IN) :: n_isotope    ! Number of isotopes
  INTEGER, INTENT(IN) :: isotope      ! isotope number
  INTEGER, INTENT(IN) :: iu_lbl       ! unit number for I/O from LbL file
  INTEGER, INTENT(IN) :: iu_monitor   ! unit number for monitoring output
  REAL (RealK), INTENT(IN) :: wavenum_min  ! in cm-1
  REAL (RealK), INTENT(IN) :: wavenum_max  ! in cm-1
  
  ! Scalar arguments with intent(out):
  INTEGER, INTENT(OUT) :: number_of_lines
  
  ! Array arguments with intent(out):
  TYPE(StrHitranRec), DIMENSION(nd_lines), INTENT(OUT) :: hitran_data
  
  ! Local parameters:
  CHARACTER (LEN = *), PARAMETER :: routine_name = "read_hitran"
  
  INTEGER, PARAMETER :: input_unit  = 10
  
  ! Local scalars:
  TYPE(StrHitranRec) :: single_record
  
  
  INTEGER :: count_lines
  INTEGER :: io_status    ! Error code for file I/O
  
  REAL :: upper_cutoff
  REAL :: lower_cutoff
  
!- End of header
  
  WRITE(*,"(a,a6)")          "Gas required: ", molecule_names(molecule)
  WRITE(iu_monitor,"(a,a6)") "Gas required: ", molecule_names(molecule)
  
  WRITE(*,"(a,f12.6,2x,f12.6)") &
    "Band limits: ", wavenum_min, wavenum_max
  WRITE(iu_monitor,"(a,f12.6,2x,f12.6)") &
    "Band limits: ", wavenum_min, wavenum_max
  
  
  
    WRITE(*,"(a)")          "Opened HITRAN data file "
    WRITE(iu_monitor,"(a)") "Opened HITRAN data file "
  
    upper_cutoff = wavenum_max + (wavenum_max-wavenum_min)/10.0
    lower_cutoff = wavenum_min - (wavenum_max-wavenum_min)/10.0
  
    number_of_lines = 0
    count_lines     = 0
  
    DO
      count_lines = count_lines + 1
      io_status   = 0
  
      READ(iu_lbl, hitran_record_frmt, IOSTAT = io_status) single_record
      IF (io_status == 0) THEN
        IF (single_record % frequency <= upper_cutoff) THEN
          IF (single_record % frequency >= lower_cutoff) THEN
            IF (single_record % mol_num == molecule) THEN
!             The strengths of lines in HITRAN are scaled by the
!             terrestrial abundance. For most applications we read
!             all lines for all isotopes. The following test should
!             be uncommented only when it is desired to read one 
!             specific isotope.
!             IF (single_record % iso_num == isotope) THEN              
                number_of_lines              = number_of_lines + 1
                hitran_data(number_of_lines) = single_record
!             ENDIF
            ENDIF
          ENDIF
        ELSE
          EXIT
        ENDIF
      ELSE IF (io_status < 0) THEN
        EXIT
      ELSE
        WRITE(*,"(a,i6)")          &
             "Error reading HITRAN file at line: ", count_lines 
        WRITE(iu_monitor,"(a,i6)") &
             "Error reading HITRAN file at line: ", count_lines 
        EXIT
      ENDIF
  
    ENDDO
  
END SUBROUTINE read_hitran
