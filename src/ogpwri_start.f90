
SUBROUTINE ogpwri_start()
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! WRITE OPTIMIZATION LOG IN PYGMO FORMAT
! ======================================================================
! INP | VARDER           | I*4 | DERIVATIVES COMPUTATION MODE
!     |                  |     | -> 0: VALUES ONLY
!     |                  |     | -> 1: USER DEFINED
!     |                  |     | -> 2: NUMERIC WITH DOUBLE DIFFERENCING
!     |                  |     | -> 3: NUMERIC WITH SINGLE DIFFERENCING
! ======================================================================
! 2023/01/25 | W. MARTENS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
! ======================================================================
   WRITE (Loglup,'("OPTGRA plugin for pagmo/pygmo:")')
   IF ( Varder==0 ) THEN
      WRITE (Loglup,'("")')
   ELSEIF ( Varder==1 .OR. Varder==-1 ) THEN
      WRITE (Loglup,'("    User-defined gradients")')
   ELSEIF ( Varder==2 ) THEN
      WRITE (Loglup,'("    Numerical gradients by double differencing")')
   ELSEIF ( Varder==3 ) THEN
      WRITE (Loglup,'("    Numerical gradients by single differencing")')
   ENDIF

   IF ( Optmet==3 ) THEN
      WRITE (Loglup,'("    Conjugate gradient method")')
   ELSEIF ( Optmet==2 ) THEN
      WRITE (Loglup,'("    Spectral conjugate gradient method")')
   ELSEIF ( Optmet==1 ) THEN
      WRITE (Loglup,'("    Modified spectral conjugate gradient method")')
   ELSEIF ( Optmet==0 ) THEN
      WRITE (Loglup,'("    Steepest descent method")')
   ENDIF

   WRITE (Loglup,'("")')

! ======================================================================
END SUBROUTINE ogpwri_start
