SUBROUTINE ogomet(Metopt)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEFINE OPTIMISATION CONTROL PARAMETERS
! ======================================================================
! INP | METOPT | I*4 | OPTIMISATION METHOD
!     |        |     | 3: CONJUGATE GRADIENT METHOD
!     |        |     | 2: SPETRAL CONJUGATE GRADIENT METHOD
!     |        |     | 1: MODIFIED SPETRAL CONJUGATE GRADIENT METHOD
!     |        |     | 0: STEEPEST DESCENT METHOD
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   INTEGER(4) Metopt
! ======================================================================
   Optmet = Metopt
! ======================================================================
END SUBROUTINE ogomet
