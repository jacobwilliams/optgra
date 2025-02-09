SUBROUTINE ogdist(Maxvar,Sndvar)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEFINE OPTIMISATION CONTROL PARAMETERS
! ======================================================================
! INP | ITEMAX           | I*4 | MAXIMUM NUMBER OF ITERATIONS
! ----------------------------------------------------------------------
! INP | MAXVAR           | R*8 | MAXIMUM DISANCE PER ITERATION
!     |                  |     | -> SCALED
! ----------------------------------------------------------------------
! INP | SNDVAR           | R*8 | PERTURBATION FOR SND ORDER DERIVATIVES
!     |                  |     | -> SCALED
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   REAL(8) Maxvar
   REAL(8) Sndvar
! ======================================================================
   Varmax = Maxvar
   Varsnd = Sndvar
! ======================================================================
END SUBROUTINE ogdist
