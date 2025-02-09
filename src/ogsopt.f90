SUBROUTINE ogsopt(Optsen)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! LINEAR OPTIMISATION MODE
! ======================================================================
! INP | OPTSEN           | I*4 | SENSITIVITY OPTIMISATION MODE
!     |                  |     | ->  0: NO
!     |                  |     | -> -1: INITIALISATION
!     |                  |     | -> +1: WITH CONSTRAINT CALCULATION
!     |                  |     | -> +2: WITH CONSTRAINT BIAS
!     |                  |     | -> +3: WITH CONSTRAINT CALC / NO OPTIM
!     |                  |     | -> +4: WITH CONSTRAINT BIAS / NO OPTIM
! ======================================================================
! 2021/03/30 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   INTEGER(4) Optsen
! ======================================================================
   Senopt = Optsen
! ======================================================================
END SUBROUTINE ogsopt
