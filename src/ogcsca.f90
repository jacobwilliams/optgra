SUBROUTINE ogcsca(Scacon)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEFINE CONSTRAINT + MERIT CONVERGENCE THRESHOLDS
! ======================================================================
! INP | SCACON(NUMCON+1) | R*8 | CONSTRAINTS CONVER THRESHOLD (1:NUMCON)
!     |                  |     | MERIT       CONVER THRESHOLD (1+NUMCON)
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   REAL(8) Scacon(Numcon+1)
! ======================================================================
   INTEGER(4) con
! ======================================================================
   DO con = 1 , Numcon + 1
      Consca(con) = Scacon(con)
   ENDDO
! ======================================================================
END SUBROUTINE ogcsca
