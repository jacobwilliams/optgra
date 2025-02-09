SUBROUTINE ogderi(Dervar,Pervar)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEFINE COMPUTATION OF DERIVATIVES
! ======================================================================
! INP | DERVAR           | I*4 | DERIVATIVES COMPUTATION MODE
!     |                  |     | -> 1: USER DEFINED
!     |                  |     | -> 2: NUMERIC WITH DOUBLE DIFFERENCING
!     |                  |     | -> 3: NUMERIC WITH SINGLE DIFFERENCING
! ----------------------------------------------------------------------
! INP | PERVAR(NUMVAR)   | R*8 | VARIABLES PERTURBATION FOR DERIVATIVES
!     |                  |     | -> NOT SCALED
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   INTEGER(4) Dervar
   REAL(8) Pervar(Numvar)
! ======================================================================
   INTEGER(4) var
! ======================================================================
   Varder = Dervar
! ----------------------------------------------------------------------
   DO var = 1 , Numvar
      Varper(var) = Pervar(var)
   ENDDO
! ======================================================================
END SUBROUTINE ogderi
