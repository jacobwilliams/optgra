SUBROUTINE ogcstr(Strcon,Lencon)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEFINE CONSTRAINT + MERIT STRING
! ======================================================================
! INP | STRCON(NUMCON+1) | C80 | CONIABLES NAME STRING
! ----------------------------------------------------------------------
! INP | LENCON(NUMCON+1) | I*4 | CONIABLES NAME LENGTH
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   CHARACTER*80 Strcon(Numcon+1)
   INTEGER(4) Lencon(Numcon+1)
! ======================================================================
   INTEGER(4) con , len
! ======================================================================
   DO con = 1 , Numcon + 1
      len = min(Lencon(con),80)
      Constr(con) = Strcon(con)
      Conlen(con) = len
   ENDDO
! ======================================================================
END SUBROUTINE ogcstr
