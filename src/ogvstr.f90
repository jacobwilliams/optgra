SUBROUTINE ogvstr(Strvar,Lenvar)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEFINE VARIABLE STRING
! ======================================================================
! INP | STRVAR(NUMVAR)   | C80 | VARIABLES NAME STRING
! ----------------------------------------------------------------------
! INP | LENVAR(NUMVAR)   | I*4 | VARIABLES NAME LENGTH
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   CHARACTER*80 Strvar(Numvar)
   INTEGER(4) Lenvar(Numvar)
! ======================================================================
   INTEGER(4) var , len
! ======================================================================
   DO var = 1 , Numvar
      len = min(Lenvar(var),80)
      Varstr(var) = Strvar(var)
      Varlen(var) = len
   ENDDO
! ======================================================================
END SUBROUTINE ogvstr
