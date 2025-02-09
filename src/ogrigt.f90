
SUBROUTINE ogrigt(Actinp,Actout)
! ======================================================================
! RIGHT-MULTIPLIES VECTOR LOWER TRIANGULAR MATRIX OBTAINED BY REDUCTION
! AND SUBSEQUENT INVERSION OF DERIVATIVES OF ACTIVE CONSTRAINTS
! ======================================================================
! INP | ACTINP(NUMCON)   | R*8 | VECTOR INITAL
! ----------------------------------------------------------------------
! OUT | ACTOUT(NUMCON)   | R*8 | VECTOR FINAL (MAY BE SAME AS ACTINP)
! ======================================================================
! SUBROUTINES CALLED: NONE
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   REAL(8) Actinp(Numcon)
   REAL(8) Actout(Numcon)
! ======================================================================
   INTEGER(4) row , col , act
   REAL(8) val
! =====================================================================
   DO col = Numact , 1 , -1
      val = Actinp(col)
      DO act = Numact , col + 1 , -1
         row = Actcon(act)
         val = val - Conred(row,col)*Actout(act)
      ENDDO
      row = Actcon(col)
      Actout(col) = val/Conred(row,col)
   ENDDO
! ======================================================================
END SUBROUTINE ogrigt
