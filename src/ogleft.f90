SUBROUTINE ogleft(Actinp,Actout)
! ======================================================================
! LEFT-MULTIPLIES VECTOR LOWER TRIANGULAR MATRIX OBTAINED BY REDUCTION
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
! ======================================================================
   DO act = 1 , Numact
      row = Actcon(act)
      val = Actinp(act)
      DO col = 1 , act - 1
         val = val - Conred(row,col)*Actout(col)
      ENDDO
      Actout(act) = val/Conred(row,act)
   ENDDO
! ======================================================================
END SUBROUTINE ogleft
