SUBROUTINE oggsst(Varsen,Quasen,Consen,Actsen,Dersen,Actsav,Consav,Redsav,Dersav,Actnum)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
! Function to get sensitivity state data, necessary for serialization.
! Do not use this directly except in serialization routines
! ======================================================================
! OUT | VARSEN(NUMVAR)   | I*4 | STORED VARIABLES VALUE
! OUT | QUASEN(NUMCON+1) | R*8 | STORED CONSTRAINTS CORRECTION VECTOR
! OUT | CONSEN(NUMCON+1) | R*8 | STORED CONSTRAINTS VALUE
! OUT | ACTSEN(NUMCON+1) | R*8 | STORED CONSTRAINTS ACTIVE
! OUT | DERSEN(NUMCON+1, | R*8 | STORED DERIVATIVE
!                NUMVAR) |     |
! OUT | ACTSAV(NUMCON+1) | I*4 | STORED ACTIVE CONSTRAINTS
! OUT | CONSAV(NUMCON+4) | I*4 | STORED ACTIVE CONSTRAINTS
! OUT | REDSAV(NUMCON+3, | R*8 | STORED DERIVATIVE
!                NUMVAR) |     |
! OUT | DERSAV(NUMCON+3, | R*8 | STORED DERIVATIVE
!                NUMVAR) |     |
! OUT | ACTNUM           | I*4 | NUMBER OF ACTIVE CONSTRAINTS
! ======================================================================
! 2021/07/19 | M. von Looz | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   REAL(8) Varsen(Numvar)
   REAL(8) Quasen(Numcon+1)
   REAL(8) Consen(Numcon+1)
   INTEGER(4) Actsen(Numcon+1)
   REAL(8) Dersen(Numcon+1,Numvar)
   INTEGER(4) Actsav(Numcon+1)
   INTEGER(4) Consav(Numcon+4)
   REAL(8) Redsav(Numcon+3,Numvar)
   REAL(8) Dersav(Numcon+3,Numvar)
! ======================================================================
   INTEGER(4) Actnum
   INTEGER(4) var , con
! ======================================================================
! Variable values saved for sensitivity
! ----------------------------------------------------------------------
   Actnum = Numact

   DO var = 1 , Numvar
      Varsen(var) = Senvar(var)
   ENDDO

   DO con = 1 , Numcon + 1
      Quasen(con) = Senqua(con)
      Consen(con) = Sencon(con)
      Actsen(con) = Senact(con)
      DO var = 1 , Numvar
         Dersen(con,var) = Sender(con,var)
      ENDDO
   ENDDO
! ======================================================================
! Temporary status saved of which constraints are active
! ----------------------------------------------------------------------
   DO con = 1 , Numcon + 1
      Actsav(con) = Actcon(con)
   ENDDO

   DO con = 1 , Numcon + 4
      Consav(con) = Conact(con)
   ENDDO

   DO con = 1 , Numcon + 3
      DO var = 1 , Numvar
         Redsav(con,var) = Conred(con,var)
         Dersav(con,var) = Conder(con,var)
      ENDDO
   ENDDO
! ======================================================================
END SUBROUTINE oggsst
