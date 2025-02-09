SUBROUTINE ogssst(Varsen,Quasen,Consen,Actsen,Dersen,Actsav,Consav,Redsav,Dersav,Actnum)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
! Function to get sensitivity state data, necessary for serialization.
! Do not use this directly except in serialization routines
! ======================================================================
! INP | VARSEN(NUMVAR)   | I*4 | STORED VARIABLES VALUE
! INP | QUASEN(NUMCON+1) | R*8 | STORED CONSTRAINTS CORRECTION VECTOR
! INP | CONSEN(NUMCON+1) | R*8 | STORED CONSTRAINTS VALUE
! INP | ACTSEN(NUMCON+1) | R*8 | STORED CONSTRAINTS ACTIVE
! INP | DERSEN(NUMCON+1, | R*8 | STORED DERIVATIVE
!                NUMVAR) |     |
! INP | ACTSAV(NUMCON+1) | I*4 | STORED ACTIVE CONSTRAINTS
! INP | CONSAV(NUMCON+4) | I*4 | STORED ACTIVE CONSTRAINTS
! INP | REDSAV(NUMCON+3, | R*8 | STORED DERIVATIVE
!                NUMVAR) |     |
! INP | DERSAV(NUMCON+3, | R*8 | STORED DERIVATIVE
!                NUMVAR) |     |
! INP | ACTNUM           | I*4 | NUMBER OF ACTIVE CONSTRAINTS
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
   Numact = Actnum

   DO var = 1 , Numvar
      Senvar(var) = Varsen(var)
   ENDDO

   DO con = 1 , Numcon + 1
      Senqua(con) = Quasen(con)
      Sencon(con) = Consen(con)
      Senact(con) = Actsen(con)
      DO var = 1 , Numvar
         Sender(con,var) = Dersen(con,var)
      ENDDO
   ENDDO
! ======================================================================
! Temporary status saved of which constraints are active
! ----------------------------------------------------------------------
   DO con = 1 , Numcon + 1
      Actcon(con) = Actsav(con)
   ENDDO

   DO con = 1 , Numcon + 4
      Conact(con) = Consav(con)
   ENDDO

   DO con = 1 , Numcon + 3
      DO var = 1 , Numvar
         Conred(con,var) = Redsav(con,var)
      ENDDO
   ENDDO

   DO con = 1 , Numcon
      DO var = 1 , Numvar
         Conder(con,var) = Dersav(con,var)
      ENDDO
   ENDDO
! ======================================================================
END SUBROUTINE ogssst
