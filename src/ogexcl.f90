SUBROUTINE ogexcl(Exc)
! ======================================================================
! REMOVE CONSTRAINT TO ACTIVE SET AND REDUCES DERIVATIVES
! ======================================================================
! INP | EXC              | I*4 | CONSTRAINT TO BE REMOVED
!     |                  |     | SEQUENCE NUMBER IN ACTIVE LIST
! ======================================================================
! SUBROUTINES CALLED: NONE
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   INTEGER(4) Exc
! ======================================================================
   REAL(8) val , bet , gam
   INTEGER(4) row , col , act , con
   CHARACTER str*256
! ======================================================================
! ADJUST LIST OF ACTIVE CONSTRAINTS
! ----------------------------------------------------------------------
   con = Actcon(Exc)
   Conact(con) = 0
   Numact = Numact - 1
   DO act = Exc , Numact
      con = Actcon(act+1)
      Actcon(act) = con
      Conact(con) = Conact(con) - 1
   ENDDO
! ======================================================================
! REDUCE FOR SUBSEQUENT CONSTRAINTS
! ----------------------------------------------------------------------
   DO act = Exc , Numact
      con = Actcon(act)
      val = 0D0
      DO col = act , act + 1
         val = val + Conred(con,col)**2
      ENDDO
      val = dsqrt(val)
      IF ( Conred(con,act)>0D0 ) val = -val
      IF ( dabs(val)<1D-15 ) THEN
         WRITE (Loglun,*) "OGEXCL-ERROR: CONSTRAINTS SINGULAR"
         CALL ogwrit(2,str)
         WRITE (Loglun,*) "VAL=" , val
         CALL ogwrit(2,str)
         STOP
      ENDIF
      Conred(con,act) = Conred(con,act) - val
      bet = 1D0/(val*Conred(con,act))
      DO row = 1 , Numcon + 3
         IF ( Conact(row)>act .OR. Conact(row)<=0 ) THEN
            gam = 0D0
            DO col = act , act + 1
               IF ( Conred(row,col)/=0D0 ) gam = gam + Conred(row,col)*Conred(con,col)
            ENDDO
            IF ( gam/=0D0 ) THEN
               gam = gam*bet
               DO col = act , act + 1
                  Conred(row,col) = Conred(row,col) + Conred(con,col)*gam
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      Conred(con,act) = val
      DO col = act + 1 , act + 1
         Conred(con,col) = 0D0
      ENDDO
   ENDDO
! ======================================================================
END SUBROUTINE ogexcl
