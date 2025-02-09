SUBROUTINE ogincl(Inc)
! ======================================================================
! ADDS CONSTRAINT TO ACTIVE SET AND REDUCES DERIVATIVES
! ======================================================================
! INP | INC              | I*4 | CONSTRAINT TO BE INCLUDED
! ======================================================================
! SUBROUTINES CALLED: NONE
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   INTEGER(4) Inc
! ======================================================================
   REAL(8) val , fac , gam , sav , max
   INTEGER(4) row , col , ind , lst
   CHARACTER str*256
! ======================================================================
! GENERAL
! ----------------------------------------------------------------------
   Numact = Numact + 1
! ======================================================================
! PERMUTATION TO GET ZERO DERIVATIVES AT END FOR NEW ACTIVE CONSTRAINT
! ----------------------------------------------------------------------
   lst = Numvar
   DO col = Numvar , Numact , -1
      IF ( Conred(Inc,col)==0D0 ) THEN
         IF ( col/=lst ) THEN
            DO row = 1 , Numcon + 3
               IF ( Conact(row)<=0 ) THEN
                  sav = Conred(row,col)
                  Conred(row,col) = Conred(row,lst)
                  Conred(row,lst) = sav
               ENDIF
            ENDDO
         ENDIF
         lst = lst - 1
      ENDIF
   ENDDO
! ======================================================================
! PERMUTATION TO GET MAXIMUM PIVOT
! ----------------------------------------------------------------------
   ind = Numact
   max = dabs(Conred(Inc,ind))
   DO col = Numact + 1 , lst
      val = dabs(Conred(Inc,col))
      IF ( val>max ) THEN
         ind = col
         max = val
      ENDIF
   ENDDO
! ----------------------------------------------------------------------
   IF ( ind/=Numact ) THEN
      DO row = 1 , Numcon + 3
         IF ( Conact(row)<=0 ) THEN
            sav = Conred(row,ind)
            Conred(row,ind) = Conred(row,Numact)
            Conred(row,Numact) = sav
         ENDIF
      ENDDO
   ENDIF
! ======================================================================
! UPDATE LIST OF ACTIVE CONSTRAINTS
! ----------------------------------------------------------------------
   Actcon(Numact) = Inc
   Conact(Inc) = Numact
! ======================================================================
! REDUCE FOR NEW ACTIVE CONSTRAINT
! ----------------------------------------------------------------------
   IF ( dabs(Conred(Inc,Numact))<1D-12 ) THEN
      WRITE (str,*) "OGINCL-WARNING: CONSTRAINT SINGULAR"
      CALL ogwrit(2,str)
      WRITE (str,*) "INC=" , Inc
      CALL ogwrit(2,str)
      WRITE (str,*) "PIV=" , Conred(Inc,Numact)
      CALL ogwrit(2,str)
      Numact = Numact - 1
      Conact(Inc) = 0
      RETURN
   ENDIF
! ----------------------------------------------------------------------
   val = dsqrt(sum(Conred(Inc,Numact:lst)**2))
   IF ( Conred(Inc,Numact)>0D0 ) val = -val
! ----------------------------------------------------------------------
   Conred(Inc,Numact) = Conred(Inc,Numact) - val
! ----------------------------------------------------------------------
   sav = Conred(Inc,Numact)
   fac = 1D0/sav
   Conred(Inc,Numact:lst) = Conred(Inc,Numact:lst)*fac
! ----------------------------------------------------------------------
   fac = sav/val
   DO row = 1 , Numcon + 3
      IF ( Conact(row)<=0 ) THEN
         gam = dot_product(Conred(row,Numact:lst),Conred(Inc,Numact:lst))
         IF ( gam/=0D0 ) THEN
            gam = gam*fac
            Conred(row,Numact:lst) = Conred(row,Numact:lst) + Conred(Inc,Numact:lst)*gam
         ENDIF
      ENDIF
   ENDDO
! ----------------------------------------------------------------------
   Conred(Inc,Numact) = val
   Conred(Inc,Numact+1:lst) = 0D0
! ======================================================================
END SUBROUTINE ogincl
