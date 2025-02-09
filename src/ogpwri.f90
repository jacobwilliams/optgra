SUBROUTINE ogpwri(Objval,Numvio,Convio)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! WRITE OPTIMIZATION LOG IN PYGMO FORMAT
! ======================================================================
! INP | OBJVAL           | R*8 | OBJECTIVE VALUE
! ----------------------------------------------------------------------
! INP | NUMVIO           | I*4 | NUMBER OF VIOLATED CONSTRAINTS
! ----------------------------------------------------------------------
! INP | CONVIO           | R*8 | TOTAL CONSTRAINT VIOLATION
! ======================================================================
! 2023/01/25 | W. MARTENS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   CHARACTER feas*2 , fmt*24
   REAL(8) Objval , Convio
   INTEGER(4) Numvio
! ======================================================================
   IF ( Verbos==0 ) RETURN
! Print header
   IF ( Fevals==0 ) CALL ogpwri_start()
! Increase counter for cost function evaluations
   Fevals = Fevals + 1
! Every 50 lines print the column names.
   IF ( mod(real(Fevals-1D0)/real(Verbos),50D0)==0D0 ) WRITE (Loglup,'(A10,A15,A15,A15,A2)') "objevals:" , "objval:" ,             &
      & "violated:" , "viol. norm:"
   IF ( Verbos/=0 .AND. mod(Fevals,Verbos)==0D0 ) THEN
      IF ( Convio>0D0 ) THEN
         feas = " i"
      ELSE
         feas = "  "
      ENDIF

! Write the log line (different format depending on violation size)
      IF ( Convio==0D0 ) THEN
         fmt = '(I10,F15.4,I15,I15,A2)'
         WRITE (Loglup,fmt) Fevals , Objval , Numvio , int(Convio) , feas
      ELSEIF ( Convio>1D-3 ) THEN
         fmt = '(I10,F15.4,I15,F15.6,A2)'
         WRITE (Loglup,fmt) Fevals , Objval , Numvio , Convio , feas
      ELSE
         fmt = '(I10,F15.4,I15,E15.6,A2)'
         WRITE (Loglup,fmt) Fevals , Objval , Numvio , Convio , feas
      ENDIF
   ENDIF

! Write final summary
   IF ( Pygfla/=0 ) CALL ogpwri_end(Objval,Numvio,Convio)
! ======================================================================
END SUBROUTINE ogpwri
