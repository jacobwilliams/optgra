
SUBROUTINE ogpwri_end(Objval,Numvio,Convio)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! WRITE OPTIMIZATION END RESULT IN PYGMO FORMAT
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
   REAL(8) Objval , Convio
   INTEGER(4) Numvio
! ======================================================================
   IF ( Pygfla==0 ) RETURN
! Write termination message
   WRITE (Loglup,'("")')
   WRITE (Loglup,'("Final values after iteration        ", I10:)') Numite
   WRITE (Loglup,'("Final objective value:              ", F10.4)') Objval
   WRITE (Loglup,'("Final constraint violation:         ", F10.4)') Convio
   WRITE (Loglup,'("Final num. of violated constraints: ",I10)') Numvio
   IF ( Pygfla==1 ) THEN
      WRITE (Loglup,'("Successful termination: Optimal solution found.")')
   ELSEIF ( Pygfla==2 ) THEN
      WRITE (Loglup,'("Successful termination: Constraints matched.")')
   ELSEIF ( Pygfla==3 ) THEN
      WRITE (Loglup,'("Not converged.")')
   ELSEIF ( Pygfla==4 ) THEN
      WRITE (Loglup,'("Problem appears infeasible.")')
   ENDIF
   WRITE (Loglup,'("")')

! ======================================================================
END SUBROUTINE ogpwri_end
