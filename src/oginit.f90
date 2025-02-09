SUBROUTINE oginit(Varnum,Connum)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! ALLOCATION OF ARRAYS AND INITIALISATION OF PARAMETERS
! ======================================================================
! INP | VARNUM           | I*4 | NUMBER OF VARIABLES
! ----------------------------------------------------------------------
! INP | CONNUM           | I*4 | NUMBER OF CONSTRAINTS
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   INTEGER(4) Varnum
   INTEGER(4) Connum
! ======================================================================
   INTEGER(4) var , con
! ======================================================================
! VARIABLES
! ----------------------------------------------------------------------
   Numvar = Varnum
! ----------------------------------------------------------------------
   ALLOCATE (Varval(Numvar))
   ALLOCATE (Vartyp(Numvar))
   ALLOCATE (Varsca(Numvar))
   ALLOCATE (Varstr(Numvar))
   ALLOCATE (Varlen(Numvar))
   ALLOCATE (Varref(Numvar))
   ALLOCATE (Vardes(Numvar))
   ALLOCATE (Vargrd(Numvar))
   ALLOCATE (Vardir(Numvar))
   ALLOCATE (Funvar(Numvar))
   ALLOCATE (Senvar(Numvar))
! ----------------------------------------------------------------------
   DO var = 1 , Numvar
      Varval(var) = 0D0
      Vartyp(var) = 0
      Varsca(var) = 1D0
      Varstr(var) = ""
      Varlen(var) = 0
      Varref(var) = 0D0
      Vardes(var) = 0D0
      Vargrd(var) = 0D0
      Vardir(var) = 0D0
      Funvar(var) = 0D0
      Senvar(var) = 0D0
   ENDDO
! ======================================================================
! CONSTRAINTS
! ----------------------------------------------------------------------
   Numcon = Connum
! ----------------------------------------------------------------------
   ALLOCATE (Conval(Numcon+1))
   ALLOCATE (Contyp(Numcon+1))
   ALLOCATE (Conpri(Numcon+1))
   ALLOCATE (Consca(Numcon+1))
   ALLOCATE (Constr(Numcon+1))
   ALLOCATE (Conlen(Numcon+1))
   ALLOCATE (Conref(Numcon+1))
   ALLOCATE (Senqua(Numcon+1))
   ALLOCATE (Sencon(Numcon+1))
   ALLOCATE (Sendel(Numcon+1))
   ALLOCATE (Senact(Numcon+1))
! ----------------------------------------------------------------------
   DO con = 1 , Numcon + 1
      Conval(con) = 0D0
      Contyp(con) = 0
      Conpri(con) = 1
      Consca(con) = 1D0
      Constr(con) = ""
      Conlen(con) = 0
      Conref(con) = 0D0
      Senqua(con) = 0D0
      Sencon(con) = 0D0
      Sendel(con) = 0D0
      Senact(con) = 0
   ENDDO
! ======================================================================
! CONTROL
! ----------------------------------------------------------------------
   Optmet = 2
   Maxite = 10
   Corite = 10
   Optite = 10
   Divite = 10
   Cnvite = 10
   Varmax = 10D0
   Varsnd = 1D0
! ======================================================================
! DERIVATIVES
! ----------------------------------------------------------------------
   Varder = 1
! ----------------------------------------------------------------------
   ALLOCATE (Varper(Numvar))
! ----------------------------------------------------------------------
   DO var = 1 , Numvar
      Varper(var) = 1D-03
   ENDDO
! ======================================================================
! LOG FILE
! ----------------------------------------------------------------------
   Loglun = 6
   Loglev = 1
! ======================================================================
! PYGMO LOG FILE
! ----------------------------------------------------------------------
   Loglup = 7
   Loglev = 0
! ======================================================================
! MATLAB CONSOLE
! ----------------------------------------------------------------------
   Matlev = 0
! ======================================================================
! TABLE FILE
! ----------------------------------------------------------------------
   Tablun = 6
   Tablev = 0
! ======================================================================
! LINEAR OPTIMISATION MODE
! ----------------------------------------------------------------------
   Senopt = 0
! ======================================================================
! WORKING VECTORS
! ----------------------------------------------------------------------
   ALLOCATE (Actcon(Numcon+1))
   ALLOCATE (Confix(Numcon))
   ALLOCATE (Conact(Numcon+4))
   ALLOCATE (Conder(Numcon+3,Numvar))
   ALLOCATE (Conred(Numcon+3,Numvar))
   ALLOCATE (Sender(Numcon+3,Numvar))
   ALLOCATE (Conopt(Numcon+1))
! ----------------------------------------------------------------------
   Numact = 0
   Actcon = 0
   Conact = 0
   Confix = 0
   Conder = 0D0
   Conred = 0D0
   Conopt = 0
! ======================================================================
END SUBROUTINE oginit
