SUBROUTINE ogsens(Consta,Concon,Convar,Varcon,Varvar)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
! ======================================================================
! OUT | CONSTA(NUMCON)   | I*4 | CONSTRAINT STATUS (0=PAS 1=ACT)
! OUT | CONCON(NUMCON+1, | R*8 | SENSITIVITY OF CONTRAINTS+MERIT W.R.T.
!     |        NUMCON)   |     |                ACTIVE CONSTRAINTS
! OUT | CONVAR(NUMCON+1, | R*8 | SENSITIVITY OF CONTRAINTS+MERIT W.R.T.
!     |        NUMVAR)   |     |                PARAMETERS
! OUT | VARCON(NUMVAR  , | R*8 | SENSITIVITY OF VARIABLES W.R.T.
!     |        NUMCON)   |     |                ACTIVE CONSTRAINTS
! OUT | VARVAR(NUMVAR  , | R*8 | SENSITIVITY OF VARIABLES W.R.T.
!     |        NUMVAR)   |     |                PARAMETERS
!     |                  |     | -> NOT SCALED
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   INTEGER(4) Consta(Numcon)
   REAL(8) Concon(Numcon+1,Numcon)
   REAL(8) Convar(Numcon+1,Numvar)
   REAL(8) Varcon(Numvar,Numcon)
   REAL(8) Varvar(Numvar,Numvar)
! ======================================================================
   REAL(8) val , sca
   INTEGER(4) var , con , act , par , ind , typ
! ======================================================================
! CONVERGED
! ----------------------------------------------------------------------
   Consta = 0
   DO act = 1 , Numact
      con = Actcon(act)
      Consta(con) = 1
   ENDDO
! ======================================================================
! SENSITIVITY OF CONTRAINTS W.R.T. ACTIVE CONSTRAINTS
! ----------------------------------------------------------------------
   Concon = 0D0
   DO con = 1 , Numcon + 1
      IF ( Conact(con)>0 ) Concon(con,con) = 1D0
      IF ( Conact(con)>0 ) CYCLE
      Conref = Conred(con,1:Numact)
      CALL ogrigt(Conref,Conref)
      DO act = 1 , Numact
         ind = Actcon(act)
         Concon(con,ind) = -Conref(act)
      ENDDO
   ENDDO
! ======================================================================
! SENSITIVITY OF CONSTRAINTS W.R.T. PARAMETERS
! ----------------------------------------------------------------------
   Convar = 0D0
   DO con = 1 , Numcon + 1
      IF ( Conact(con)>0 ) CYCLE
      DO var = 1 , Numvar
         IF ( Vartyp(var)==0 ) CYCLE
         val = Sender(con,var)
         DO act = 1 , Numact
            ind = Actcon(act)
            val = val + Concon(con,ind)*Sender(ind,var)
         ENDDO
         Convar(con,var) = val
      ENDDO
   ENDDO
! ======================================================================
! SENSITIVITY OF VARIABLES W.R.T. ACTIVE CONSTRAINTS
! ----------------------------------------------------------------------
   Varcon = 0D0
   DO var = 1 , Numvar
      IF ( Vartyp(var)/=0 ) CYCLE
      DO act = 1 , Numact
         con = Actcon(act)
         Conref(act) = Conder(con,var)
      ENDDO
      CALL ogleft(Conref,Conref)
      CALL ogrigt(Conref,Conref)
      DO act = 1 , Numact
         con = Actcon(act)
         Varcon(var,con) = -Conref(act)
      ENDDO
   ENDDO
! ======================================================================
! SENSITIVITY OF VARIABLES W.R.T. PARAMETERS
! ----------------------------------------------------------------------
   Varvar = 0D0
   DO par = 1 , Numvar
      Varvar(par,par) = 1D0
      IF ( Vartyp(par)/=1 ) CYCLE
      DO var = 1 , Numvar
         IF ( Vartyp(var)/=0 ) CYCLE
         val = 0D0
         DO act = 1 , Numact
            con = Actcon(act)
            val = val + Varcon(var,con)*Sender(con,par)
         ENDDO
         Varvar(var,par) = val
      ENDDO
   ENDDO
! ======================================================================
! DESCALE SENSITIVITY
! ----------------------------------------------------------------------
   DO con = 1 , Numcon + 1
      typ = Contyp(con)
      sca = Consca(con)
      IF ( typ<0 ) sca = -sca
      Convar(con,1:Numvar) = Convar(con,1:Numvar)*sca
      Concon(con,1:Numcon) = Concon(con,1:Numcon)*sca
      IF ( con>Numcon ) CYCLE
      Varcon(1:Numvar,con) = Varcon(1:Numvar,con)/sca
      Concon(1:Numcon+1,con) = Concon(1:Numcon+1,con)/sca
   ENDDO
! ----------------------------------------------------------------------
   DO var = 1 , Numvar
      sca = Varsca(var)
      Varcon(var,1:Numcon) = Varcon(var,1:Numcon)*sca
      Varvar(var,1:Numvar) = Varvar(var,1:Numvar)*sca
      Convar(1:Numcon+1,var) = Convar(1:Numcon+1,var)/sca
      Varvar(1:Numvar,var) = Varvar(1:Numvar,var)/sca
   ENDDO
! ======================================================================
END SUBROUTINE ogsens
