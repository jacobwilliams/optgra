SUBROUTINE ogeval(Valvar,Valcon,Dervar,Dercon,calval,calder)
! ======================================================================
! COMPUTES SCALED CONTRAINTS+MERIT AND DERIVATIVES
! FROM     SCALED VARIABLES
! ======================================================================
! INP | DERVAR           | I*4 | DERIVATIVES COMPUTATION MODE
!     |                  |     | -> 0: VALUES ONLY
!     |                  |     | -> 1: USER DEFINED
!     |                  |     | -> 2: NUMERIC WITH DOUBLE DIFFERENCING
!     |                  |     | -> 3: NUMERIC WITH SINGLE DIFFERENCING
! ----------------------------------------------------------------------
! INP | CALVAL           | EXT | FUNCTION FOR VALUES
!     |                  |     | -> CALVAL (VALVAR, VALCON)
!     |                  |     | -> INPUT AND OUTPUT NOT SCALED
! ----------------------------------------------------------------------
! INP | CALDER           | EXT | FUNCTION FOR VALUES AND DERIVATIVES
!     |                  |     | -> CALDER (VALVAR, VALCON, DERCON)
!     |                  |     | -> INPUT AND OUTPUT NOT SCALED
! ======================================================================
! SUBROUTINES CALLED: CALVAL, CALDER
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   REAL(8) Valvar(Numvar)
   REAL(8) Valcon(Numcon+1)
   INTEGER(4) Dervar
   REAL(8) Dercon(Numcon+1,Numvar)
   EXTERNAL calval
   EXTERNAL calder
! ======================================================================
   INTEGER(4) var , con , cod , len , ind , numvio
   REAL(8) val , sca , fac , per , sav , der , err , conerr , convio
   CHARACTER typ*3 , sta*3 , nam*80 , str*256
   REAL(8) ggg(4,4) , bbb(4) , vvv(4) , objval
! ======================================================================
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varvec
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: convec
! ----------------------------------------------------------------------
   ALLOCATE (varvec(Numvar))
   ALLOCATE (convec(Numcon+1))
! ======================================================================
! GENERAL
! ----------------------------------------------------------------------
   WRITE (str,'()')
   CALL ogwrit(3,str)
   IF ( Dervar==0 ) THEN
      WRITE (str,'("COMPUTE RESULTS")')
      CALL ogwrit(3,str)
   ELSEIF ( Dervar==1 .OR. Dervar==-1 ) THEN
      WRITE (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES USER DEFINED")')
      CALL ogwrit(3,str)
   ELSEIF ( Dervar==2 ) THEN
      WRITE (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES BY DOUBLE DIFFERENCING")')
      CALL ogwrit(3,str)
   ELSEIF ( Dervar==3 ) THEN
      WRITE (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES BY SINGLE DIFFERENCING")')
      CALL ogwrit(3,str)
   ENDIF
! ======================================================================
! WRITE VARIABLES
! ----------------------------------------------------------------------
   WRITE (str,'()')
   CALL ogwrit(3,str)
   WRITE (str,'("VARIABLES NOT SCALED:")')
   CALL ogwrit(3,str)
   WRITE (str,'()')
   CALL ogwrit(3,str)
   DO var = 1 , Numvar
      val = Valvar(var)
      sca = Varsca(var)
      cod = Vartyp(var)
      IF ( cod==0 ) typ = "FRE"
      IF ( cod==1 ) typ = "PAR"
      nam = Varstr(var)
      len = Varlen(var)
      WRITE (str,'("VAR/VAL/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') var , val*sca , sca , typ , nam(1:len)
      CALL ogwrit(3,str)
   ENDDO
! ======================================================================
! DE-SCALE VARIABLES
! ----------------------------------------------------------------------
   DO var = 1 , Numvar
      varvec(var) = Valvar(var)*Varsca(var)
   ENDDO
! ======================================================================
! GET RESULTS
! GET DERIVATIVES IF USER DEFINED
! ----------------------------------------------------------------------
   IF ( Dervar==0 ) THEN
      CALL calval(varvec,Valcon,0)
   ELSEIF ( Dervar==1 .OR. Dervar==-1 ) THEN
      CALL calval(varvec,Valcon,1)
   ELSEIF ( Dervar==2 ) THEN
      CALL calval(varvec,Valcon,1)
   ELSEIF ( Dervar==3 ) THEN
      CALL calval(varvec,Valcon,1)
   ENDIF
! ======================================================================
   IF ( 1==2 ) THEN
      ggg(1,1) = +1D+01
      ggg(2,1) = +1D+00
      ggg(3,1) = +2D+00
      ggg(4,1) = +3D+00
      ggg(1,2) = ggg(2,1)
      ggg(2,2) = +1D+01
      ggg(3,2) = +4D+00
      ggg(4,2) = +5D+00
      ggg(1,3) = ggg(3,1)
      ggg(2,3) = ggg(3,2)
      ggg(3,3) = +1D+01
      ggg(4,3) = +6D+00
      ggg(1,4) = ggg(4,1)
      ggg(2,4) = ggg(4,2)
      ggg(3,4) = ggg(4,3)
      ggg(4,4) = +1D+01
! ----------------------------------------------------------------------
      bbb(1) = +1D+01
      bbb(2) = +1D+01
      bbb(3) = +1D+01
      bbb(4) = +1D+01
! ----------------------------------------------------------------------
      CALL mul2m(ggg,4,1,1,4,Valvar,4,1,1,4,vvv,4,1,1,1)
      CALL mul2m(Valvar,1,1,1,1,vvv,4,1,1,4,Valcon,1,1,1,1)
      CALL mul2m(bbb,1,1,1,1,Valvar,4,1,1,4,vvv,1,1,1,1)
      CALL mulvs(Valcon,0.5D0,Valcon,1)
      CALL sum2v(Valcon,vvv,Valcon,1)
      CALL mulvs(Valcon,-1D0,Valcon,1)
      WRITE (str,*) "VALCON=" , (Valcon(ind),ind=1,1)
      CALL ogwrit(3,str)
   ENDIF
! ----------------------------------------------------------------------
! SCALE RESULTS
! ----------------------------------------------------------------------
   DO con = 1 , Numcon + 1
      convec(con) = Valcon(con)
      sca = Consca(con)
      cod = Contyp(con)
      IF ( cod==-1 ) sca = -sca
      Valcon(con) = Valcon(con)/sca
   ENDDO
! ======================================================================
! WRITE RESULTS
! ----------------------------------------------------------------------
   WRITE (str,'()')
   CALL ogwrit(3,str)
   WRITE (str,'("RESULTS NOT SCALED:")')
   CALL ogwrit(3,str)
   WRITE (str,'()')
   CALL ogwrit(3,str)
   conerr = 0D0     ! total constraint error (scaled to constr. threshod)
   convio = 0D0     ! total constaint error norm (unscaled)
   ind = 0     ! index of largest constraint violation
   fac = 0D0     ! value of largest constraint violation
   numvio = 0     ! number of violated constraints
   DO con = 1 , Numcon + 1
      val = Valcon(con)
      sca = Consca(con)
      cod = Contyp(con)
      sta = "   "
      err = 0D0
      IF ( cod==-1 ) sca = -sca
      IF ( cod==-2 ) typ = "DER"
      IF ( cod==-1 .AND. con<=Numcon ) typ = "LTE"
      IF ( cod==0 .AND. con<=Numcon ) typ = "EQU"
      IF ( cod==1 .AND. con<=Numcon ) typ = "GTE"
      IF ( cod==-1 .AND. con>Numcon ) typ = "MIN"
      IF ( cod==1 .AND. con>Numcon ) typ = "MAX"
      IF ( cod==0 .AND. con<=Numcon .AND. dabs(val)>1D0 ) THEN
         sta = "VIO"
         err = dabs(val)
         numvio = numvio + 1
      ENDIF
      IF ( cod/=0 .AND. con<=Numcon .AND. cod/=-2 .AND. -val>1D0 ) THEN
         sta = "VIO"
         err = dabs(val)
         numvio = numvio + 1
      ENDIF
      conerr = conerr + err
      convio = convio + (err*sca)**2
      IF ( err>fac ) ind = con
      IF ( err>fac ) fac = err
      nam = Constr(con)
      len = Conlen(con)
      WRITE (str,'("CON/VAL/SCA/TYP/STA/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A3,1X,A)') con , val*sca , sca , typ , sta , nam(1:len)
      CALL ogwrit(3,str)
   ENDDO
   WRITE (str,'()')
   CALL ogwrit(3,str)
   WRITE (str,'("CONSTRAINT ERROR.:",2(1X,D13.6),I6)') conerr , fac , ind
   CALL ogwrit(3,str)
   WRITE (str,'()')
   CALL ogwrit(3,str)
! write pygmo-style log output
   objval = -Valcon(Numcon+1)
   convio = dsqrt(convio)
   CALL ogpwri(objval,numvio,convio,Dervar)
! ======================================================================
! NO DERIVATIVES
! ----------------------------------------------------------------------
   IF ( Dervar==0 ) THEN
      RETURN
   ELSEIF ( Dervar==1 .OR. Dervar==-1 ) THEN
      CALL calder(varvec,convec,Dercon)
   ENDIF
! ----------------------------------------------------------------------
   IF ( 1==2 ) THEN
      CALL mul2m(Valvar,1,1,1,1,ggg,4,1,1,4,Dercon,1,1,1,4)
      CALL sum2v(Dercon,bbb,Dercon,4)
      CALL mulvs(Dercon,-1D0,Dercon,4)
      WRITE (str,*) "DERCON=" , (Dercon(1,ind),ind=1,4)
      CALL ogwrit(3,str)
   ENDIF
! ======================================================================
! WRITE DERIVATIVES
! ----------------------------------------------------------------------
   WRITE (str,'()')
   CALL ogwrit(3,str)
   WRITE (str,'("DERIVATIVES SCALED:")')
   CALL ogwrit(3,str)
   WRITE (str,'()')
   CALL ogwrit(3,str)
! ----------------------------------------------------------------------
   DO var = 1 , Numvar
! ----------------------------------------------------------------------
! WRITE VARIABLE
! ----------------------------------------------------------------------
      val = Valvar(var)
      sca = Varsca(var)
      cod = Vartyp(var)
      IF ( cod==0 ) typ = "FRE"
      IF ( cod==1 ) typ = "PAR"
      nam = Varstr(var)
      len = Varlen(var)
      WRITE (str,'("VAR/VAL/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') var , val*sca , sca , typ , nam(1:len)
      CALL ogwrit(4,str)
      WRITE (str,'()')
      CALL ogwrit(4,str)
! ----------------------------------------------------------------------
! DERIVATIVES BY DOUBLE DIFFERENCING
! ----------------------------------------------------------------------
      IF ( Dervar==2 ) THEN
! ----------------------------------------------------------------------
         per = Varper(var)
         sav = varvec(var)
         varvec(var) = sav + per
         CALL calval(varvec,Dercon(1,var),0)
         varvec(var) = sav - per
         CALL calval(varvec,convec,0)
         fac = 0.5D0/per
         DO con = 1 , Numcon + 1
            Dercon(con,var) = (Dercon(con,var)-convec(con))*fac
         ENDDO
         varvec(var) = sav
! ----------------------------------------------------------------------
      ENDIF
! ----------------------------------------------------------------------
! DERIVATIVES BY SINGLE DIFFERENCING
! ----------------------------------------------------------------------
      IF ( Dervar==3 ) THEN
! ----------------------------------------------------------------------
         per = Varper(var)
         sav = varvec(var)
         varvec(var) = sav + per
         CALL calval(varvec,Dercon(1,var),0)
         fac = 1.0D0/per
         DO con = 1 , Numcon + 1
            Dercon(con,var) = (Dercon(con,var)-convec(con))*fac
         ENDDO
         varvec(var) = sav
! ----------------------------------------------------------------------
      ENDIF
! ----------------------------------------------------------------------
! SCALE DERIVATIVES
! ----------------------------------------------------------------------
      DO con = 1 , Numcon + 1
         fac = Varsca(var)/Consca(con)
         IF ( Contyp(con)==-1 ) fac = -fac
         Dercon(con,var) = Dercon(con,var)*fac
      ENDDO
! ----------------------------------------------------------------------
! WRITE DERIVATIVES
! ======================================================================
      DO con = 1 , Numcon + 1
         der = Dercon(con,var)
         IF ( der/=0D0 ) THEN
            sca = Consca(con)
            cod = Contyp(con)
            IF ( cod==-1 ) sca = -sca
            IF ( cod==-2 ) typ = "DER"
            IF ( cod==-1 .AND. con<=Numcon ) typ = "LTE"
            IF ( cod==0 .AND. con<=Numcon ) typ = "EQU"
            IF ( cod==1 .AND. con<=Numcon ) typ = "GTE"
            IF ( cod==-1 .AND. con>Numcon ) typ = "MIN"
            IF ( cod==1 .AND. con>Numcon ) typ = "MAX"
            nam = Constr(con)
            len = Conlen(con)
            WRITE (str,'("CON/DER/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') con , der*sca/Varsca(var) , sca , typ ,          &
                 & nam(1:len)
            CALL ogwrit(4,str)
         ENDIF
      ENDDO
      WRITE (str,'()')
      CALL ogwrit(4,str)
! ----------------------------------------------------------------------
   ENDDO
! ======================================================================
   DEALLOCATE (varvec)
   DEALLOCATE (convec)
! ======================================================================
END SUBROUTINE ogeval
