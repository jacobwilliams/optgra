SUBROUTINE ogexec(Valvar,Valcon,Finopt,Finite,Calval,Calder)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN
! ======================================================================
! I/O | VALVAR(NUMVAR)   | R*8 | VARIABLES VALUE
!     |                  |     | -> NOT SCALED
! ----------------------------------------------------------------------
! OUT | VALCON(NUMCON+1) | R*8 | CONSTRAINTS VALUE (1:NUMCON)
!     |                  |     | MERIT       VALUE (1+NUMCON)
!     |                  |     | -> NOT SCALED
! ----------------------------------------------------------------------
! OUT | FINOPT           | I*4 | TERMINATION STATUS
!     |                  |     | -> 1=    MATCHED &     OPTIMAL
!     |                  |     | -> 2=    MATCHED & NOT OPTIMAL
!     |                  |     | -> 3=NOT MATCHED & NOT OPTIMAL
!     |                  |     | -> 4=NOT FEASIBL & NOT OPTIMAL
! ----------------------------------------------------------------------
! INP | CALVAL           | EXT | FUNCTION FOR VALUES
!     |                  |     | -> CALDER (VALVAR, VALCON)
!     |                  |     | -> INPUT AND OUTPUT NOT SCALED
! ----------------------------------------------------------------------
! INP | CALDER           | EXT | FUNCTION FOR VALUES AND DERIVATIVES
!     |                  |     | -> CALDER (VALVAR, VALCON, CONDER)
!     |                  |     | -> INPUT AND OUTPUT NOT SCALED
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   REAL(8) Valvar(Numvar)
   REAL(8) Valcon(Numcon+1)
   INTEGER(4) Finopt , finish , Finite , itecor , iteopt
   EXTERNAL Calval
   EXTERNAL Calder
! ======================================================================
   INTEGER(4) var , con , typ , len , num , numvio
   REAL(8) val , sca , red , der , fac , old , convio
   CHARACTER str*256 , nam*256
! ----------------------------------------------------------------------
   INTEGER(4) numequ , itediv , itecnv
   REAL(8) varacc , cosnew , cosold , varsav , meamer
   REAL(8) conerr , desnor , norerr , meaerr
! ======================================================================
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varsum
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varcor
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: concor
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
! ----------------------------------------------------------------------
         ALLOCATE (varsum(Numvar))
         ALLOCATE (varcor(Numvar))
         ALLOCATE (concor(Numcon+1))
! ======================================================================
! GENERAL
! ----------------------------------------------------------------------
!     LOGLEV = 2
         CALL ogwrit(2,"")
         CALL ogwrit(2,"OPTGRA START")
         CALL ogwrit(2,"")
         Finopt = 3
         itecor = 0
         iteopt = 0
         meaerr = 0D0
         meamer = 0D0
         itediv = 0
         itecnv = 0
         Conopt = 0
         concor = 0D0
         varcor = 0D0
         desnor = 0D0
! ----------------------------------------------------------------------
         Varstp = Varsnd
         Numite = 0
         cosnew = 0D0
         DO var = 1 , Numvar
            varsum(var) = 0D0
            varcor(var) = 0D0
         ENDDO
! ======================================================================
! EQUALTIY CONSTRAINTS IN ACTIVE SET
! ----------------------------------------------------------------------
         Numact = 0
         DO con = 1 , Numcon
!          IF (CONTYP(CON) .EQ. -2) CYCLE
            nam = Constr(con)
            len = Conlen(con)
            WRITE (str,*) "CON/PRI=" , con , Conpri(con) , " " , nam(1:len)
            CALL ogwrit(3,str)
            Conact(con) = 0
            IF ( Consca(con)>=1D+09 ) Contyp(con) = -2
            IF ( Contyp(con)==0 ) THEN
            ELSEIF ( Contyp(con)==-2 ) THEN
               Conact(con) = -2
            ENDIF
         ENDDO
         numequ = Numact
         Conact(Numcon+1) = -3
         Conact(Numcon+2) = -3
! ======================================================================
! SCALE VARIABLES
! ----------------------------------------------------------------------
         DO var = 1 , Numvar
            Varval(var) = Valvar(var)/Varsca(var)
         ENDDO
! ======================================================================
! HEADER FOR TABLE
! ----------------------------------------------------------------------
         IF ( Tablev>=1 ) WRITE (Tablun,'("ITER",1X,"OPT",1X,1000(1X,I10))') (var,var=1,Numvar) , (con,con=1,Numcon)
         spag_nextblock_1 = 2
      CASE (2)
         SPAG_Loop_1_1: DO
! ======================================================================
!      IF (NUMITE .GE. 52) MATLEV = 3
!      IF (NUMITE .GE. 55) MATLEV = 2
! ======================================================================
! NEW ITERATION
! ----------------------------------------------------------------------
            IF ( Numite>=Corite .AND. itecor==0 ) THEN
               Finopt = 3
               Finite = Numite
               CALL ogwrit(1,"")
               WRITE (str,'("OPTGRA: Converged: not ITERAT=",2I4,2D11.3)') Numite , Maxite , conerr , desnor
               CALL ogwrit(1,str)

! Final Pygmo output
! TODO: can this final fitness call be avoided (just for output)?
               Pygfla = 3
                      ! pygmo flag in COMMON: no covergence

               CALL ogeval(Varval,Conval,Varder,Conder(1:Numcon+1,:),Calval,Calder)
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Numite>=Maxite .OR. (Numite-itecor>=Optite-1 .AND. itecor/=0) ) THEN
               Finopt = 2
               Finite = iteopt
               Varval = varcor
               Conval = concor
               CALL ogwrit(1,"")
               WRITE (str,'("OPTGRA: Converged: mat ITERAT=",2I4,2D11.3)') Numite , Maxite , conerr , desnor
               CALL ogwrit(1,str)
               CALL ogpwri_end(Numite,-Valcon(Numcon+1),numvio,convio,1)

! Final Pygmo output
! TODO: can this final fitness call be avoided (just for output)?
               Pygfla = 2
                      ! pygmo flag in COMMON: constraints matched

               CALL ogeval(Varval,Conval,Varder,Conder(1:Numcon+1,:),Calval,Calder)
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
! ----------------------------------------------------------------------
            Numite = Numite + 1
! ----------------------------------------------------------------------
            CALL ogwrit(3,"")
            WRITE (str,'("ITERAT=",I5)') Numite
            CALL ogwrit(2,str)
! ======================================================================
! GET VALUES AND GRADIENTS
! ======================================================================
            IF ( Senopt<=0 ) THEN
               CALL ogeval(Varval,Conval,Varder,Conder(1:Numcon+1,:),Calval,Calder)
            ELSEIF ( Senopt==+1 .OR. Senopt==+3 ) THEN
               Varval = Senvar
               CALL ogeval(Varval,Conval,0,Conder(1:Numcon+1,:),Calval,Calder)
            ELSEIF ( Senopt==+2 .OR. Senopt==+4 ) THEN
               Varval = Senvar
               DO con = 1 , Numcon + 1
                  sca = Consca(con)
                  IF ( Contyp(con)==-1 ) sca = -sca
                  Conval(con) = Sencon(con) - Sendel(con)/sca
               ENDDO
            ENDIF
            IF ( Senopt==-1 ) THEN
               Senvar = Varval
               Sencon = Conval
            ENDIF
! ======================================================================
            IF ( Varder==-1 .AND. Senopt<=0 ) THEN
               Conred(1:Numcon+1,:) = Conder(1:Numcon+1,:)
               CALL ogeval(Varval,Conval,2,Conder(1:Numcon+1,:),Calval,Calder)
               WRITE (str,'("GRADIENT CHECK")')
               CALL ogwrit(1,str)
               DO var = 1 , Numvar
                  DO con = 1 , Numcon + 1
                     fac = Varsca(var)/Consca(con)
                     fac = 1D0
                     der = Conder(con,var)*fac
                     red = Conred(con,var)*fac
                     IF ( abs(der)<1D-6 .AND. abs(red)<1D-6 ) CYCLE
                     IF ( dabs(der-red)<1D-2 ) CYCLE
                     IF ( der/=0D0 ) THEN
                        fac = red/der
                     ELSE
                        fac = 0D0
                     ENDIF
                     IF ( dabs(fac-1D0)<1D-2 ) CYCLE
                     WRITE (str,'("VAR/CON/ANA/NUM/A2N=",2I4,3(1X,D13.6))') var , con , red , der , fac
                     CALL ogwrit(1,str)
                     nam = Varstr(var)
                     len = Varlen(var)
                     WRITE (str,'("      VAR=",A,1X,D13.6)') nam(1:len) , Varval(var)*Varsca(var)
                     CALL ogwrit(1,str)
                     nam = Constr(con)
                     len = Conlen(con)
                     WRITE (str,'("      CON=",A,1X,D13.6)') nam(1:len) , Conval(con)*Consca(con)
                     CALL ogwrit(1,str)
                  ENDDO
               ENDDO
!          CONDER(1:NUMCON+1,:) = CONRED(1:NUMCON+1,:)
!          GOTO 9999
            ENDIF
! ======================================================================
            Sender = Conder
            DO var = 1 , Numvar
               IF ( Vartyp(var)/=1 ) CYCLE
!          WRITE (STR,*) "VAR=",VAR,VARVAL(VAR)*VARSCA(VAR)
!          CALL OGWRIT (2,STR)
               Conder(1:Numcon+1,var) = 0D0
            ENDDO
! ======================================================================
            IF ( Numite==1 ) THEN
               Vargrd = Varval
            ELSE
               Vargrd = Varref
            ENDIF
            Varref = Varval
            Conref = Conval
! ======================================================================
            varacc = 0D0
! ======================================================================
            cosold = cosnew
            cosnew = Conval(Numcon+1)
            CALL ogwrit(3,"")
            WRITE (str,'("OPTGRA: VALCOS=",D15.8,1X,D15.8)') cosnew , cosnew - cosold
            CALL ogwrit(3,str)
! ======================================================================
! CORRECTION PART
! ----------------------------------------------------------------------
            CALL ogcorr(varacc,finish,conerr,norerr,Calval,Calder)
! ----------------------------------------------------------------------
            IF ( Tablev>=1 ) WRITE (Tablun,'(I4,1X,"COR",1X,1000(1X,D10.3))') Numite , (Varval(var),var=1,Numvar) ,                &
                                  & (Conval(con),con=1,Numcon)
! ----------------------------------------------------------------------
            IF ( Senopt/=0 ) THEN
               IF ( finish/=0 ) EXIT SPAG_Loop_1_1
               Finopt = 0
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
! ----------------------------------------------------------------------
            IF ( finish==0 ) THEN
               Numact = 0
               old = meaerr
               itediv = itediv + 1
               num = min(itediv,Divite)
               meaerr = (meaerr*(num-1)+norerr)/num
!         WRITE (STR,*) "MEAERR=",MEAERR
!         CALL OGWRIT (2,STR)
               IF ( itediv<Divite .OR. meaerr<=old ) CYCLE
               finish = -1
            ENDIF
! ----------------------------------------------------------------------
            IF ( finish==-1 ) THEN
               Finopt = 4
               Finite = Numite
               CALL ogwrit(1,"")
               WRITE (str,'("OPTGRA: Converged: unf ITERAT=",2I4,2D11.3)') Numite , Maxite , conerr , desnor
               CALL ogwrit(1,str)

! Final Pygmo output
! TODO: can this final fitness call be avoided (just for output)?
               Pygfla = 4
                      ! pygmo flag in COMMON: infeasible

               CALL ogeval(Varval,Conval,Varder,Conder(1:Numcon+1,:),Calval,Calder)
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
! ----------------------------------------------------------------------
            itediv = 0
            iteopt = Numite
            IF ( itecor==0 .OR. concor(Numcon+1)<Conval(Numcon+1) ) THEN
               varcor = Varval
               concor = Conval
            ENDIF
            IF ( itecor==0 ) itecor = Numite
! ----------------------------------------------------------------------
            old = meamer
            itecnv = itecnv + 1
            num = min(itecnv,Cnvite)
            meamer = (meamer*(num-1)+concor(Numcon+1))/num
!      WRITE (STR,*) "MEAMER=",ITECNV,NUM,MEAMER,OLD,OLD/MEAMER
!      CALL OGWRIT (-1,STR)
            IF ( itecnv>=Cnvite .AND. meamer<old ) THEN
               Finopt = 2
               Finite = iteopt
               Varval = varcor
               Conval = concor
               CALL ogwrit(1,"")
               WRITE (str,'("OPTGRA: Converged: mat ITERAT=",2I4,2D11.3)') Numite , Maxite , conerr , desnor
               CALL ogwrit(1,str)

! Final Pygmo output
! TODO: can this final fitness call be avoided (just for output)?
               Pygfla = 2
                      ! pygmo flag in COMMON: matched

               CALL ogeval(Varval,Conval,Varder,Conder(1:Numcon+1,:),Calval,Calder)
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            EXIT SPAG_Loop_1_1
         ENDDO SPAG_Loop_1_1
! ======================================================================
! OPTIMISATION PART
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
         IF ( Senopt<+3 ) THEN
            varsav = Varmax
            Varmax = Varmax*10D-1
            CALL ogopti(varacc,numequ,finish,desnor,Calval)
            Varmax = varsav
         ENDIF
! ----------------------------------------------------------------------
         IF ( Senopt/=0 ) THEN
            CALL ogwrit(1,"")
            IF ( finish==0 ) THEN
               Finopt = 0
               CALL ogwrit(1,"OPTGRA sensitivity converged: not")
            ELSE
               Finopt = 1
               CALL ogwrit(1,"OPTGRA sensitivity converged: yes")
            ENDIF
            spag_nextblock_1 = 3
            CYCLE SPAG_DispatchLoop_1
         ENDIF
! ----------------------------------------------------------------------
         IF ( finish==0 ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
! ======================================================================
! NOT CONVERGED
! ----------------------------------------------------------------------
         IF ( varacc/=0D0 ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
! ======================================================================
! CONVERGED
! ----------------------------------------------------------------------
         Finopt = 1
         Finite = Numite
         CALL ogwrit(1,"")
         WRITE (str,'("OPTGRA: Converged: yes ITERAT=",2I4,2D11.3)') Numite , Maxite , conerr , desnor
         CALL ogwrit(1,str)
         CALL ogwrit(3,"")

! Final Pygmo output
! TODO: can this final fitness call be avoided (just for output)?
         Pygfla = 1
                  ! covergence
         CALL ogeval(Varval,Conval,Varder,Conder(1:Numcon+1,:),Calval,Calder)
         spag_nextblock_1 = 3
      CASE (3)

!      WRITE (STR,*) "DIF=",NORM2(VARVAL-VARREF)
!      CALL OGWRIT (1,STR)
! ======================================================================
! ======================================================================
! DESCALE VARIABLES
! ----------------------------------------------------------------------
!      CALL OGWMAT (3)
!      CALL OGEVAL (VARVAL, VALCON, 0, CONDER, CALVAL, CALDER)
         DO var = 1 , Numvar
            Valvar(var) = Varval(var)*Varsca(var)
         ENDDO
!      IF (SENOPT .NE. 0) THEN
!          CALL OGEVAL (VARVAL, VALCON, 0, CONDER, CALVAL, CALVAL)
!      ENDIF
! ======================================================================
! DESCALE VALUES
! ----------------------------------------------------------------------
         DO con = 1 , Numcon + 1
            typ = Contyp(con)
            sca = Consca(con)
            IF ( typ==-1 ) sca = -sca
            Valcon(con) = Conval(con)*sca
         ENDDO
! ======================================================================
         CALL ogwrit(3,"")
         CALL ogwrit(3,"STATUS OF CONSTRAINTS:")
         CALL ogwrit(3,"")
         CALL ogwrit(3," ACT  PAS  NON COST___VAL CONSTRAINT")
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            nam = Constr(con)
            len = Conlen(con)
            val = Conval(con)
            IF ( Conact(con)>0 ) THEN
               WRITE (str,'( I4,5X,6X,D10.3,1X,A)') con , val , nam(1:len)
               CALL ogwrit(3,str)
            ELSEIF ( Conact(con)==0 ) THEN
               WRITE (str,'( 5X,I4,6X,D10.3,1X,A)') con , val , nam(1:len)
               CALL ogwrit(3,str)
            ELSEIF ( Conact(con)<0 ) THEN
               WRITE (str,'(10X,I4,1X,D10.3,1X,A)') con , val , nam(1:len)
               CALL ogwrit(3,str)
            ENDIF
         ENDDO

         DEALLOCATE (varsum)
         DEALLOCATE (varcor)
         DEALLOCATE (concor)
! ----------------------------------------------------------------------
         CALL ogwrit(2,"")
         CALL ogwrit(2,"OPTGRA END")
         CALL ogwrit(2,"")
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
! ======================================================================
END SUBROUTINE ogexec
