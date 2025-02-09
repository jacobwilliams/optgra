SUBROUTINE ogcorr(Varacc,Finish,Toterr,Norerr,Calval,Calder)
! ======================================================================
! CORRECTION PART
! ======================================================================
! SUBROUTINES CALLED: OGRIGT, OGLEFT, OGEXCL, OGINCL, OGEVAL
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   REAL(8) Varacc
   INTEGER(4) Finish
   EXTERNAL Calval , Calder
! ======================================================================
   INTEGER(4) coritr , numfff , minpri , maxpri , curpri
   REAL(8) cornor , foldis , cstval
   REAL(8) conerr , Toterr , Norerr
   REAL(8) corinv , varvio , conmax , normax
   INTEGER(4) conind , norind , inelop , maxitr
! ----------------------------------------------------------------------
   INTEGER(4) con , var , act , ind , len , cos , stp
   INTEGER(4) typ , cor , pri , vio , fff
   REAL(8) val , fac , upr , del , co2 , co1 , co0 , de2 , dis
   REAL(8) eps , err , dlt , sca , dif
   REAL(8) exc
   CHARACTER str*256 , nam*256
! ======================================================================
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: cosact
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varvec
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varsav
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varcor
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: corvec
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: consav
   INTEGER(4) , DIMENSION(:) , ALLOCATABLE :: conttt
   INTEGER(4) , DIMENSION(:) , ALLOCATABLE :: concor
   INTEGER(4) , DIMENSION(:) , ALLOCATABLE :: coninc
   INTEGER(4) , DIMENSION(:) , ALLOCATABLE :: conhit
   INTEGER(4) , DIMENSION(:) , ALLOCATABLE :: fffcon
   INTEGER(4) , DIMENSION(:) , ALLOCATABLE :: prisav
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
! ----------------------------------------------------------------------
         ALLOCATE (cosact(Numvar))
         ALLOCATE (varvec(Numvar))
         ALLOCATE (varsav(Numvar))
         ALLOCATE (varcor(Numvar))
         ALLOCATE (corvec(Numvar))
         ALLOCATE (consav(Numcon+1))
         ALLOCATE (conttt(Numcon+1))
         ALLOCATE (concor(Numcon))
         ALLOCATE (coninc(Numcon))
         ALLOCATE (conhit(Numcon))
         ALLOCATE (fffcon(Numcon))
         ALLOCATE (prisav(Numcon))
! ======================================================================
! CORRECTION PART
! ----------------------------------------------------------------------
         coritr = 0
         maxitr = 10
         IF ( Senopt/=0 ) maxitr = 1
! ======================================================================
         cos = Numcon + 1
         vio = Numcon + 2
         stp = 0
         Finish = 0
         eps = 1D-03
         dlt = 1D-06
         varvio = Varmax*1D+03
         numfff = Numact
         fffcon = Actcon
         conttt = Contyp
         prisav = Conpri
         concor = 0
! ----------------------------------------------------------------------
         minpri = 1000
         maxpri = -1000
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            minpri = min(minpri,Conpri(con))
            maxpri = max(maxpri,Conpri(con))
         ENDDO
         ind = 0
         IF ( numfff>0 ) ind = 1
         IF ( Senopt>0 ) ind = 1
         minpri = minpri - ind
         CALL ogwrit(3,"")
         CALL ogwrit(3,"PRIORITISE CONSTRAINTS")
         CALL ogwrit(3,"")
         IF ( Senopt<=0 ) THEN
            DO fff = 1 , numfff
               con = fffcon(fff)
               nam = Constr(con)
               len = Conlen(con)
               typ = Contyp(con)
               WRITE (str,'("PRI",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
               CALL ogwrit(3,str)
               Conpri(con) = minpri
            ENDDO
         ENDIF
         IF ( Senopt>0 ) THEN
            DO con = 1 , Numcon
               IF ( Contyp(con)==-2 ) CYCLE
               IF ( Senact(con)<=0 ) CYCLE
               nam = Constr(con)
               len = Conlen(con)
               typ = Contyp(con)
               WRITE (str,'("PRI",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
               CALL ogwrit(3,str)
               Conpri(con) = minpri
            ENDDO
         ENDIF
         CALL ogwrit(3,"")
         spag_nextblock_1 = 2
      CASE (2)
! ======================================================================
! Evaluation loop
! ----------------------------------------------------------------------
         WRITE (str,'("CORITR=",I2,1X,I2)') coritr , maxitr
         CALL ogwrit(3,str)
         CALL ogwrit(3,"")
! ======================================================================
! Inequality loop
! ----------------------------------------------------------------------
         inelop = 2
         IF ( numfff>0 ) inelop = 1
         IF ( Senopt>0 ) inelop = 1
         varsav = Varval
         consav = Conval
         spag_nextblock_1 = 3
      CASE (3)
! ----------------------------------------------------------------------
         WRITE (str,'("INELOP=",I2)') inelop
         CALL ogwrit(3,str)
         CALL ogwrit(3,"")
! ----------------------------------------------------------------------
         Varval = varsav
         Conval = consav
         Contyp = conttt
         IF ( inelop==1 ) THEN
            DO fff = 1 , numfff
               con = fffcon(fff)
               nam = Constr(con)
               len = Conlen(con)
               typ = Contyp(con)
               WRITE (str,'("TYP",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
               CALL ogwrit(3,str)
               Contyp(con) = 0
            ENDDO
         ENDIF
         IF ( inelop==1 .AND. Senopt>0 ) THEN
            DO con = 1 , Numcon
               IF ( Contyp(con)==-2 ) CYCLE
               IF ( Senact(con)<=0 ) CYCLE
               nam = Constr(con)
               len = Conlen(con)
               typ = Contyp(con)
               WRITE (str,'("FIX",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
               CALL ogwrit(2,str)
               Contyp(con) = 0
            ENDDO
         ENDIF
         Numact = 0
         Conact = 0
         Conred(1:Numcon+1,:) = Conder(1:Numcon+1,:)
         Conred(Numcon+2,:) = Vardir
! ----------------------------------------------------------------------
! CHECK CONSTRAINTS
! ----------------------------------------------------------------------
         conerr = 0D0
         Norerr = 0D0
         conind = 0
         norind = 0
         conmax = 0D0
         normax = 0D0
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            typ = Contyp(con)
            val = Conval(con)
            IF ( val<-1D0 ) THEN
               err = dabs(val)
            ELSEIF ( typ/=0 ) THEN
               err = 0D0
            ELSEIF ( val>1D0 ) THEN
               err = dabs(val)
            ELSE
               err = 0D0
            ENDIF
            conerr = conerr + err
            IF ( err>conmax ) THEN
               conind = con
               conmax = err
            ENDIF
            fac = 0D0
            DO var = 1 , Numvar
               IF ( Vartyp(var)==1 ) CYCLE
               fac = fac + Conder(con,var)**2
            ENDDO
            fac = dsqrt(fac)
            IF ( err==0D0 ) THEN
            ELSEIF ( fac/=0D0 ) THEN
               err = err/fac
            ELSE
               CALL ogwrit(0,"")
               CALL ogwrit(0,"ERROR: CONSTRAINT CAN NOT BE SATISFIED")
               WRITE (str,'("CON/VAL= ",I5,1X,D13.6)') con , val
               CALL ogwrit(0,str)
               Finish = 0
               spag_nextblock_1 = 7
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            Norerr = Norerr + err
            IF ( err>normax ) THEN
               norind = con
               normax = err
            ENDIF
         ENDDO
! ----------------------------------------------------------------------
         Toterr = conerr
         CALL ogwrit(3,"")
         WRITE (str,'("NUMFFF/TOTERR/NORERR/COSVAL=",I4,3(1X,D13.6))') numfff , Toterr , Norerr , Conval(Numcon+1)
         CALL ogwrit(2,str)
         CALL ogwrit(3,"")
         WRITE (str,'("MAXIM TOTAL ERROR.: ",D13.6,I6)') conmax , conind
         CALL ogwrit(3,str)
         WRITE (str,'("MAXIM NORM  ERROR.: ",D13.6,I6)') normax , norind
         CALL ogwrit(3,str)
! ----------------------------------------------------------------------
         IF ( Senopt>0 .AND. coritr==maxitr ) THEN
            IF ( cstval==0D0 ) THEN
               Finish = 1
            ELSE
               Finish = 0
            ENDIF
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( coritr==0 .AND. conerr==0D0 ) THEN
            Numact = numfff
            Actcon = fffcon
            Finish = 1
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( coritr/=0 .AND. conerr==0D0 ) THEN
            Finish = 1
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ELSEIF ( coritr==maxitr ) THEN
            Finish = 0
            CALL ogwrit(3,"")
            WRITE (str,'("CORITR=",I2)') coritr
            CALL ogwrit(3,str)
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ELSE
            Finish = 0
            coritr = coritr + 1
         ENDIF

! ======================================================================
! Priority loop
! ----------------------------------------------------------------------
         curpri = minpri
         spag_nextblock_1 = 4
      CASE (4)
! ----------------------------------------------------------------------
         WRITE (str,'("CURPRI=",I2,1X,I2)') curpri , maxpri
         CALL ogwrit(3,str)
         CALL ogwrit(3,"")
! ======================================================================
! MINIMUM NORM CORRECTION
! ----------------------------------------------------------------------
         CALL ogwrit(3,"")
         WRITE (str,'("CORRECTION OF CONSTRAINTS")')
         CALL ogwrit(3,str)
         CALL ogwrit(3,"")
         WRITE (str,*) "INELOP/CURPRI=" , inelop , curpri
         CALL ogwrit(3,str)
! ----------------------------------------------------------------------
         conerr = 0D0
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            conhit(con) = 0
            IF ( Conval(con)<-eps ) THEN
               conerr = conerr + dabs(Conval(con))
            ELSEIF ( Contyp(con)/=0 ) THEN
            ELSEIF ( Conval(con)>+eps ) THEN
               conerr = conerr + dabs(Conval(con))
            ENDIF
         ENDDO
! ----------------------------------------------------------------------
         CALL ogwrit(3,"")
         WRITE (str,'("LINEAR ERROR.: ",D13.6)') conerr
         CALL ogwrit(3,str)
! ----------------------------------------------------------------------
         CALL ogwrit(3,"")
         WRITE (str,'(" ACT  PAS  MOV",'//' " COST___VAL COST___GRD",'//' " DIST___DEL CONSTRAINT")')
         CALL ogwrit(3,str)
         spag_nextblock_1 = 5
      CASE (5)
! ======================================================================
! Move loop
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            coninc(con) = 0
            pri = Conpri(con)
            val = Conval(con)
            typ = Contyp(con)
            act = Conact(con)
            cor = concor(con)
            IF ( act>0 ) CYCLE
            IF ( pri>curpri ) THEN
               Conact(con) = -1
            ELSEIF ( val<-dlt ) THEN
               Conact(con) = -1
            ELSEIF ( val>+dlt ) THEN
               Conact(con) = -1
            ELSE
               Conact(con) = 0
            ENDIF
            IF ( pri>curpri ) THEN
               concor(con) = 0
            ELSEIF ( val<-dlt ) THEN
               concor(con) = -1
            ELSEIF ( typ/=0 ) THEN
               concor(con) = 0
            ELSEIF ( val>+dlt ) THEN
               concor(con) = +1
            ELSE
               concor(con) = 0
            ENDIF
!          IF (ACT .NE. CONACT(CON) .OR. COR .NE. CONCOR(CON)) THEN
!              NAM = CONSTR(CON)
!              LEN = CONLEN(CON)
!              WRITE (STR,'(5X,5X,I4,23X,D10.3,1X,A,4I4)')
!     &        CON, CONVAL(CON),NAM(1:LEN),
!     &        CONACT(CON),CONCOR(CON), ACT, COR
!              CALL OGWRIT (3,STR)
!          ENDIF
         ENDDO
! ======================================================================
! STEEPEST ASCENT VECTOR
! ======================================================================
! MERIT VALUE AND DERIVATIVES
! ----------------------------------------------------------------------
         cstval = 0D0
         varvec = 0D0
         Conred(vio,:) = 0D0
! ----------------------------------------------------------------------
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            IF ( Conpri(con)>curpri ) CYCLE
            IF ( concor(con)==0 ) CYCLE
            fac = concor(con)
            cstval = cstval - Conval(con)*fac
            Conred(vio,:) = Conred(vio,:) - Conred(con,:)*fac
            varvec = varvec - Conder(con,:)*fac
         ENDDO
         SPAG_Loop_1_1: DO
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEEPEST ASCENT VECTOR
! ----------------------------------------------------------------------
            corvec = Conred(vio,:)
! ----------------------------------------------------------------------
            cornor = dsqrt(sum(corvec(Numact+1:Numvar)**2))
! ----------------------------------------------------------------------
! MERIT PARTIAL W.R.T. CONSTRAINTS
! ----------------------------------------------------------------------
            CALL ogrigt(corvec,cosact)
! ----------------------------------------------------------------------
! CONSTRAINT REMOVAL
! ----------------------------------------------------------------------
            ind = 0
            exc = 1D-12
            upr = exc
            DO act = 1 , Numact
               con = Actcon(act)
               IF ( Contyp(con)==0 ) CYCLE
               val = cosact(act)
               IF ( val<=exc ) CYCLE
               IF ( val<upr ) CYCLE
!         IF (VAL .GE. UPR .AND. UPR .GT. 0D0) CYCLE
               upr = val
               ind = act
            ENDDO
! ----------------------------------------------------------------------
            IF ( ind/=0 ) THEN
               con = Actcon(ind)
               nam = Constr(con)
               len = Conlen(con)
               WRITE (str,'(5X,I4,5X,3(1X,D10.3),1X,A)') con , cstval , cornor , upr , nam(1:len)
               CALL ogwrit(3,str)
               CALL ogexcl(ind)
               IF ( coninc(con)>=5 ) THEN
                  WRITE (str,'("OGCORR-WARNING: CONSTRAINT INCLUDED")')
                  CALL ogwrit(1,str)
                  WRITE (str,'("CON/INC/UPR=",2I4,1X,D10.3)') con , coninc(con) , upr
                  CALL ogwrit(1,str)
               ENDIF
               IF ( coninc(con)>=20 ) THEN
                  Finish = 0
                  spag_nextblock_1 = 7
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               CYCLE
            ENDIF
! ----------------------------------------------------------------------
! NORMALISE STEEPEST ASCEND VECTOR
! ----------------------------------------------------------------------
            IF ( cstval<-cornor*varvio ) EXIT SPAG_Loop_1_1
! ----------------------------------------------------------------------
            corinv = 1D0/cornor
            corvec(Numact+1:Numvar) = corvec(Numact+1:Numvar)*corinv
! ----------------------------------------------------------------------
! CONSTRAINT INCLUSION
! ----------------------------------------------------------------------
            ind = 0
            upr = 0D0
! ----------------------------------------------------------------------
            DO con = 1 , Numcon
               IF ( Contyp(con)==-2 ) CYCLE
               IF ( Conpri(con)>curpri ) CYCLE
               IF ( Conact(con)/=0 ) CYCLE
               del = dot_product(Conred(con,Numact+1:Numvar),Conred(vio,Numact+1:Numvar))
               val = dabs(del)*Varmax/cornor
               IF ( val<eps ) CYCLE
               fac = dot_product(Conred(con,1:Numvar),Conred(con,1:Numvar))
               del = del/dsqrt(fac)
               IF ( del<upr ) THEN
                  upr = del
                  ind = con
               ENDIF
               IF ( Contyp(con)/=0 ) CYCLE
               IF ( concor(con)/=0 ) CYCLE
               del = -del
               IF ( del<upr ) THEN
                  upr = del
                  ind = con
               ENDIF
            ENDDO
! ----------------------------------------------------------------------
            IF ( ind/=0 ) THEN
               con = ind
               nam = Constr(con)
               len = Conlen(con)
               WRITE (str,'(I4,5X,5X,3(1X,D10.3),1X,A)') con , cstval , cornor , upr , nam(1:len)
               CALL ogwrit(3,str)
               CALL ogincl(ind)
               coninc(con) = coninc(con) + 1
               CYCLE
            ENDIF
            EXIT SPAG_Loop_1_1
         ENDDO SPAG_Loop_1_1
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
         DO var = 1 , Numvar
            val = varvec(var)
            DO act = 1 , Numact
               con = Actcon(act)
               val = val - Conder(con,var)*cosact(act)
            ENDDO
            varvec(var) = val*corinv
         ENDDO
! ----------------------------------------------------------------------
         varcor = Varval - Varref
         co2 = dot_product(varvec,varvec)
         co1 = dot_product(varvec,varcor)*0.5D0
         co0 = dot_product(varcor,varcor) - Varmax**2
         de2 = co1**2 - co2*co0
         IF ( de2>=0D0 .AND. co2/=0D0 ) THEN
            dis = (dsqrt(de2)-co1)/co2
         ELSE
            dis = 0D0
         ENDIF
! ----------------------------------------------------------------------
         DO var = 1 , Numvar
            fac = varvec(var)
            IF ( fac==0D0 ) CYCLE
            dif = Varval(var) - Varref(var)
            sca = Varmax*1D-0
            val = (dif+sca)/fac
            fac = (dif-sca)/fac
            IF ( fac>val ) val = fac
            IF ( val<dis ) dis = val
         ENDDO
         IF ( dis<0D0 ) dis = 0D0
! ----------------------------------------------------------------------
         foldis = dis
! ======================================================================
! OPTIMISE DIRETION OF STEPPEST ASCENT
! ======================================================================
         IF ( cstval==0D0 ) THEN
! ----------------------------------------------------------------------
!          WRITE (STR,'("CNV=",3(1X,D10.3))') CSTVAL/CORNOR/VARVIO
!          CALL OGWRIT (3,STR)
            WRITE (str,'(4X,5X,5X,3(1X,D10.3))') cstval , cornor , foldis
            CALL ogwrit(3,str)
            CALL ogwrit(3,"")
            IF ( curpri>=maxpri ) THEN
               WRITE (str,'("MAXPRI=",I3)') maxpri
               CALL ogwrit(3,str)
! ======================================================================
! ======================================================================
! MATCHED INEQUALITY CONSTRAINTS + MINIMUM CORRECTION NORM
! ----------------------------------------------------------------------
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
! ======================================================================
! ======================================================================
               IF ( Senopt<=0 ) CALL ogeval(Varval,Conval,0,Conder,Calval,Calder)
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ELSE
               WRITE (str,'("CURPRI=",I3)') curpri
               CALL ogwrit(3,str)
               curpri = curpri + 1
               spag_nextblock_1 = 4
               CYCLE SPAG_DispatchLoop_1
            ENDIF
! ----------------------------------------------------------------------
         ELSEIF ( cstval<-cornor*varvio ) THEN
! ----------------------------------------------------------------------
            WRITE (str,'("CNV=",3(1X,D10.3))') cstval/cornor/varvio
            CALL ogwrit(2,str)
            CALL ogwrit(3,"")
            WRITE (str,'("CSTVAL=",3D10.3)') cstval , cornor , varvio
            CALL ogwrit(3,str)
            IF ( inelop==1 ) THEN
               WRITE (str,'("INELOP=",I3)') inelop
               CALL ogwrit(3,str)
               inelop = inelop + 1
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ELSE
               Finish = 0
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ENDIF
! ----------------------------------------------------------------------
         ENDIF
! ======================================================================
! IF CONSTRAINT IS HIT
! ----------------------------------------------------------------------
         ind = 0
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            IF ( Conact(con)/=-1 ) CYCLE
            val = dot_product(Conred(con,Numact+1:Numvar),corvec(Numact+1:Numvar))
            IF ( val==0D0 ) CYCLE
            val = -Conval(con)/val
            IF ( val<=0D0 ) CYCLE
            IF ( val>=foldis ) CYCLE
            foldis = val
            ind = con
         ENDDO
! ======================================================================
! UPDATE VARIABLES, CONSTRAINTS AND COST FUNCTION
! ----------------------------------------------------------------------
         varvec = varvec*foldis
! ----------------------------------------------------------------------
         Varacc = Varacc + foldis
         Varval = Varval + varvec
! ----------------------------------------------------------------------
         DO con = 1 , Numcon + 1
            val = dot_product(corvec(Numact+1:Numvar),Conred(con,Numact+1:Numvar))
            Conval(con) = Conval(con) + val*foldis
         ENDDO
! ----------------------------------------------------------------------
         cstval = cstval + foldis*cornor
! ======================================================================
! MAXIMUM TRAVEL DISTANCE REACHED: NEXT ITERATION
! ----------------------------------------------------------------------
         IF ( ind==0 ) THEN
            WRITE (str,'("CNV=",3(1X,D10.3))') cstval/cornor/varvio
            CALL ogwrit(3,str)
            WRITE (str,'(4X,5X,5X,3(1X,D10.3))') cstval , cornor , foldis
            CALL ogwrit(3,str)
            IF ( inelop==1 ) THEN
               WRITE (str,'("INELOP=",I3)') inelop
               CALL ogwrit(3,str)
               inelop = inelop + 1
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ELSE
               Numact = 0
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
! ======================================================================
! CONSTRAINT HIT: UPDATE CONSTRAINTS + CORRECT
! ----------------------------------------------------------------------
         con = ind
         nam = Constr(con)
         len = Conlen(con)
         val = Conval(con)
         WRITE (str,'(5X,5X,I4,3(1X,D10.3),1X,A)') con , cstval , cornor , foldis , nam(1:len)
         CALL ogwrit(3,str)
         IF ( conhit(con)>=20 ) THEN
            WRITE (str,'("OGCORR-WARNING: CONSTRAINT HIT")')
            CALL ogwrit(1,str)
            WRITE (str,'("CON/HIT=",2I4)') con , conhit(con)
            CALL ogwrit(1,str)
            Finish = 0
            spag_nextblock_1 = 7
            CYCLE SPAG_DispatchLoop_1
         ENDIF
! ----------------------------------------------------------------------
         conhit(con) = conhit(con) + 1
         spag_nextblock_1 = 5
         CYCLE SPAG_DispatchLoop_1
      CASE (6)
! ======================================================================
! ======================================================================
! MATCHED INEQUALITY CONSTRAINTS + MINIMUM CORRECTION NORM
! ----------------------------------------------------------------------
         CALL ogwrit(3,"")
         WRITE (str,'("CSTVAL:",D13.6)') cstval
         CALL ogwrit(3,str)
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
! ----------------------------------------------------------------------
         CALL ogwrit(3,"")
         CALL ogwrit(3,"STATUS OF VIOLATED CONSTRAINTS:")
         CALL ogwrit(3,"")
         CALL ogwrit(3," CON COST___VAL CONSTRAINT")
         conerr = 0D0
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            nam = Constr(con)
            len = Conlen(con)
            val = Conval(con)
            IF ( val<-dlt ) THEN
               conerr = conerr + dabs(val)
               WRITE (str,'( I4,D11.3,1X,A)') con , val , nam(1:len)
               CALL ogwrit(3,str)
            ELSEIF ( Contyp(con)/=0 ) THEN
            ELSEIF ( val>dlt ) THEN
               conerr = conerr + dabs(val)
               WRITE (str,'( I4,D11.3,1X,A)') con , val , nam(1:len)
               CALL ogwrit(3,str)
            ENDIF
         ENDDO
! ----------------------------------------------------------------------
         CALL ogwrit(3,"")
         WRITE (str,'("LINEAR ERROR.: ",D13.6)') conerr
         CALL ogwrit(3,str)
         spag_nextblock_1 = 7
      CASE (7)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
         Contyp = conttt
         Conpri = prisav
! ----------------------------------------------------------------------
         DEALLOCATE (cosact)
         DEALLOCATE (varvec)
         DEALLOCATE (varsav)
         DEALLOCATE (varcor)
         DEALLOCATE (corvec)
         DEALLOCATE (consav)
         DEALLOCATE (conttt)
         DEALLOCATE (concor)
         DEALLOCATE (coninc)
         DEALLOCATE (conhit)
         DEALLOCATE (fffcon)
         DEALLOCATE (prisav)
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
! ======================================================================
END SUBROUTINE ogcorr
