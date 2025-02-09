SUBROUTINE ogopti(Varacc,Numequ,Finish,Desnor,Calval)
! ======================================================================
! OPTIMISATION PART
! ======================================================================
! I/O | VARACC           | R*8 | ITERATION SCALED DISTANCE ACCUMULATED
! ----------------------------------------------------------------------
! INP | NUMEQU           | I*4 | NUMBER OF EQUALITY CONSTRAINTS
! ----------------------------------------------------------------------
! OUT | FINISH           | I*4 | 0=LIMIT 1=OPTIM
! ----------------------------------------------------------------------
! INP | CALVAL           | EXT | FUNCTION FOR VALUES
!     |                  |     | CALDER (VARVAL, CONVAL)
! ======================================================================
! SUBROUTINES CALLED: INVRGT, INVLFT, ACTEXC, ACTINC, OGEVAL
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   REAL(8) Varacc
! ----------------------------------------------------------------------
   INTEGER(4) Numequ
   INTEGER(4) Finish
   EXTERNAL Calval
! ======================================================================
   INTEGER(4) staflg , faccnt , numcor
   REAL(8) Desnor , foldis , cosimp , cornor , quacor , refdis
   REAL(8) co0 , co1 , co2 , nor
   REAL(8) cosco2 , cosco1
   REAL(8) maxdis , norprv
! ----------------------------------------------------------------------
   INTEGER(4) con , var , cos , act , ind , len , inc
   INTEGER(4) nnn , typ , des , prv , met
   REAL(8) val , max , det , ccc , dis
   REAL(8) fac , del , exc , eps , imp
   REAL(8) bet , tht
   CHARACTER str*256 , nam*256
! ======================================================================
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: cosact
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varvec
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varwrk
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: corvec
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: desder
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: desprv
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: varprv
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: convec
   REAL(8) , DIMENSION(:) , ALLOCATABLE :: conqua
   INTEGER(4) , DIMENSION(:) , ALLOCATABLE :: concor
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
! ----------------------------------------------------------------------
         ALLOCATE (cosact(Numcon))
         ALLOCATE (varvec(Numvar))
         ALLOCATE (varwrk(Numvar))
         ALLOCATE (corvec(Numcon))
         ALLOCATE (desder(Numvar))
         ALLOCATE (desprv(Numvar))
         ALLOCATE (varprv(Numvar))
         ALLOCATE (convec(Numcon+1))
         ALLOCATE (conqua(Numcon+1))
         ALLOCATE (concor(Numcon+1))
! ======================================================================
! OPTIMISATION PART
! ----------------------------------------------------------------------
         cos = Numcon + 1
         des = Numcon + 2
         prv = Numcon + 3
         faccnt = 0
         Conred(1:cos,:) = Conder(1:cos,:)
         Conred(des,:) = Vardes
         desder = Vardes
         Conred(prv,:) = Funvar
         Conder(prv,:) = Funvar
         CALL ogwrit(3,"")
         CALL ogwrit(3,"OPTIMISATION PART")
! ----------------------------------------------------------------------
!      WRITE (STR,'("NUMACT = ",I4)') NUMACT
!      CALL OGWRIT (3,STR)
         numcor = Numact
         concor = Actcon
         Numact = 0
         Conact(1:des) = 0
         spag_nextblock_1 = 2
      CASE (2)
!      DO COR = 1, NUMCOR
!          CON = CONCOR(COR)
!          NAM = CONSTR(CON)
!          LEN = CONLEN(CON)
!          WRITE (STR,'("ACT = ",I4,5X,1X,A)') CON, NAM(1:LEN)
!          CALL OGWRIT (3,STR)
!          CALL OGINCL (CON)
!      ENDDO
! ======================================================================
! ======================================================================
! VECTOR OF STEEPEST ASCENT
! ----------------------------------------------------------------------
         CALL ogwrit(3,"")
         CALL ogwrit(3,"VECTOR OF STEEPEST ASCENT")
         CALL ogwrit(3,"")
         CALL ogwrit(3,"REMOVE AND INCLUDE CONSTRAINTS:")
         CALL ogwrit(3,"")
         CALL ogwrit(3," REM  INC CONSTRAINT")
! ----------------------------------------------------------------------
! REMOVE PASSIVE INEQUALITY CONSTRAINTS
! ----------------------------------------------------------------------
         DO act = Numact , 1 , -1
            con = Actcon(act)
            IF ( Conval(con)<=1D0 ) CYCLE
            nam = Constr(con)
            len = Conlen(con)
            WRITE (str,'(I4,5X,1X,A)') con , nam(1:len)
            CALL ogwrit(2,str)
            CALL ogexcl(act)
            Conact(con) = -1
         ENDDO
! ----------------------------------------------------------------------
! INCLUDE VIOLATED INEQUALITY CONSTRAINTS AND SELECT PASSIVE ONES
! SELECT PASSIVE INEQUALITY CONSTRAINTS
! ----------------------------------------------------------------------
         DO con = 1 , Numcon
            IF ( Contyp(con)==-2 ) CYCLE
            IF ( Conact(con)>0 ) THEN
            ELSEIF ( Contyp(con)==0 ) THEN
               Conact(con) = 0
            ELSEIF ( Conval(con)<-1D0 ) THEN
               Conact(con) = 0
            ELSEIF ( Conval(con)<=+1D0 ) THEN
               Conact(con) = 0
            ELSE
               Conact(con) = -1
            ENDIF
         ENDDO
! ======================================================================
         nnn = 1
         SPAG_Loop_1_1: DO
            nnn = nnn + 1
            IF ( nnn>999 ) THEN
               Finish = 0
               WRITE (str,*) "NNN=" , nnn
               CALL ogwrit(2,str)
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
! ======================================================================
! DERIVATIVES OF MERIT W.R.T. ACTIVE CONSTRAINTS
! ----------------------------------------------------------------------
            CALL ogrigt(-Conred(cos,1:Numact),cosact)
            Desnor = dsqrt(sum(Conred(cos,Numact+1:Numvar)**2))
! ----------------------------------------------------------------------
! CONSTRAINT REMOVAL
! ----------------------------------------------------------------------
            ind = 0
            exc = -1D-12
            max = exc
            DO act = 1 , Numact
               con = Actcon(act)
               IF ( Contyp(con)==0 ) CYCLE
               val = cosact(act)
               fac = dot_product(Conred(con,1:Numvar),Conred(con,1:Numvar))
               fac = dsqrt(fac)
               val = val*fac
               IF ( val>=exc ) CYCLE
               IF ( val>max ) CYCLE
               max = val
               ind = act
            ENDDO
! ----------------------------------------------------------------------
            IF ( ind/=0 ) THEN
               con = Actcon(ind)
               nam = Constr(con)
               len = Conlen(con)
               WRITE (str,'(I4,5X,3(1X,D10.3),1X,A)') con , Desnor , max , Varmax , nam(1:len)
               CALL ogwrit(3,str)
               CALL ogexcl(ind)
               CYCLE
            ENDIF
! ----------------------------------------------------------------------
! CONSTRAINT INCLUSION
! ----------------------------------------------------------------------
            IF ( Desnor/=0D0 ) THEN
! ----------------------------------------------------------------------
               inc = 0
               eps = 1D-03
               max = -1D10
               max = 0D0
! ----------------------------------------------------------------------
               DO con = 1 , Numcon
                  IF ( Contyp(con)==-2 ) CYCLE
                  IF ( Conact(con)/=0 ) CYCLE
                  del = dot_product(Conred(con,Numact+1:Numvar),Conred(cos,Numact+1:Numvar))/Desnor
                  val = dabs(del)*Varmax
                  IF ( val<eps ) CYCLE
                  fac = dot_product(Conred(con,1:Numvar),Conred(con,1:Numvar))
                  fac = dsqrt(fac)
                  del = del/fac
                  IF ( del<0D0 .AND. del<max ) THEN
                     max = del
                     inc = con
                  ENDIF
                  IF ( Contyp(con)/=0 ) CYCLE
                  del = -del
                  IF ( del<0D0 .AND. del<max ) THEN
                     max = del
                     inc = con
                  ENDIF
               ENDDO
! ----------------------------------------------------------------------
               IF ( inc/=0 ) THEN
                  con = inc
                  nam = Constr(con)
                  len = Conlen(con)
                  WRITE (str,'(5X,I4,3(1X,D10.3),1X,A)') con , Desnor , max , Varmax , nam(1:len)
                  CALL ogwrit(3,str)
                  CALL ogincl(inc)
                  CYCLE
               ENDIF
            ENDIF
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! MATCHED INEQUALITY CONSTRAINTS + STEEPEST ASCENT VECTOR NORM
! ----------------------------------------------------------------------
            CALL ogwrit(3,"")
            CALL ogwrit(3,"STATUS OF MATCHED INEQUALITY CONSTRAINTS:")
            CALL ogwrit(3,"")
            CALL ogwrit(3," ACT  PAS CONSTRAINT")
            DO con = 1 , Numcon
               IF ( Contyp(con)==-2 ) CYCLE
               nam = Constr(con)
               len = Conlen(con)
               IF ( Contyp(con)==0 ) THEN
               ELSEIF ( Conact(con)>0 ) THEN
                  WRITE (str,'(I4,5X,1X,A)') con , nam(1:len)
                  CALL ogwrit(3,str)
               ELSEIF ( Conact(con)==0 ) THEN
                  WRITE (str,'(5X,I4,1X,A)') con , nam(1:len)
                  CALL ogwrit(3,str)
               ENDIF
            ENDDO
            CALL ogwrit(3,"")
            WRITE (str,'("STEEPEST ASCENT NORM: ",D13.6)') Desnor
            CALL ogwrit(3,str)
            WRITE (str,'("MAXIMUM DISTANCE....: ",D13.6)') Varmax
            CALL ogwrit(3,str)
! ======================================================================
            Finish = 0
! ======================================================================
! IF CONVERGENCE
! ----------------------------------------------------------------------
            cosimp = Desnor*Varmax
            IF ( dabs(cosimp)<=1D0 ) THEN
               foldis = 0D0
               Finish = 1
               WRITE (str,'("FINAL...............:",1X,D13.6,'//'11X,1(1X,D10.3),1X,D16.9)') foldis , cosimp , Conval(cos) + cosimp
               CALL ogwrit(2,str)
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
! ======================================================================
! IF CONSTRAINT IS HIT
! ----------------------------------------------------------------------
            DO var = 1 , Numvar
               val = Conder(cos,var)
               DO act = 1 , Numact
                  con = Actcon(act)
                  val = val + Conder(con,var)*cosact(act)
               ENDDO
               Vardes(var) = val
            ENDDO
! ----------------------------------------------------------------------
            ind = 0
            dis = 1D10
            DO con = 1 , Numcon
               IF ( Contyp(con)==-2 ) CYCLE
               IF ( Conact(con)/=-1 ) CYCLE
               val = dot_product(Conred(con,Numact+1:Numvar),Conred(cos,Numact+1:Numvar))
               IF ( val==0D0 ) CYCLE
               val = -Conval(con)/val*Desnor
               IF ( val<=0D0 ) CYCLE
               IF ( val>=dis ) CYCLE
               dis = val
               ind = con
            ENDDO
! ----------------------------------------------------------------------
            IF ( ind/=0 ) THEN
               val = dsqrt(sum((Varval-Varref+Vardes*dis/Desnor)**2))
               IF ( val>Varmax ) ind = 0
            ENDIF
! ----------------------------------------------------------------------
            IF ( ind/=0 ) THEN
               IF ( Confix(ind)<=0 ) THEN
                  IF ( val>Varmax*1D-1 ) ind = 0
               ENDIF
            ENDIF
! ----------------------------------------------------------------------
            IF ( ind/=0 ) THEN
               imp = dis*Desnor
               con = ind
               tht = 1D0
               bet = 0D0
               nam = Constr(con)
               len = Conlen(con)
               WRITE (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,11X,1(1X,D10.3),1X,D16.9,22X,1X,I4,1X,A)') dis , imp ,           &
                    & Conval(cos) + cosimp , con , nam(1:len)
               CALL ogwrit(2,str)
               Varacc = Varacc + dis
               Varval = Varval + dis*Vardes/Desnor
               DO con = 1 , Numcon + 1
                  val = dot_product(Conred(con,Numact+1:Numvar),Conred(cos,Numact+1:Numvar))
                  Conval(con) = Conval(con) + val*dis/Desnor
               ENDDO
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            EXIT SPAG_Loop_1_1
         ENDDO SPAG_Loop_1_1
         SPAG_Loop_1_2: DO
! ======================================================================
! ----------------------------------------------------------------------
            CALL ogrigt(-Conred(cos,1:Numact),cosact)
            DO var = 1 , Numvar
               val = Conder(cos,var)
               DO act = 1 , Numact
                  con = Actcon(act)
                  val = val + Conder(con,var)*cosact(act)
               ENDDO
               desprv(var) = val
            ENDDO
            Desnor = dsqrt(sum(desprv**2))
            WRITE (str,'("DESNOR=",D13.6)') Desnor
!      CALL OGWRIT (2,STR)
! ----------------------------------------------------------------------
            CALL ogrigt(-Conred(prv,1:Numact),cosact)
            DO var = 1 , Numvar
               val = Conder(prv,var)
               DO act = 1 , Numact
                  con = Actcon(act)
                  val = val + Conder(con,var)*cosact(act)
               ENDDO
               varprv(var) = val
            ENDDO
            norprv = dsqrt(sum(varprv**2))
            WRITE (str,'("NORPRV=",D13.6)') norprv
!      CALL OGWRIT (2,STR)
! ----------------------------------------------------------------------
            CALL ogrigt(-Conred(des,1:Numact),cosact)
            DO var = 1 , Numvar
               val = desder(var)
               DO act = 1 , Numact
                  con = Actcon(act)
                  val = val + Conder(con,var)*cosact(act)
               ENDDO
               Vardir(var) = val
            ENDDO
            nor = dsqrt(sum(Vardir**2))
            WRITE (str,'("NOR=",D13.6)') nor
!      CALL OGWRIT (2,STR)
! ----------------------------------------------------------------------
! MET = 3: CONJUGATE GRADIENT METHOD
! MET = 2: SPETRAL CONJUGATE GRADIENT METHOD
! MET = 1: MODIFIED SPETRAL CONJUGATE GRADIENT METHOD
! MET = 0: STEEPEST DESCENT METHOD
! ----------------------------------------------------------------------
            met = Optmet
            tht = 1D0
            bet = 0D0
            IF ( met==0 ) THEN
            ELSEIF ( met==1 ) THEN
               varvec = desprv - varprv
               IF ( norprv**2>1D-12 ) THEN
                  tht = -dot_product(Vardir,varvec)/norprv**2
                  bet = Desnor**2/norprv**2
               ENDIF
            ELSEIF ( met==2 ) THEN
               varvec = desprv - varprv
               varwrk = Varref - Vargrd
               val = dot_product(varwrk,varvec)
               fac = dot_product(Vardir,varvec)
               IF ( dabs(val)>1D-12 .AND. dabs(fac)>1D-12 ) THEN
                  tht = -dot_product(varwrk,varwrk)/val
                  varwrk = -varvec*tht - varwrk
                  bet = dot_product(varwrk,desprv)/fac
               ENDIF
            ELSEIF ( met==3 ) THEN
               IF ( norprv/=0D0 ) THEN
                  tht = 1D0
                  bet = Desnor**2/norprv**2
               ENDIF
            ENDIF
!      WRITE (STR,'("THT=",D13.6)') THT
!      CALL OGWRIT (3,STR)
!      WRITE (STR,'("BET=",D13.6)') BET
!      CALL OGWRIT (3,STR)
! ----------------------------------------------------------------------
            eps = 1D-03
!      WRITE (STR,*) "THT/BET=",THT,BET
!      CALL OGWRIT (2,STR)
            Vardes = tht*desprv + bet*Vardir
            Desnor = dsqrt(sum(Vardes**2))
            nor = Desnor
            DO con = 1 , Numcon
               IF ( Contyp(con)==-2 ) CYCLE
               IF ( Conact(con)/=0 ) CYCLE
               del = dot_product(Conder(con,1:Numvar),Vardes(1:Numvar))/nor
               val = dabs(del)*Varmax
               IF ( val<eps ) CYCLE
               fac = dot_product(Conder(con,1:Numvar),Conder(con,1:Numvar))
               del = del/dsqrt(fac)
               nam = Constr(con)
               len = Conlen(con)
               typ = Contyp(con)
               IF ( del<0D0 ) THEN
!          CALL OGINCL (CON)
                  act = Conact(con)
!          WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)')
!     &           CON,ACT,+NOR,VAL,DEL,NAM(1:LEN)
!          CALL OGWRIT (2,STR)
                  IF ( act/=0 ) CYCLE SPAG_Loop_1_2
                  bet = 0D0
                  tht = 1D0
               ENDIF
               IF ( Contyp(con)/=0 ) CYCLE
               del = -del
               IF ( del<0D0 ) THEN
!          CALL OGINCL (CON)
                  act = Conact(con)
!          WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)')
!     &           CON,ACT,-NOR,VAL,DEL,NAM(1:LEN)
!          CALL OGWRIT (2,STR)
                  IF ( act/=0 ) CYCLE SPAG_Loop_1_2
                  bet = 0D0
                  tht = 1D0
                  del = dot_product(Conder(con,1:Numvar),Conder(cos,1:Numvar))
                  val = dabs(del)*Varmax/Desnor
!          WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)')
!     &           CON,TYP,-DESNOR,VAL,DEL,NAM(1:LEN)
!          CALL OGWRIT (2,STR)
               ENDIF
            ENDDO
            cosco1 = dot_product(Vardes,Conder(cos,1:Numvar))/Desnor
            IF ( cosco1<0D0 ) THEN
               WRITE (str,*) "COSCO1=" , cosco1
               CALL ogwrit(2,str)
               bet = 0D0
               tht = 1D0
            ENDIF
            Vardes = tht*desprv + bet*Vardir
            Desnor = dsqrt(sum(Vardes**2))
            EXIT SPAG_Loop_1_2
         ENDDO SPAG_Loop_1_2
         SPAG_Loop_1_3: DO
!      WRITE (STR,*) "THT/BET=",THT,BET
!      CALL OGWRIT (2,STR)
! ======================================================================
! SECOND ORDER DERIVATIVES BY NUMERICAL DIFFERENCING
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
            CALL ogwrit(3,"")
            WRITE (str,'("SECOND ORDER CORRECTION")')
            CALL ogwrit(3,str)
            Vardes = tht*desprv + bet*Vardir
            Desnor = dsqrt(sum(Vardes**2))
! ----------------------------------------------------------------------
! MAXIMUM TRAVEL DISTANCE
! ----------------------------------------------------------------------
            dis = Varmax
            varvec = Varval - Varref
            co0 = dot_product(varvec,varvec) - dis**2
            IF ( co0>=0D0 ) THEN
               dis = 0D0
            ELSE
               co1 = dot_product(Vardes,varvec)
               co2 = Desnor**2
               det = co1**2 - co0*co2
               dis = (dsqrt(det)-co1)/co2
            ENDIF
            dis = dis*Desnor
            maxdis = dis
! ======================================================================
! COMPUTE SECOND ORDER EFFECTS
! ----------------------------------------------------------------------
            nnn = 0
            IF ( Senopt>=+1 ) THEN
               DO con = 1 , Numcon
                  IF ( Contyp(con)==-2 ) CYCLE
                  act = Conact(con)
                  ind = Senact(con)
                  IF ( act==-1 .AND. ind==-1 ) CYCLE
                  IF ( act==0 .AND. ind==0 ) CYCLE
                  IF ( act>0 .AND. ind>0 ) CYCLE
                  nam = Constr(con)
                  len = Conlen(con)
                  WRITE (str,'(I4,1X,I4,1X,I4,1X,A)') con , act , ind , nam(1:len)
                  CALL ogwrit(2,str)
                  nnn = 1
               ENDDO
            ENDIF
            IF ( Senopt<=0 .OR. nnn==1 ) THEN
               fac = Varstp/Desnor
               varvec = Varref + Vardes*Varstp/Desnor
               CALL ogeval(varvec,convec,0,Conder,Calval,Calval)
               conqua = matmul(Conder(1:cos,1:Numvar),Vardes(1:Numvar))
               conqua = 2D0*(convec-Conref-conqua*fac)/fac**2
            ENDIF
            IF ( Senopt==-1 ) THEN
               Senqua = conqua
            ELSEIF ( Senopt>=+1 .AND. nnn==0 ) THEN
               conqua = Senqua
            ENDIF
! ======================================================================
! COMPUTE CORRECTION VECTOR
! ----------------------------------------------------------------------
            DO act = 1 , Numact
               con = Actcon(act)
               corvec(act) = conqua(con)
            ENDDO
            CALL ogleft(corvec,corvec)
! ----------------------------------------------------------------------
            cornor = dsqrt(sum(corvec(1:Numact)**2))*0.5D0/Desnor/Desnor
            CALL ogwrit(3,"")
            WRITE (str,'("STEEPEST ASCENT  NORM: ",D13.6)') Desnor
            CALL ogwrit(3,str)
            WRITE (str,'("ACCUMULATED  DISTANCE: ",D13.6)') Varacc
            CALL ogwrit(3,str)
! ======================================================================
! GET EXTREMUM
! ----------------------------------------------------------------------
            cosco1 = dot_product(Vardes,Conder(cos,1:Numvar))/Desnor
            cosco2 = conqua(cos) - dot_product(Conred(cos,1:Numact),corvec(1:Numact))
            cosco2 = cosco2*0.5D0/Desnor/Desnor
            WRITE (str,*) "COSCO2/COSCO1=" , cosco2 , cosco1
            CALL ogwrit(3,str)
            IF ( cosco1<0D0 ) CALL ogwrit(2,str)
! ----------------------------------------------------------------------
            foldis = 0D0
            quacor = cornor*foldis*foldis
            cosimp = foldis*(cosco1+foldis*cosco2)
            CALL ogwrit(3,"")
            WRITE (str,'(    "STEEPEST ASCENT FOLLOW",'//'  5X,"DISTANCE",'//'  1X,"CORRECTION",'//'  2X,"MERIT_DEL",'//           &
                  &'  6X,"MERIT_VALUE")')
            CALL ogwrit(3,str)
            WRITE (str,'("INITIAL.............:",1X,D13.6,'//'  2(1X,D10.3),1X,D16.9)') foldis , quacor , cosimp , Conval(cos)     &
                 & + cosimp
            CALL ogwrit(3,str)
! ======================================================================
            IF ( cosco2<0D0 ) THEN
               foldis = -0.5D0*cosco1/cosco2
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               WRITE (str,'("MERIT MAXIMUM.......:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') foldis , quacor , cosimp , Conval(cos)    &
                    & + cosimp
               CALL ogwrit(3,str)
               staflg = 1
            ELSEIF ( cosco2>0D0 ) THEN
               foldis = Varmax
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               WRITE (str,'("MERIT CONVEX........:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') foldis , quacor , cosimp , Conval(cos)    &
                    & + cosimp
               CALL ogwrit(3,str)
               staflg = 2
            ELSE
               foldis = Varmax
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               WRITE (str,'("MERIT LINEAR........:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') foldis , quacor , cosimp , Conval(cos)    &
                    & + cosimp
               CALL ogwrit(3,str)
               staflg = 2
            ENDIF
! ======================================================================
! IF MAXIMUM DISTANCE IS HIT
! ----------------------------------------------------------------------
            IF ( foldis>Varmax ) THEN
               foldis = Varmax
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               WRITE (str,'("MAXIMUM DISTANCE....:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') foldis , quacor , cosimp , Conval(cos)    &
                    & + cosimp
               CALL ogwrit(3,str)
               staflg = 2
            ENDIF
! ======================================================================
! IF CONVERGENCE
! ----------------------------------------------------------------------
            IF ( dabs(cosimp)<=1D0 ) THEN
               WRITE (str,'("FINAL...............:",1X,D13.6,'//'  2(1X,D10.3),1X,D16.9,2D11.3)') foldis , quacor , cosimp ,       &
                    & Conval(cos) + cosimp , tht , bet
               CALL ogwrit(2,str)
               IF ( tht/=1D0 .OR. bet/=0D0 ) THEN
                  tht = 1D0
                  bet = 0D0
                  CYCLE
               ENDIF
               foldis = 0D0
               Finish = 1
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
! ======================================================================
! IF REMAINING DISTANCE IS HIT
! ----------------------------------------------------------------------
            IF ( foldis>maxdis ) THEN
               foldis = maxdis
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               WRITE (str,'("REMAINING DISTANCE..:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') foldis , quacor , cosimp , Conval(cos)    &
                    & + cosimp
               CALL ogwrit(3,str)
               staflg = 2
            ENDIF
! ======================================================================
! IF CONSTRAINT IS HIT
! ----------------------------------------------------------------------
            ind = 0
            DO con = 1 , Numcon
               IF ( Contyp(con)==-2 ) CYCLE
               IF ( Conact(con)/=-1 ) CYCLE
               co2 = conqua(con) - dot_product(Conred(con,1:Numact),corvec(1:Numact))
               co1 = dot_product(Conder(con,1:Numvar),Vardes(1:Numvar))
               co0 = Conval(con)*2D0
               IF ( co2/=0D0 ) THEN
                  det = co1**2 - co2*co0
                  IF ( det<0D0 ) CYCLE
                  det = dsqrt(det)
                  val = 1D10
                  fac = (-co1+det)/co2
                  IF ( fac>0D0 .AND. fac<val ) val = fac
                  fac = (-co1-det)/co2
                  IF ( fac>0D0 .AND. fac<val ) val = fac
               ELSEIF ( co1/=0D0 ) THEN
                  val = -co0/co1*0.5D0
               ELSE
                  CYCLE
               ENDIF
               val = val*Desnor
               IF ( val>0D0 .AND. val<foldis ) THEN
                  foldis = val
                  ind = con
               ENDIF
            ENDDO
! ----------------------------------------------------------------------
            IF ( ind/=0 ) THEN
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               con = ind
               nam = Constr(con)
               len = Conlen(con)
               WRITE (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,2(1X,D10.3),1X,D16.9,1X,I4,1X,A)') foldis , quacor , cosimp ,    &
                    & Conval(cos) + cosimp , con , nam(1:len)
               CALL ogwrit(3,str)
               staflg = 3
            ENDIF
! ======================================================================
! UPDATE
! ----------------------------------------------------------------------
            refdis = foldis
            EXIT SPAG_Loop_1_3
         ENDDO SPAG_Loop_1_3
         SPAG_Loop_1_4: DO
            WRITE (str,'("FINAL...............:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') foldis , quacor , cosimp , Conval(cos)       &
                 & + cosimp
            CALL ogwrit(3,str)
! ----------------------------------------------------------------------
            fac = foldis/Desnor
! ----------------------------------------------------------------------
! VARIABLE DELTA
! ----------------------------------------------------------------------
            CALL ogrigt(corvec,cosact)
            dis = 0D0
            DO var = 1 , Numvar
               val = 0D0
               DO act = 1 , Numact
                  ind = Actcon(act)
                  val = val - cosact(act)*Conder(ind,var)
               ENDDO
               varvec(var) = fac*(Vardes(var)+(val*fac*0.5D0))
               dis = dis + varvec(var)*varvec(var)
            ENDDO
            dis = dsqrt(dis)
! ----------------------------------------------------------------------
            WRITE (str,*) "REFDIS=" , refdis
            CALL ogwrit(3,str)
            WRITE (str,*) "FOLDIS=" , foldis
            CALL ogwrit(3,str)
            WRITE (str,*) "DIS=" , dis
            CALL ogwrit(3,str)
            IF ( dis>refdis*1.2D0 .AND. Senopt>0 ) THEN
               faccnt = faccnt + 1
               IF ( faccnt>=10 ) EXIT SPAG_Loop_1_4
               foldis = foldis*0.5D0
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               CYCLE
            ENDIF
! ----------------------------------------------------------------------
! UPDATE VARIABLES
! ----------------------------------------------------------------------
            Varacc = Varacc + foldis
            Varval = Varval + varvec
            ccc = dsqrt(sum((Varval-Varref)**2)) - Varmax**2
            IF ( ccc>=0D0 ) THEN
               WRITE (str,*) "CCC > 0" , ccc
               CALL ogwrit(3,str)
               staflg = 2
            ENDIF
! ======================================================================
! MAXIMUM REACHED: NEXT ITERATION
! ----------------------------------------------------------------------
            IF ( staflg==1 ) THEN
               WRITE (str,'("MERIT MAXIMUM.......:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9,2D11.3)') foldis , quacor , cosimp ,         &
                    & Conval(cos) + cosimp , tht , bet
               CALL ogwrit(2,str)
               IF ( Senopt>0 ) Finish = 1
               EXIT SPAG_Loop_1_4
            ENDIF
! ======================================================================
! MAXIMUM TRAVEL DISTANCE REACHED: NEXT ITERATION
! ----------------------------------------------------------------------
            IF ( staflg==2 ) THEN
               WRITE (str,'("REMAINING DISTANCE..:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9,2D11.3)') foldis , quacor , cosimp ,         &
                    & Conval(cos) + cosimp , tht , bet
               CALL ogwrit(2,str)
               EXIT SPAG_Loop_1_4
            ENDIF
! ======================================================================
! CONSTRAINT HIT: UPDATE CONSTRAINT + CORRECT
! ----------------------------------------------------------------------
            IF ( staflg==3 ) THEN
               nam = Constr(con)
               len = Conlen(con)
               WRITE (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,2(1X,D10.3),1X,D16.9,2D11.3,1X,I4,1X,A)') foldis , quacor ,      &
                    & cosimp , Conval(cos) + cosimp , tht , bet , con , nam(1:len)
               CALL ogwrit(2,str)
               convec = conqua - matmul(Conred(1:cos,1:Numact),corvec(1:Numact))
               convec = convec*fac*0.5D0
               convec = convec + matmul(Conder(1:cos,1:Numvar),Vardes(1:Numvar))
               Conval = Conval + convec*fac
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            EXIT SPAG_Loop_1_4
         ENDDO SPAG_Loop_1_4
         spag_nextblock_1 = 3
      CASE (3)
! ======================================================================
! ----------------------------------------------------------------------
         Funvar = desprv
         Confix = Conact(1:Numcon)
         IF ( Senopt==-1 ) Senact = Conact(1:Numcon)
! ----------------------------------------------------------------------
         DEALLOCATE (cosact)
         DEALLOCATE (varvec)
         DEALLOCATE (varwrk)
         DEALLOCATE (corvec)
         DEALLOCATE (desder)
         DEALLOCATE (desprv)
         DEALLOCATE (varprv)
         DEALLOCATE (convec)
         DEALLOCATE (conqua)
         DEALLOCATE (concor)
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
! ======================================================================
END SUBROUTINE ogopti
