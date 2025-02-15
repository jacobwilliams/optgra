!****************************************************************************************************
!>
!  Near-linear optimisation tool tailored for s/c trajectory design.

!TODO:
! * all these SET functions should be eliminated and all vars set in the init function
! * remove the STOP statement in one subroutine.
! * add an optgra class. move all methods into this class.
! * replace global module variables with class variables.
! * see about removing the SPAG DISPATCH loop thing which I don't like (study original code).
! * make real kind changable via compiler directive (real32, real64, real128)
! * get example working.

module optgra_module

   use iso_fortran_env, only: wp => real64, ip => int32

   implicit none

   abstract interface
      subroutine calval_f(varvec,Valcon,i)
         !! FUNCTION FOR VALUES
         !! INPUT AND OUTPUT NOT SCALED
         import :: wp
         implicit none
         real(wp), dimension(:), intent(in) :: varvec !! size is Numvar
         real(wp), dimension(:), intent(out) :: Valcon !! size is Numcon+1
         integer, intent(in) :: i !! JW: THIS IS NOT DOCUMENTED ?
      end subroutine calval_f
      subroutine calder_f(varvec,convec,Dercon)
         !! FUNCTION FOR VALUES AND DERIVATIVES
         !! INPUT AND OUTPUT NOT SCALED
         import :: wp
         implicit none
         real(wp), dimension(:), intent(in) :: varvec !! size is Numvar
         real(wp), dimension(:), intent(out) :: convec !! size is Numcon+1
         real(wp), dimension(:,:), intent(out) :: Dercon !! size is Numcon+1,Numvar
      end subroutine calder_f

   end interface

   !.... this was the include file (minus the common blocks)
   !     should be moving into a class to make it threadsafe .... TODO  -JW
! ======================================================================
! INCLUDE FILE OGDATA
! ======================================================================
! CONTAINS PARAMETERS AND COMMON DATA OF THE OPTIMISATION
! ======================================================================
      integer(ip),   PARAMETER :: MAXSTR = 80
! ======================================================================
      integer(ip)                            :: NUMVAR = 0
      real(wp),      DIMENSION(:  ), allocatable :: VARVAL
      integer(ip),   DIMENSION(:  ), allocatable :: VARTYP
      real(wp),      DIMENSION(:  ), allocatable :: VARSCA
      CHARACTER(len=maxstr),    DIMENSION(:  ), allocatable :: VARSTR
      integer(ip),   DIMENSION(:  ), allocatable :: VARLEN
      real(wp),      DIMENSION(:  ), allocatable :: VARREF
      real(wp),      DIMENSION(:  ), allocatable :: VARDES
      real(wp),      DIMENSION(:  ), allocatable :: VARGRD
      real(wp),      DIMENSION(:  ), allocatable :: VARDIR
      real(wp),      DIMENSION(:  ), allocatable :: FUNVAR
      real(wp),      DIMENSION(:  ), allocatable :: SENVAR
! ----------------------------------------------------------------------
      integer(ip)                            :: NUMCON = 0
      real(wp),      DIMENSION(:  ), allocatable :: CONVAL
      integer(ip),   DIMENSION(:  ), allocatable :: CONTYP
      integer(ip),   DIMENSION(:  ), allocatable :: CONPRI
      real(wp),      DIMENSION(:  ), allocatable :: CONSCA
      CHARACTER(len=maxstr),    DIMENSION(:  ), allocatable :: CONSTR
      integer(ip),   DIMENSION(:  ), allocatable :: CONLEN
      real(wp),      DIMENSION(:  ), allocatable :: CONREF
      real(wp),      DIMENSION(:  ), allocatable :: SENQUA
      real(wp),      DIMENSION(:  ), allocatable :: SENCON
      real(wp),      DIMENSION(:  ), allocatable :: SENDEL
      integer(ip),   DIMENSION(:  ), allocatable :: SENACT
! ----------------------------------------------------------------------
      integer(ip)                            :: OPTMET = 0
      integer(ip)                            :: MAXITE = 0
      integer(ip)                            :: CORITE = 0
      integer(ip)                            :: OPTITE = 0
      integer(ip)                            :: DIVITE = 0
      integer(ip)                            :: CNVITE = 0
      real(wp)                               :: VARMAX = 0
      real(wp)                               :: VARSND = 0
      real(wp)                               :: VARSTP = 0
! ----------------------------------------------------------------------
      integer(ip)                            :: VARDER = 0
      real(wp),      DIMENSION(:  ), allocatable :: VARPER
! ----------------------------------------------------------------------
      integer(ip)                            :: LOGLUN = 0  ! log file unit
      integer(ip)                            :: LOGLEV = 0  ! log level
! ----------------------------------------------------------------------
      integer(ip)                            :: LOGLUP = 0  ! pygmo log file unit
      integer(ip)                            :: VERBOS = 0  ! pygmo verbosity
      integer(ip)                            :: FEVALS = 0  ! pygmo: number of const fun evals
      integer(ip)                            :: PYGFLA = 0  ! pygmo: flag indicating status of optimisation
      integer(ip)                            :: NUMITE = 0  ! number of iterations
! ----------------------------------------------------------------------
      integer(ip)                            :: MATLEV = 0
! ----------------------------------------------------------------------
      integer(ip)                            :: TABLUN = 0
      integer(ip)                            :: TABLEV = 0
! ----------------------------------------------------------------------
      integer(ip)                            :: SENOPT = 0
! ----------------------------------------------------------------------
      integer(ip)                            :: NUMACT = 0
      integer(ip),   DIMENSION(:  ), allocatable :: ACTCON
      integer(ip),   DIMENSION(:  ), allocatable :: CONFIX
      integer(ip),   DIMENSION(:  ), allocatable :: CONACT
      real(wp),      DIMENSION(:,:), allocatable :: CONDER
      real(wp),      DIMENSION(:,:), allocatable :: CONRED
      real(wp),      DIMENSION(:,:), allocatable :: SENDER
      integer(ip)                            :: CONVER = 0
      integer(ip),   DIMENSION(:  ), allocatable :: CONOPT
! ----------------------------------------------------------------------

contains
!****************************************************************************************************

SUBROUTINE mul2m(A1,M1,K1,L1,N1,A2,M2,K2,L2,N2,A,M,K,L,N)

   !! Matrix multiply.
   !!
   !! `A(K:K+N1,L:L+N) = A1(K1:K1+N1,L1:L1+N2) * A2(K2:K2+N2,L2:L2+N3)`

   IMPLICIT NONE

   integer,intent(in) :: m1, m2, m, k, k1, k2, l, l1 , l2 , n , n1 , n2
   real(wp),intent(out) :: A(M,*)
   real(wp),intent(in) :: A1(M1,*)
   real(wp),intent(in) :: A2(M2,*)

   real(wp) :: f1 , f2
   integer(ip) :: i , i1 , i2 , ic , ir

   DO i1 = K , K + N1 - 1
      DO i = L , L + N - 1
         A(i1,i) = 0D0
      ENDDO
   ENDDO

   DO i1 = 0 , N1 - 1
      DO i2 = 0 , N2 - 1
         IF ( K1>=0 ) THEN
            f1 = A1(i1+K1,i2+L1)
         ELSE
            f1 = A1(i2-K1,i1+L1)
         ENDIF
         IF ( f1/=0D0 ) THEN
            DO i = 0 , N - 1
               IF ( K2>=0 ) THEN
                  f2 = A2(i2+K2,i+L2)
               ELSE
                  f2 = A2(i-K2,i2+L2)
               ENDIF
               IF ( f2/=0D0 ) THEN
                  f2 = f2*f1
                  ic = i1 + K
                  ir = i + L
                  A(ic,ir) = A(ic,ir) + f2
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO

END SUBROUTINE mul2m

SUBROUTINE mulvs(X,A,Z,Kd)
   !! Scalar Vector multiply.
   !!
   !! `Z (1:KD) = X (1:KD) * A`
   IMPLICIT NONE
   real(wp),intent(in) :: A !! SCALAR
   real(wp),intent(in) :: X(*) !! VECTOR
   real(wp),intent(out) :: Z(*) !! VECTOR
   integer(ip),intent(in) :: Kd !! NUMBER OF ELEMENTS TO BE USED
   integer(ip) :: i
   DO i = 1 , Kd
      Z(i) = X(i)*A
   ENDDO
END SUBROUTINE mulvs

SUBROUTINE ogcdel(Delcon)
   !! DEFINE CONSTRAINT + MERIT CONVERGENCE THRESHOLDS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW
   IMPLICIT NONE

   real(wp),intent(in) :: Delcon(Numcon+1) !! CONSTRAINTS DELTAS (1:NUMCON)

   integer(ip) :: con

   DO con = 1 , Numcon
      Sendel(con) = Delcon(con)
   ENDDO
END SUBROUTINE ogcdel

SUBROUTINE ogclos()
   !! DEALLOCATION OF ARRAYS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   ! VARIABLES
   DEALLOCATE (Varval)
   DEALLOCATE (Vartyp)
   DEALLOCATE (Varsca)
   DEALLOCATE (Varstr)
   DEALLOCATE (Varlen)
   DEALLOCATE (Varref)
   DEALLOCATE (Vardes)
   DEALLOCATE (Vargrd)
   DEALLOCATE (Vardir)
   DEALLOCATE (Funvar)
   DEALLOCATE (Senvar)

   ! CONSTRAINTS
   DEALLOCATE (Conval)
   DEALLOCATE (Contyp)
   DEALLOCATE (Conpri)
   DEALLOCATE (Consca)
   DEALLOCATE (Constr)
   DEALLOCATE (Conlen)
   DEALLOCATE (Conref)
   DEALLOCATE (Senqua)
   DEALLOCATE (Sencon)
   DEALLOCATE (Sendel)
   DEALLOCATE (Senact)

   ! DERIVATIVES
   DEALLOCATE (Varper)

   ! WORKING VECTORS
   DEALLOCATE (Actcon)
   DEALLOCATE (Confix)
   DEALLOCATE (Conact)
   DEALLOCATE (Conder)
   DEALLOCATE (Conred)
   DEALLOCATE (Sender)
   DEALLOCATE (Conopt)

END SUBROUTINE ogclos

SUBROUTINE ogcorr(Varacc,Finish,Toterr,Norerr,Calval,Calder)

   !! CORRECTION PART
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp) Varacc
   integer(ip) Finish
   procedure(Calval_f) :: Calval
   procedure(Calder_f) :: Calder

   integer(ip) coritr , numfff , minpri , maxpri , curpri
   real(wp) cornor , foldis , cstval
   real(wp) conerr , Toterr , Norerr
   real(wp) corinv , varvio , conmax , normax
   integer(ip) conind , norind , inelop , maxitr
! ----------------------------------------------------------------------
   integer(ip) con , var , act , ind , len , cos , stp
   integer(ip) typ , cor , pri , vio , fff
   real(wp) val , fac , upr , del , co2 , co1 , co0 , de2 , dis
   real(wp) eps , err , dlt , sca , dif
   real(wp) exc
   CHARACTER(len=256) :: str , nam
! ======================================================================
   real(wp) , DIMENSION(:) , ALLOCATABLE :: cosact
   real(wp) , DIMENSION(:) , ALLOCATABLE :: varvec
   real(wp) , DIMENSION(:) , ALLOCATABLE :: varsav
   real(wp) , DIMENSION(:) , ALLOCATABLE :: varcor
   real(wp) , DIMENSION(:) , ALLOCATABLE :: corvec
   real(wp) , DIMENSION(:) , ALLOCATABLE :: consav
   integer(ip) , DIMENSION(:) , ALLOCATABLE :: conttt
   integer(ip) , DIMENSION(:) , ALLOCATABLE :: concor
   integer(ip) , DIMENSION(:) , ALLOCATABLE :: coninc
   integer(ip) , DIMENSION(:) , ALLOCATABLE :: conhit
   integer(ip) , DIMENSION(:) , ALLOCATABLE :: fffcon
   integer(ip) , DIMENSION(:) , ALLOCATABLE :: prisav
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
               err = abs(val)
            ELSEIF ( typ/=0 ) THEN
               err = 0D0
            ELSEIF ( val>1D0 ) THEN
               err = abs(val)
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
            fac = sqrt(fac)
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
               conerr = conerr + abs(Conval(con))
            ELSEIF ( Contyp(con)/=0 ) THEN
            ELSEIF ( Conval(con)>+eps ) THEN
               conerr = conerr + abs(Conval(con))
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
            cornor = sqrt(sum(corvec(Numact+1:Numvar)**2))
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
               val = abs(del)*Varmax/cornor
               IF ( val<eps ) CYCLE
               fac = dot_product(Conred(con,1:Numvar),Conred(con,1:Numvar))
               del = del/sqrt(fac)
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
            dis = (sqrt(de2)-co1)/co2
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
               conerr = conerr + abs(val)
               WRITE (str,'( I4,D11.3,1X,A)') con , val , nam(1:len)
               CALL ogwrit(3,str)
            ELSEIF ( Contyp(con)/=0 ) THEN
            ELSEIF ( val>dlt ) THEN
               conerr = conerr + abs(val)
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

SUBROUTINE ogcpri(Pricon)

   !! DEFINE CONTRAINT PRIORITY
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Pricon(Numcon+1) !! CONSTRAINTS PRIORITY (1:NUMCON)
                                              !! -> 1-N

   integer(ip) con

   DO con = 1 , Numcon + 1
      Conpri(con) = Pricon(con)
   ENDDO

END SUBROUTINE ogcpri

SUBROUTINE ogcsca(Scacon)

   !! DEFINE CONSTRAINT + MERIT CONVERGENCE THRESHOLDS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp),intent(in) :: Scacon(Numcon+1) !! CONSTRAINTS CONVER THRESHOLD (1:NUMCON)
                                           !! MERIT       CONVER THRESHOLD (1+NUMCON)

   integer(ip) con

   DO con = 1 , Numcon + 1
      Consca(con) = Scacon(con)
   ENDDO

END SUBROUTINE ogcsca

SUBROUTINE ogcstr(Strcon,Lencon)

   !! DEFINE CONSTRAINT + MERIT STRING
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   character(len=maxstr),intent(in) :: Strcon(Numcon+1) !! CONIABLES NAME STRING
   integer(ip),intent(in) :: Lencon(Numcon+1) !! CONIABLES NAME LENGTH

   integer(ip) con , len

   DO con = 1 , Numcon + 1
      len = min(Lencon(con),maxstr)
      Constr(con) = Strcon(con)
      Conlen(con) = len
   ENDDO

END SUBROUTINE ogcstr

SUBROUTINE ogctyp(Typcon)

   !! DEFINE CONTRAINT + MERIT TYPE
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Typcon(Numcon+1) !! CONSTRAINTS TYPE (1:NUMCON)
                                              !! -> 1=GTE -1=LTE 0=EQU -2=DERIVED DATA
                                              !! MERIT       TYPE (1+NUMCON)
                                              !! -> 1=MAX -1=MIN

   integer(ip) con

   DO con = 1 , Numcon + 1
      Contyp(con) = Typcon(con)
   ENDDO

END SUBROUTINE ogctyp

SUBROUTINE ogderi(Dervar,Pervar)

   !! DEFINE COMPUTATION OF DERIVATIVES
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Dervar !! DERIVATIVES COMPUTATION MODE
                                    !!  * 1: USER DEFINED
                                    !!  * 2: NUMERIC WITH DOUBLE DIFFERENCING
                                    !!  * 3: NUMERIC WITH SINGLE DIFFERENCING
   real(wp),intent(in) :: Pervar(Numvar) !! VARIABLES PERTURBATION FOR DERIVATIVES
                                         !! -> NOT SCALED

   integer(ip) var

   Varder = Dervar

   DO var = 1 , Numvar
      Varper(var) = Pervar(var)
   ENDDO

END SUBROUTINE ogderi

SUBROUTINE ogdist(Maxvar,Sndvar)

   !! DEFINE OPTIMISATION CONTROL PARAMETERS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp),intent(in) :: Maxvar !! MAXIMUM DISTANCE PER ITERATION
                                 !!  -> SCALED
   real(wp),intent(in) :: Sndvar !! PERTURBATION FOR 2ND ORDER DERIVATIVES
                                 !!  -> SCALED

   Varmax = Maxvar
   Varsnd = Sndvar

END SUBROUTINE ogdist

SUBROUTINE ogeval(Valvar,Valcon,Dervar,Dercon,calval,calder)

   !! COMPUTES SCALED CONTRAINTS+MERIT AND DERIVATIVES
   !! FROM     SCALED VARIABLES
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp),intent(in) :: Valvar(Numvar)
   real(wp),intent(out) :: Valcon(Numcon+1)
   integer(ip),intent(in) :: Dervar !! DERIVATIVES COMPUTATION MODE
                                    !!  * 0: VALUES ONLY
                                    !!  * 1: USER DEFINED
                                    !!  * 2: NUMERIC WITH DOUBLE DIFFERENCING
                                    !!  * 3: NUMERIC WITH SINGLE DIFFERENCING
   real(wp),intent(out) :: Dercon(Numcon+1,Numvar)
   procedure(Calval_f) :: Calval !! FUNCTION FOR VALUES
                                 !! -> CALVAL (VALVAR, VALCON)
                                 !! -> INPUT AND OUTPUT NOT SCALED
   procedure(Calder_f) :: Calder !! FUNCTION FOR VALUES AND DERIVATIVES
                                 !! -> CALDER (VALVAR, VALCON, DERCON)
                                 !! -> INPUT AND OUTPUT NOT SCALED

   integer(ip) var , con , cod , len , ind , numvio
   real(wp) val , sca , fac , per , sav , der , err , conerr , convio
   character(len=3) :: typ
   character(len=3) :: sta
   character(len=maxstr) :: nam
   character(len=256) :: str

   real(wp) ggg(4,4) , bbb(4) , vvv(4) , objval
! ======================================================================
   real(wp) , DIMENSION(:) , ALLOCATABLE :: varvec
   real(wp) , DIMENSION(:) , ALLOCATABLE :: convec
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
      IF ( cod==0 .AND. con<=Numcon .AND. abs(val)>1D0 ) THEN
         sta = "VIO"
         err = abs(val)
         numvio = numvio + 1
      ENDIF
      IF ( cod/=0 .AND. con<=Numcon .AND. cod/=-2 .AND. -val>1D0 ) THEN
         sta = "VIO"
         err = abs(val)
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
   convio = sqrt(convio)
   ! CALL ogpwri(objval,numvio,convio,Dervar)  ! JW : this routine doesn't have a Dervar dummy arg
   CALL ogpwri(objval,numvio,convio) ! JW : replaced with this
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
         CALL calval(varvec,Dercon(1:,var),0)  ! JW added : here
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
         CALL calval(varvec,Dercon(1:,var),0)  ! JW added : here
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

SUBROUTINE ogexcl(Exc)

   !! REMOVE CONSTRAINT TO ACTIVE SET AND REDUCES DERIVATIVES
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Exc !! CONSTRAINT TO BE REMOVED
                                 !! SEQUENCE NUMBER IN ACTIVE LIST

   real(wp) val , bet , gam
   integer(ip) row , col , act , con
   CHARACTER(len=256) :: str
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
      val = sqrt(val)
      IF ( Conred(con,act)>0D0 ) val = -val
      IF ( abs(val)<1D-15 ) THEN
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

SUBROUTINE ogexec(Valvar,Valcon,Finopt,Finite,Calval,Calder)
!! Main routine.
!!
!! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp),intent(inout) :: Valvar(Numvar) !! VARIABLES VALUE
                                            !! -> NOT SCALED
   real(wp),intent(out) :: Valcon(Numcon+1) !! CONSTRAINTS VALUE (1:NUMCON)
                                            !! MERIT       VALUE (1+NUMCON)
                                            !! -> NOT SCALED
   integer(ip),intent(out) :: Finopt !! TERMINATION STATUS
                                     !!
                                     !!  * 1=    MATCHED &     OPTIMAL
                                     !!  * 2=    MATCHED & NOT OPTIMAL
                                     !!  * 3=NOT MATCHED & NOT OPTIMAL
                                     !!  * 4=NOT FEASIBL & NOT OPTIMAL
   integer(ip),intent(out) :: Finite
   procedure(Calval_f) :: Calval !! FUNCTION FOR VALUES
                                 !! -> CALDER (VALVAR, VALCON)
                                 !! -> INPUT AND OUTPUT NOT SCALED
   procedure(Calder_f) :: Calder !! FUNCTION FOR VALUES AND DERIVATIVES
                                 !! -> CALDER (VALVAR, VALCON, CONDER)
                                 !! -> INPUT AND OUTPUT NOT SCALED

   integer(ip) :: finish , itecor , iteopt
   integer(ip) :: var , con , typ , len , num , numvio
   real(wp) :: val , sca , red , der , fac , old , convio
   CHARACTER(len=256) :: str , nam

   integer(ip) numequ , itediv , itecnv
   real(wp) varacc , cosnew , cosold , varsav , meamer
   real(wp) conerr , desnor , norerr , meaerr

   real(wp) , DIMENSION(:) , ALLOCATABLE :: varsum
   real(wp) , DIMENSION(:) , ALLOCATABLE :: varcor
   real(wp) , DIMENSION(:) , ALLOCATABLE :: concor
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
         IF ( Tablev>=1 ) WRITE (Tablun,'("ITER",1X,"OPT",1X,*(1X,I10))') (var,var=1,Numvar) , (con,con=1,Numcon)
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
               !CALL ogpwri_end(Numite,-Valcon(Numcon+1),numvio,convio,1) ! JW : another argument missmatch
               CALL ogpwri_end(-Valcon(Numcon+1),numvio,convio) ! JW : replaced with this

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
                     IF ( abs(der-red)<1D-2 ) CYCLE
                     IF ( der/=0D0 ) THEN
                        fac = red/der
                     ELSE
                        fac = 0D0
                     ENDIF
                     IF ( abs(fac-1D0)<1D-2 ) CYCLE
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
            IF ( Tablev>=1 ) WRITE (Tablun,'(I4,1X,"COR",1X,*(1X,D10.3))') Numite , (Varval(var),var=1,Numvar) , (Conval(con),con=1,Numcon)
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
            CALL ogopti(varacc,numequ,finish,desnor,Calval,Calder)
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

SUBROUTINE oggsst(Varsen,Quasen,Consen,Actsen,Dersen,Actsav,Consav,Redsav,Dersav,Actnum)

   !! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
   !!
   !! Function to get sensitivity state data, necessary for serialization.
   !! Do not use this directly except in serialization routines
   !!
   !! 2021/07/19 | M. von Looz | NEW

   IMPLICIT NONE

   real(wp),intent(out) :: Varsen(Numvar) !! STORED VARIABLES VALUE
   real(wp),intent(out) :: Quasen(Numcon+1) !! STORED CONSTRAINTS CORRECTION VECTOR
   real(wp),intent(out) :: Consen(Numcon+1) !! STORED CONSTRAINTS VALUE
   integer(ip),intent(out) :: Actsen(Numcon+1) !! STORED CONSTRAINTS ACTIVE
   real(wp),intent(out) :: Dersen(Numcon+1,Numvar) !! STORED DERIVATIVE
   integer(ip),intent(out) :: Actsav(Numcon+1) !! STORED ACTIVE CONSTRAINTS
   integer(ip),intent(out) :: Consav(Numcon+4) !! STORED ACTIVE CONSTRAINTS
   real(wp),intent(out) :: Redsav(Numcon+3,Numvar) !! STORED DERIVATIVE
   real(wp),intent(out) :: Dersav(Numcon+3,Numvar) !! STORED DERIVATIVE
   integer(ip),intent(out) :: Actnum !! NUMBER OF ACTIVE CONSTRAINTS

   integer(ip) var , con

   ! Variable values saved for sensitivity
   Actnum = Numact

   DO var = 1 , Numvar
      Varsen(var) = Senvar(var)
   ENDDO

   DO con = 1 , Numcon + 1
      Quasen(con) = Senqua(con)
      Consen(con) = Sencon(con)
      Actsen(con) = Senact(con)
      DO var = 1 , Numvar
         Dersen(con,var) = Sender(con,var)
      ENDDO
   ENDDO

   ! Temporary status saved of which constraints are active
   DO con = 1 , Numcon + 1
      Actsav(con) = Actcon(con)
   ENDDO

   DO con = 1 , Numcon + 4
      Consav(con) = Conact(con)
   ENDDO

   DO con = 1 , Numcon + 3
      DO var = 1 , Numvar
         Redsav(con,var) = Conred(con,var)
         Dersav(con,var) = Conder(con,var)
      ENDDO
   ENDDO

END SUBROUTINE oggsst

SUBROUTINE ogincl(Inc)

   !! ADDS CONSTRAINT TO ACTIVE SET AND REDUCES DERIVATIVES
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Inc !! CONSTRAINT TO BE INCLUDED

   real(wp) val , fac , gam , sav , max
   integer(ip) row , col , ind , lst
   CHARACTER(len=256) :: str

   ! GENERAL
   Numact = Numact + 1

   ! PERMUTATION TO GET ZERO DERIVATIVES AT END FOR NEW ACTIVE CONSTRAINT
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

   ! PERMUTATION TO GET MAXIMUM PIVOT
   ind = Numact
   max = abs(Conred(Inc,ind))
   DO col = Numact + 1 , lst
      val = abs(Conred(Inc,col))
      IF ( val>max ) THEN
         ind = col
         max = val
      ENDIF
   ENDDO

   IF ( ind/=Numact ) THEN
      DO row = 1 , Numcon + 3
         IF ( Conact(row)<=0 ) THEN
            sav = Conred(row,ind)
            Conred(row,ind) = Conred(row,Numact)
            Conred(row,Numact) = sav
         ENDIF
      ENDDO
   ENDIF

   ! UPDATE LIST OF ACTIVE CONSTRAINTS
   Actcon(Numact) = Inc
   Conact(Inc) = Numact

   ! REDUCE FOR NEW ACTIVE CONSTRAINT
   IF ( abs(Conred(Inc,Numact))<1D-12 ) THEN
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

   val = sqrt(sum(Conred(Inc,Numact:lst)**2))
   IF ( Conred(Inc,Numact)>0D0 ) val = -val

   Conred(Inc,Numact) = Conred(Inc,Numact) - val

   sav = Conred(Inc,Numact)
   fac = 1D0/sav
   Conred(Inc,Numact:lst) = Conred(Inc,Numact:lst)*fac

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

   Conred(Inc,Numact) = val
   Conred(Inc,Numact+1:lst) = 0D0

END SUBROUTINE ogincl

SUBROUTINE oginit(Varnum,Connum)

   !! ALLOCATION OF ARRAYS AND INITIALISATION OF PARAMETERS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Varnum !! NUMBER OF VARIABLES
   integer(ip),intent(in) :: Connum !! NUMBER OF CONSTRAINTS

   integer(ip) var , con

   ! VARIABLES
   Numvar = Varnum

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

   ! CONSTRAINTS
   Numcon = Connum

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

   ! CONTROL
   Optmet = 2
   Maxite = 10
   Corite = 10
   Optite = 10
   Divite = 10
   Cnvite = 10
   Varmax = 10D0
   Varsnd = 1D0

   ! DERIVATIVES
   Varder = 1
   ALLOCATE (Varper(Numvar))
   DO var = 1 , Numvar
      Varper(var) = 1D-03
   ENDDO

   ! LOG FILE
   Loglun = 6
   Loglev = 1

   ! PYGMO LOG FILE
   Loglup = 7
   Loglev = 0

   ! MATLAB CONSOLE
   Matlev = 0

   ! TABLE FILE
   Tablun = 6
   Tablev = 0

   ! LINEAR OPTIMISATION MODE
   Senopt = 0

   ! WORKING VECTORS
   ALLOCATE (Actcon(Numcon+1))
   ALLOCATE (Confix(Numcon))
   ALLOCATE (Conact(Numcon+4))
   ALLOCATE (Conder(Numcon+3,Numvar))
   ALLOCATE (Conred(Numcon+3,Numvar))
   ALLOCATE (Sender(Numcon+3,Numvar))
   ALLOCATE (Conopt(Numcon+1))
   Numact = 0
   Actcon = 0
   Conact = 0
   Confix = 0
   Conder = 0D0
   Conred = 0D0
   Conopt = 0

END SUBROUTINE oginit

SUBROUTINE ogiter(Itemax,Itecor,Iteopt,Itediv,Itecnv)

   !! DEFINE OPTIMISATION CONTROL PARAMETERS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Itemax !! MAXIMUM NUMBER OF ITERATIONS
   integer(ip),intent(in) :: Itecor
   integer(ip),intent(in) :: Iteopt
   integer(ip),intent(in) :: Itediv
   integer(ip),intent(in) :: Itecnv

   Maxite = Itemax
   Corite = Itecor
   Optite = Iteopt
   Divite = Itediv
   Cnvite = Itecnv
   IF ( Corite>Maxite ) Corite = Maxite
   IF ( Optite>Maxite ) Optite = Maxite
   IF ( Divite>Corite ) Divite = Corite
   IF ( Cnvite>Optite ) Cnvite = Optite

END SUBROUTINE ogiter

SUBROUTINE ogleft(Actinp,Actout)

   !! LEFT-MULTIPLIES VECTOR LOWER TRIANGULAR MATRIX OBTAINED BY REDUCTION
   !! AND SUBSEQUENT INVERSION OF DERIVATIVES OF ACTIVE CONSTRAINTS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp),intent(in) :: Actinp(Numcon) !! VECTOR INITAL
   real(wp),intent(out) :: Actout(Numcon) !! VECTOR FINAL (MAY BE SAME AS ACTINP)

   integer(ip) row , col , act
   real(wp) val

   DO act = 1 , Numact
      row = Actcon(act)
      val = Actinp(act)
      DO col = 1 , act - 1
         val = val - Conred(row,col)*Actout(col)
      ENDDO
      Actout(act) = val/Conred(row,act)
   ENDDO

END SUBROUTINE ogleft

SUBROUTINE ogomet(Metopt)

   !! DEFINE OPTIMISATION CONTROL PARAMETERS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Metopt !! OPTIMISATION METHOD:
                                    !!  * 3: CONJUGATE GRADIENT METHOD
                                    !!  * 2: SPECTRAL CONJUGATE GRADIENT METHOD
                                    !!  * 1: MODIFIED SPECTRAL CONJUGATE GRADIENT METHOD
                                    !!  * 0: STEEPEST DESCENT METHOD

   Optmet = Metopt

END SUBROUTINE ogomet

SUBROUTINE ogopti(Varacc,Numequ,Finish,Desnor,Calval,Calder)

   !! OPTIMISATION PART
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp),intent(inout) :: Varacc !! ITERATION SCALED DISTANCE ACCUMULATED
   integer(ip),intent(in) :: Numequ !! NUMBER OF EQUALITY CONSTRAINTS
   integer(ip),intent(out) :: Finish !! 0=LIMIT 1=OPTIM
   procedure(Calval_f) :: Calval !! FUNCTION FOR VALUES CALDER (VARVAL, CONVAL)
   procedure(Calder_f) :: Calder !! JW : not originally here. added for consistent interface. [not used]

   integer(ip) staflg , faccnt , numcor
   real(wp) Desnor , foldis , cosimp , cornor , quacor , refdis
   real(wp) co0 , co1 , co2 , nor
   real(wp) cosco2 , cosco1
   real(wp) maxdis , norprv

   integer(ip) con , var , cos , act , ind , len , inc
   integer(ip) nnn , typ , des , prv , met
   real(wp) val , max , det , ccc , dis
   real(wp) fac , del , exc , eps , imp
   real(wp) bet , tht
   CHARACTER(len=256) :: str , nam

   real(wp) , DIMENSION(:) , ALLOCATABLE :: cosact
   real(wp) , DIMENSION(:) , ALLOCATABLE :: varvec
   real(wp) , DIMENSION(:) , ALLOCATABLE :: varwrk
   real(wp) , DIMENSION(:) , ALLOCATABLE :: corvec
   real(wp) , DIMENSION(:) , ALLOCATABLE :: desder
   real(wp) , DIMENSION(:) , ALLOCATABLE :: desprv
   real(wp) , DIMENSION(:) , ALLOCATABLE :: varprv
   real(wp) , DIMENSION(:) , ALLOCATABLE :: convec
   real(wp) , DIMENSION(:) , ALLOCATABLE :: conqua
   integer(ip) , DIMENSION(:) , ALLOCATABLE :: concor
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
            Desnor = sqrt(sum(Conred(cos,Numact+1:Numvar)**2))
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
               fac = sqrt(fac)
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
                  val = abs(del)*Varmax
                  IF ( val<eps ) CYCLE
                  fac = dot_product(Conred(con,1:Numvar),Conred(con,1:Numvar))
                  fac = sqrt(fac)
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
            IF ( abs(cosimp)<=1D0 ) THEN
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
               val = sqrt(sum((Varval-Varref+Vardes*dis/Desnor)**2))
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
            Desnor = sqrt(sum(desprv**2))
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
            norprv = sqrt(sum(varprv**2))
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
            nor = sqrt(sum(Vardir**2))
            WRITE (str,'("NOR=",D13.6)') nor
!      CALL OGWRIT (2,STR)
! ----------------------------------------------------------------------
! MET = 3: CONJUGATE GRADIENT METHOD
! MET = 2: SPECTRAL CONJUGATE GRADIENT METHOD
! MET = 1: MODIFIED SPECTRAL CONJUGATE GRADIENT METHOD
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
               IF ( abs(val)>1D-12 .AND. abs(fac)>1D-12 ) THEN
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
            Desnor = sqrt(sum(Vardes**2))
            nor = Desnor
            DO con = 1 , Numcon
               IF ( Contyp(con)==-2 ) CYCLE
               IF ( Conact(con)/=0 ) CYCLE
               del = dot_product(Conder(con,1:Numvar),Vardes(1:Numvar))/nor
               val = abs(del)*Varmax
               IF ( val<eps ) CYCLE
               fac = dot_product(Conder(con,1:Numvar),Conder(con,1:Numvar))
               del = del/sqrt(fac)
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
                  val = abs(del)*Varmax/Desnor
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
            Desnor = sqrt(sum(Vardes**2))
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
            Desnor = sqrt(sum(Vardes**2))
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
               dis = (sqrt(det)-co1)/co2
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
               !CALL ogeval(varvec,convec,0,Conder,Calval,Calval)  ! JW : is this a typo ??
               CALL ogeval(varvec,convec,0,Conder,Calval,Calder)   ! JW : replaced with this
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
            cornor = sqrt(sum(corvec(1:Numact)**2))*0.5D0/Desnor/Desnor
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
            IF ( abs(cosimp)<=1D0 ) THEN
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
                  det = sqrt(det)
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
            dis = sqrt(dis)
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
            ccc = sqrt(sum((Varval-Varref)**2)) - Varmax**2
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

SUBROUTINE ogplog(Luplog,Bosver)

   !! DEFINE WRITING IN PYGMO LOG FORMAT
   !!
   !! 2023/01/25 | W. MARTENS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Luplog !! LOGICAL UNIT FOR WRITING PYGMO LOG
   integer(ip),intent(in) :: Bosver !! VERBOSITY LEVEL:
                                    !!  * 0=NO OUTPUT
                                    !!  * 1 OUTPUT EVERY ITERATION
                                    !!  * 2 OUTPUT EVERY 2ND ITERATION
                                    !!  * N OUTPUT EVERY NTH ITERATION

   Loglup = Luplog
   Verbos = Bosver
   Fevals = 0     ! initialize number of cost fun evaluations
   Pygfla = 0     ! pygmo output status flag: 0: continue iterating, 1: final output

END SUBROUTINE ogplog

SUBROUTINE ogpwri(Objval,Numvio,Convio)

   !! WRITE OPTIMIZATION LOG IN PYGMO FORMAT
   !!
   !! 2023/01/25 | W. MARTENS | NEW

   IMPLICIT NONE

   real(wp),intent(in) :: Objval !! OBJECTIVE VALUE
   integer(ip),intent(in) :: Numvio !! NUMBER OF VIOLATED CONSTRAINTS
   real(wp),intent(in) :: Convio !! TOTAL CONSTRAINT VIOLATION

   CHARACTER(len=2) :: feas
   CHARACTER(len=24) :: fmt

   IF ( Verbos==0 ) RETURN
   ! Print header
   IF ( Fevals==0 ) CALL ogpwri_start()
   ! Increase counter for cost function evaluations
   Fevals = Fevals + 1
   ! Every 50 lines print the column names.
   IF ( mod(real(Fevals-1D0)/real(Verbos),50D0)==0D0 ) &
      WRITE (Loglup,'(A10,A15,A15,A15,A2)') "objevals:" , "objval:" , "violated:" , "viol. norm:"
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

END SUBROUTINE ogpwri

SUBROUTINE ogpwri_end(Objval,Numvio,Convio)

   !! WRITE OPTIMIZATION END RESULT IN PYGMO FORMAT
   !!
   !! 2023/01/25 | W. MARTENS | NEW

   IMPLICIT NONE

   real(wp),intent(in) :: Objval !! OBJECTIVE VALUE
   real(wp),intent(in) :: Convio !! TOTAL CONSTRAINT VIOLATION
   integer(ip),intent(in) :: Numvio !! NUMBER OF VIOLATED CONSTRAINTS

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

END SUBROUTINE ogpwri_end

SUBROUTINE ogpwri_start()

   !! WRITE OPTIMIZATION LOG IN PYGMO FORMAT
   !!
   !! 2023/01/25 | W. MARTENS | NEW

   IMPLICIT NONE

   ! VARDER is the DERIVATIVES COMPUTATION MODE
   !   -> 0: VALUES ONLY
   !   -> 1: USER DEFINED
   !   -> 2: NUMERIC WITH DOUBLE DIFFERENCING
   !   -> 3: NUMERIC WITH SINGLE DIFFERENCING

   WRITE (Loglup,'("OPTGRA plugin for pagmo/pygmo:")')
   IF ( Varder==0 ) THEN
      WRITE (Loglup,'("")')
   ELSEIF ( Varder==1 .OR. Varder==-1 ) THEN
      WRITE (Loglup,'("    User-defined gradients")')
   ELSEIF ( Varder==2 ) THEN
      WRITE (Loglup,'("    Numerical gradients by double differencing")')
   ELSEIF ( Varder==3 ) THEN
      WRITE (Loglup,'("    Numerical gradients by single differencing")')
   ENDIF

   IF ( Optmet==3 ) THEN
      WRITE (Loglup,'("    Conjugate gradient method")')
   ELSEIF ( Optmet==2 ) THEN
      WRITE (Loglup,'("    Spectral conjugate gradient method")')
   ELSEIF ( Optmet==1 ) THEN
      WRITE (Loglup,'("    Modified spectral conjugate gradient method")')
   ELSEIF ( Optmet==0 ) THEN
      WRITE (Loglup,'("    Steepest descent method")')
   ENDIF

   WRITE (Loglup,'("")')

END SUBROUTINE ogpwri_start

SUBROUTINE ogrigt(Actinp,Actout)

   !! RIGHT-MULTIPLIES VECTOR LOWER TRIANGULAR MATRIX OBTAINED BY REDUCTION
   !! AND SUBSEQUENT INVERSION OF DERIVATIVES OF ACTIVE CONSTRAINTS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp),intent(in) :: Actinp(Numcon) !! VECTOR INITAL
   real(wp),intent(out) :: Actout(Numcon) !! VECTOR FINAL (MAY BE SAME AS ACTINP)

   integer(ip) row , col , act
   real(wp) val

   DO col = Numact , 1 , -1
      val = Actinp(col)
      DO act = Numact , col + 1 , -1
         row = Actcon(act)
         val = val - Conred(row,col)*Actout(act)
      ENDDO
      row = Actcon(col)
      Actout(col) = val/Conred(row,col)
   ENDDO

END SUBROUTINE ogrigt

SUBROUTINE ogsens(Consta,Concon,Convar,Varcon,Varvar)

   !! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(out) :: Consta(Numcon) !! CONSTRAINT STATUS (0=PAS 1=ACT)
   real(wp),intent(out) :: Concon(Numcon+1,Numcon) !! SENSITIVITY OF CONTRAINTS+MERIT W.R.T. ACTIVE CONSTRAINTS
   real(wp),intent(out) :: Convar(Numcon+1,Numvar) !! SENSITIVITY OF CONTRAINTS+MERIT W.R.T. PARAMETERS
   real(wp),intent(out) :: Varcon(Numvar,Numcon) !! SENSITIVITY OF VARIABLES W.R.T. ACTIVE CONSTRAINTS
   real(wp),intent(out) :: Varvar(Numvar,Numvar) !! SENSITIVITY OF VARIABLES W.R.T. PARAMETERS
                                                 !! -> NOT SCALED

   real(wp) val , sca
   integer(ip) var , con , act , par , ind , typ

   ! CONVERGED
   Consta = 0
   DO act = 1 , Numact
      con = Actcon(act)
      Consta(con) = 1
   ENDDO

   ! SENSITIVITY OF CONTRAINTS W.R.T. ACTIVE CONSTRAINTS
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

   ! SENSITIVITY OF CONSTRAINTS W.R.T. PARAMETERS
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

   ! SENSITIVITY OF VARIABLES W.R.T. ACTIVE CONSTRAINTS
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

   ! SENSITIVITY OF VARIABLES W.R.T. PARAMETERS
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

   ! DESCALE SENSITIVITY
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

   DO var = 1 , Numvar
      sca = Varsca(var)
      Varcon(var,1:Numcon) = Varcon(var,1:Numcon)*sca
      Varvar(var,1:Numvar) = Varvar(var,1:Numvar)*sca
      Convar(1:Numcon+1,var) = Convar(1:Numcon+1,var)/sca
      Varvar(1:Numvar,var) = Varvar(1:Numvar,var)/sca
   ENDDO

END SUBROUTINE ogsens

SUBROUTINE ogsopt(Optsen)

   !! LINEAR OPTIMISATION MODE
   !!
   !! 2021/03/30 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Optsen !! SENSITIVITY OPTIMISATION MODE
                                    !!
                                    !!  *  0: NO
                                    !!  * -1: INITIALISATION
                                    !!  * +1: WITH CONSTRAINT CALCULATION
                                    !!  * +2: WITH CONSTRAINT BIAS
                                    !!  * +3: WITH CONSTRAINT CALC / NO OPTIM
                                    !!  * +4: WITH CONSTRAINT BIAS / NO OPTIM

   Senopt = Optsen

END SUBROUTINE ogsopt

SUBROUTINE ogssst(Varsen,Quasen,Consen,Actsen,Dersen,Actsav,Consav,Redsav,Dersav,Actnum)

   !! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
   !!
   !! Function to get sensitivity state data, necessary for serialization.
   !! Do not use this directly except in serialization routines
   !!
   !! 2021/07/19 | M. von Looz | NEW

   IMPLICIT NONE

   real(wp),intent(in)    :: Varsen(Numvar) !! STORED VARIABLES VALUE
   real(wp),intent(in)    :: Quasen(Numcon+1) !! STORED CONSTRAINTS CORRECTION VECTOR
   real(wp),intent(in)    :: Consen(Numcon+1) !! STORED CONSTRAINTS VALUE
   integer(ip),intent(in) :: Actsen(Numcon+1) !! STORED CONSTRAINTS ACTIVE
   real(wp),intent(in)    :: Dersen(Numcon+1,Numvar) !! STORED DERIVATIVE
   integer(ip),intent(in) :: Actsav(Numcon+1) !! STORED ACTIVE CONSTRAINTS
   integer(ip),intent(in) :: Consav(Numcon+4) !! STORED ACTIVE CONSTRAINTS
   real(wp),intent(in)    :: Redsav(Numcon+3,Numvar) !! STORED DERIVATIVE
   real(wp),intent(in)    :: Dersav(Numcon+3,Numvar) !! STORED DERIVATIVE

   integer(ip) Actnum
   integer(ip) var , con

   ! Variable values saved for sensitivity
   Numact = Actnum

   DO var = 1 , Numvar
      Senvar(var) = Varsen(var)
   ENDDO

   DO con = 1 , Numcon + 1
      Senqua(con) = Quasen(con)
      Sencon(con) = Consen(con)
      Senact(con) = Actsen(con)
      DO var = 1 , Numvar
         Sender(con,var) = Dersen(con,var)
      ENDDO
   ENDDO

   ! Temporary status saved of which constraints are active
   DO con = 1 , Numcon + 1
      Actcon(con) = Actsav(con)
   ENDDO

   DO con = 1 , Numcon + 4
      Conact(con) = Consav(con)
   ENDDO

   DO con = 1 , Numcon + 3
      DO var = 1 , Numvar
         Conred(con,var) = Redsav(con,var)
      ENDDO
   ENDDO

   DO con = 1 , Numcon
      DO var = 1 , Numvar
         Conder(con,var) = Dersav(con,var)
      ENDDO
   ENDDO
END SUBROUTINE ogssst

SUBROUTINE ogvsca(Scavar)

   !! DEFINE VARIABLE SCALE FACTOR
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   real(wp),intent(in) :: Scavar(Numvar) !! VARIABLES SCALE FACTOR

   integer(ip) var

   DO var = 1 , Numvar
      Varsca(var) = Scavar(var)
   ENDDO

END SUBROUTINE ogvsca

SUBROUTINE ogvstr(Strvar,Lenvar)

   !! DEFINE VARIABLE STRING
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   character(len=maxstr),intent(in) :: Strvar(Numvar) !! VARIABLES NAME STRING
   integer(ip),intent(in) :: Lenvar(Numvar) !! VARIABLES NAME LENGTH

   integer(ip) var , len

   DO var = 1 , Numvar
      len = min(Lenvar(var),maxstr)
      Varstr(var) = Strvar(var)
      Varlen(var) = len
   ENDDO

END SUBROUTINE ogvstr

SUBROUTINE ogvtyp(Typvar)

   !! DEFINE VARIABLE TYPE
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Typvar(Numvar) !! VARIABLES TYPE
                                            !!  *  0=FREE VARIABLE
                                            !!  *  1=PARAMETER FOR SENSITIVITY

   integer(ip) var

   DO var = 1 , Numvar
      Vartyp(var) = Typvar(var)
   ENDDO

END SUBROUTINE ogvtyp

SUBROUTINE ogwlog(Lunlog,Levlog)

   !! DEFINE WRITING IN LOG FILE
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Lunlog !! LOGICAL UNIT FOR WRITING LOG
   integer(ip),intent(in) :: Levlog !! LEVEL OF LOG
                                    !!  * 0=NO OUTPUT
                                    !!  * 1<ALL

   Loglun = Lunlog
   Loglev = Levlog

END SUBROUTINE ogwlog

SUBROUTINE ogwmat(Levmat)

   !! DEFINE WRITING IN MATLAB CONSOLE
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Levmat !! LEVEL OF LOG
                                    !!  * 0=NO OUTPUT
                                    !!  * 1<ALL

   Matlev = Levmat

END SUBROUTINE ogwmat

SUBROUTINE ogwrit(Lev,Str)

   !! Write a meessage to the log.
   !!
   !! 2014/07/29 | J. SCHOENMAEKERS | NEW
   !! 2025/02/09 | J. Williams | added trim()

   IMPLICIT NONE

   CHARACTER(len=*),intent(in) :: Str !! string to print
   integer(ip),intent(in) :: Lev !! only print if `Loglev` is >= this

   IF ( Lev<=Loglev ) THEN
      WRITE (Loglun,'(A)') trim(Str)
      FLUSH (Loglun)
   ENDIF

END SUBROUTINE ogwrit

SUBROUTINE ogwtab(Luntab,Levtab)

   !! DEFINE WRITING IN TABLE FILE
   !!
   !! 2008/01/16 | J. SCHOENMAEKERS | NEW

   IMPLICIT NONE

   integer(ip),intent(in) :: Luntab !! LOGICAL UNIT FOR WRITING TABLE
   integer(ip),intent(in) :: Levtab !! LEVEL OF TAB
                                    !!
                                    !!  * 0=NO OUTPUT
                                    !!  * 1<ALL

   Tablun = Luntab
   Tablev = Levtab

END SUBROUTINE ogwtab

SUBROUTINE sum2v(V1,V2,V,K)
   !! Vector addition.
   !!
   !! `V(1:K) = V1(1:K) + V2(1:K)`
   IMPLICIT NONE
   integer(ip) :: i
   real(wp),intent(in) :: V1(*) , V2(*)
   real(wp),intent(out) :: V(*)
   integer(ip),intent(in) :: K
   DO i = 1 , K
      V(i) = V1(i) + V2(i)
   ENDDO
END SUBROUTINE sum2v

!****************************************************************************************************
   end module optgra_module
!****************************************************************************************************
