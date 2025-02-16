!****************************************************************************************************
!>
!  Near-linear optimisation tool tailored for s/c trajectory design.
!
!  This is a modernization of the original Fortran 77 code from: [pyoptgra](https://github.com/esa/pyoptgra)
!
!### History
!  * Jacob Williams, Feb 2025 : created.

!TODO:
! * remove the STOP statement in one subroutine.
! * see about removing the SPAG DISPATCH loop thing which I don't like (study original code).

module optgra_module

   use iso_fortran_env, ip => int32

   implicit none

   private

#ifdef REAL32
   integer,parameter :: wp = real32   !! Real working precision [4 bytes]
#elif REAL64
   integer,parameter :: wp = real64   !! Real working precision [8 bytes]
#elif REAL128
   integer,parameter :: wp = real128  !! Real working precision [16 bytes]
#else
   integer,parameter :: wp = real64   !! Real working precision if not specified [8 bytes]
#endif

   integer,parameter,public :: optgra_wp = wp   !! Working precision

   integer(ip), PARAMETER :: MAXSTR = 80 !! max string length for var names

   type,public :: optgra
      !! Main class.

      private

      integer(ip) :: NUMVAR = 0
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

      integer(ip) :: NUMCON = 0
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

      integer(ip)  :: OPTMET = 2
      integer(ip)  :: MAXITE = 10
      integer(ip)  :: CORITE = 10
      integer(ip)  :: OPTITE = 10
      integer(ip)  :: DIVITE = 10
      integer(ip)  :: CNVITE = 10
      real(wp)     :: Varmax = 10.0_wp
      real(wp)     :: VARSND = 1.0_wp
      real(wp)     :: VARSTP = 1.0_wp

      integer(ip) :: VARDER = 1
      real(wp),      DIMENSION(:  ), allocatable :: VARPER

      integer(ip) :: LOGLUN = output_unit  !! log file unit
      integer(ip) :: LOGLEV = 1  !! log level

      integer(ip) :: LOGLUP = output_unit  !! pygmo log file unit
      integer(ip) :: VERBOS = 0  !! pygmo verbosity
      integer(ip) :: FEVALS = 0  !! pygmo: number of const fun evals
      integer(ip) :: PYGFLA = 0  !! pygmo: flag indicating status of optimisation
      integer(ip) :: NUMITE = 0  !! number of iterations

      integer(ip) :: MATLEV = 0

      integer(ip) :: TABLUN = output_unit
      integer(ip) :: TABLEV = 0

      integer(ip) :: SENOPT = 0

      integer(ip) :: NUMACT = 0
      integer(ip),   DIMENSION(:  ), allocatable :: ACTCON
      integer(ip),   DIMENSION(:  ), allocatable :: CONFIX
      integer(ip),   DIMENSION(:  ), allocatable :: CONACT
      real(wp),      DIMENSION(:,:), allocatable :: CONDER
      real(wp),      DIMENSION(:,:), allocatable :: CONRED
      real(wp),      DIMENSION(:,:), allocatable :: SENDER
      !integer(ip) :: CONVER = 0   ! not used ?
      !integer(ip),   DIMENSION(:  ), allocatable :: CONOPT  ! not used ?

      procedure(calval_f),pointer :: calval => null() !! function for values
      procedure(calder_f),pointer :: calder => null() !! function for derivatives

   contains

      private

      procedure,public :: initialize
      procedure,public :: destroy => ogclos
      procedure,public :: solve => ogexec
      procedure :: ogrigt
      procedure :: ogwrit
      procedure :: ogincl
      procedure :: ogexcl
      procedure :: ogeval
      procedure :: ogcorr
      procedure :: ogopti
      procedure :: ogleft
      procedure :: ogpwri
      procedure :: ogpwri_start
      procedure :: ogpwri_end

   end type optgra

   abstract interface
      subroutine calval_f(me,varvec,Valcon,i)
         !! FUNCTION FOR VALUES
         !! INPUT AND OUTPUT NOT SCALED
         import :: optgra,wp
         implicit none
         class(optgra),intent(inout) :: me
         real(wp), dimension(:), intent(in) :: varvec !! size is Numvar
         real(wp), dimension(:), intent(out) :: Valcon !! size is Numcon+1
         integer, intent(in) :: i !! JW: THIS IS NOT DOCUMENTED ?
      end subroutine calval_f
      subroutine calder_f(me,varvec,convec,Dercon)
         !! FUNCTION FOR VALUES AND DERIVATIVES
         !! INPUT AND OUTPUT NOT SCALED
         import :: optgra,wp
         implicit none
         class(optgra),intent(inout) :: me
         real(wp), dimension(:), intent(in) :: varvec !! size is Numvar
         real(wp), dimension(:), intent(out) :: convec !! size is Numcon+1
         real(wp), dimension(:,:), intent(out) :: Dercon !! size is Numcon+1,Numvar
      end subroutine calder_f

   end interface



contains
!****************************************************************************************************

   SUBROUTINE mul2m(A1,M1,K1,L1,N1,A2,M2,K2,L2,N2,A,M,K,L,N)

      !! Matrix multiply.
      !!
      !! `A(K:K+N1,L:L+N) = A1(K1:K1+N1,L1:L1+N2) * A2(K2:K2+N2,L2:L2+N3)`

      integer,intent(in) :: m1, m2, m, k, k1, k2, l, l1 , l2 , n , n1 , n2
      real(wp),intent(out) :: A(M,*)
      real(wp),intent(in) :: A1(M1,*)
      real(wp),intent(in) :: A2(M2,*)

      real(wp) :: f1 , f2
      integer(ip) :: i , i1 , i2 , ic , ir

      DO i1 = K , K + N1 - 1
         DO i = L , L + N - 1
            A(i1,i) = 0.0_wp
         ENDDO
      ENDDO

      DO i1 = 0 , N1 - 1
         DO i2 = 0 , N2 - 1
            IF ( K1>=0 ) THEN
               f1 = A1(i1+K1,i2+L1)
            ELSE
               f1 = A1(i2-K1,i1+L1)
            ENDIF
            IF ( f1/=0.0_wp ) THEN
               DO i = 0 , N - 1
                  IF ( K2>=0 ) THEN
                     f2 = A2(i2+K2,i+L2)
                  ELSE
                     f2 = A2(i-K2,i2+L2)
                  ENDIF
                  IF ( f2/=0.0_wp ) THEN
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

      real(wp),intent(in) :: A !! SCALAR
      real(wp),intent(in) :: X(*) !! VECTOR
      real(wp),intent(out) :: Z(*) !! VECTOR
      integer(ip),intent(in) :: Kd !! NUMBER OF ELEMENTS TO BE USED
      integer(ip) :: i

      DO i = 1 , Kd
         Z(i) = X(i)*A
      ENDDO
   END SUBROUTINE mulvs

   SUBROUTINE sum2v(V1,V2,V,K)
      !! Vector addition.
      !!
      !! `V(1:K) = V1(1:K) + V2(1:K)`

      real(wp),intent(in) :: V1(*) , V2(*)
      real(wp),intent(out) :: V(*)
      integer(ip),intent(in) :: K
      integer(ip) :: i

      DO i = 1 , K
         V(i) = V1(i) + V2(i)
      ENDDO
   END SUBROUTINE sum2v

   subroutine initialize(me,Numvar,Numcon,Calval,Calder,Delcon,Conpri,Consca,Constr,Conlen,&
      Contyp,Varder,Varper,Varmax,Varsnd,&
      Maxite,Itecor,Iteopt,Itediv,Itecnv,&
      Loglup,Verbos,Senopt,Varsca,&
      Varstr,Varlen,Vartyp,Loglun,Loglev,Matlev,&
      Tablun,Tablev,Optmet)

      !! Initialize the class. This should be the first routine called.

      class(optgra),intent(out) :: me
      integer(ip),intent(in) :: Numvar !! NUMBER OF VARIABLES
      integer(ip),intent(in) :: Numcon !! NUMBER OF CONSTRAINTS
      procedure(Calval_f) :: Calval !! FUNCTION FOR VALUES
                                    !! -> INPUT AND OUTPUT NOT SCALED
      procedure(Calder_f) :: Calder !! FUNCTION FOR VALUES AND DERIVATIVES
                                    !! -> INPUT AND OUTPUT NOT SCALED
      real(wp),intent(in) :: Delcon(Numcon+1)  !! CONSTRAINTS DELTAS
                                               !! (CONSTRAINT + MERIT CONVERGENCE THRESHOLDS)
      integer(ip),intent(in) :: Conpri(Numcon+1)   !! CONSTRAINTS PRIORITY (1:NUMCON)
                                                   !! -> 1-N
      real(wp),intent(in) :: Consca(Numcon+1)  !! CONSTRAINTS CONVER THRESHOLD (1:NUMCON)
                                               !! MERIT       CONVER THRESHOLD (1+NUMCON)
      character(len=maxstr),intent(in) :: Constr(Numcon+1) !! CONIABLES NAME STRING
      integer(ip),intent(in) :: Conlen(Numcon+1) !! CONIABLES NAME LENGTH
      integer(ip),intent(in) :: Contyp(Numcon+1) !! CONSTRAINTS TYPE (1:NUMCON)
                                                 !!
                                                 !!  *  1=GTE
                                                 !!  * -1=LTE
                                                 !!  *  0=EQU
                                                 !!  * -2=DERIVED DATA
                                                 !!
                                                 !! MERIT TYPE (1+NUMCON)
                                                 !!
                                                 !!  * 1=MAX
                                                 !!  * -1=MIN
      integer(ip),intent(in) :: Varder  !! DERIVATIVES COMPUTATION MODE
                                        !!
                                        !!  * 1: USER DEFINED
                                        !!  * 2: NUMERIC WITH DOUBLE DIFFERENCING
                                        !!  * 3: NUMERIC WITH SINGLE DIFFERENCING
      real(wp),intent(in) :: Varper(Numvar) !! VARIABLES PERTURBATION FOR DERIVATIVES
                                            !! -> NOT SCALED
      real(wp),intent(in) :: Varmax !! MAXIMUM DISTANCE PER ITERATION
                                    !!  -> SCALED
      real(wp),intent(in) :: Varsnd !! PERTURBATION FOR 2ND ORDER DERIVATIVES
                                    !!  -> SCALED
      integer(ip),intent(in) :: Maxite !! MAXIMUM NUMBER OF ITERATIONS
      integer(ip),intent(in) :: Itecor
      integer(ip),intent(in) :: Iteopt
      integer(ip),intent(in) :: Itediv
      integer(ip),intent(in) :: Itecnv
      integer(ip),intent(in) :: Loglup !! LOGICAL UNIT FOR WRITING PYGMO LOG
      integer(ip),intent(in) :: Verbos !! VERBOSITY LEVEL:
                                       !!
                                       !!  * 0=NO OUTPUT
                                       !!  * 1 OUTPUT EVERY ITERATION
                                       !!  * 2 OUTPUT EVERY 2ND ITERATION
                                       !!  * N OUTPUT EVERY NTH ITERATION
      integer(ip),intent(in) :: Senopt  !! SENSITIVITY OPTIMISATION MODE
                                        !!
                                        !!  *  0: NO
                                        !!  * -1: INITIALISATION
                                        !!  * +1: WITH CONSTRAINT CALCULATION
                                        !!  * +2: WITH CONSTRAINT BIAS
                                        !!  * +3: WITH CONSTRAINT CALC / NO OPTIM
                                        !!  * +4: WITH CONSTRAINT BIAS / NO OPTIM

      real(wp),intent(in) :: Varsca(Numvar) !! VARIABLES SCALE FACTOR
      character(len=maxstr),intent(in) :: Varstr(Numvar) !! VARIABLES NAME STRING
      integer(ip),intent(in) :: Varlen(Numvar) !! VARIABLES NAME LENGTH
      integer(ip),intent(in) :: Vartyp(Numvar) !! VARIABLES TYPE
                                               !!
                                               !!  *  0=FREE VARIABLE
                                               !!  *  1=PARAMETER FOR SENSITIVITY
      integer(ip),intent(in) :: Loglun !! LOGICAL UNIT FOR WRITING LOG
      integer(ip),intent(in) :: Loglev !! LEVEL OF LOG:
                                       !!
                                       !!  * 0=NO OUTPUT
                                       !!  * 1<ALL
      integer(ip),intent(in) :: Matlev !! LEVEL OF LOG:
                                       !!
                                       !!  * 0=NO OUTPUT
                                       !!  * 1<ALL
      integer(ip),intent(in) :: Tablun !! LOGICAL UNIT FOR WRITING TABLE
      integer(ip),intent(in) :: Tablev !! LEVEL OF TAB
                                       !!
                                       !!  * 0=NO OUTPUT
                                       !!  * 1<ALL
      integer(ip),intent(in) :: Optmet !! OPTIMISATION METHOD:
                                       !!
                                       !!  * 3: CONJUGATE GRADIENT METHOD
                                       !!  * 2: SPECTRAL CONJUGATE GRADIENT METHOD
                                       !!  * 1: MODIFIED SPECTRAL CONJUGATE GRADIENT METHOD
                                       !!  * 0: STEEPEST DESCENT METHOD

      integer(ip) :: con, var !! counter

      call me%destroy()

      ! set the functions:
      me%Calval => Calval
      me%Calder => Calder

      ! VARIABLES
      me%Numvar = Numvar

      ALLOCATE (me%Varval(me%Numvar))
      ALLOCATE (me%Vartyp(me%Numvar))
      ALLOCATE (me%Varsca(me%Numvar))
      ALLOCATE (me%Varstr(me%Numvar))
      ALLOCATE (me%Varlen(me%Numvar))
      ALLOCATE (me%Varref(me%Numvar))
      ALLOCATE (me%Vardes(me%Numvar))
      ALLOCATE (me%Vargrd(me%Numvar))
      ALLOCATE (me%Vardir(me%Numvar))
      ALLOCATE (me%Funvar(me%Numvar))
      ALLOCATE (me%Senvar(me%Numvar))

      DO var = 1 , me%Numvar
         me%Varval(var) = 0.0_wp
         me%Vartyp(var) = 0
         me%Varsca(var) = 1.0_wp
         me%Varstr(var) = ""
         me%Varlen(var) = 0
         me%Varref(var) = 0.0_wp
         me%Vardes(var) = 0.0_wp
         me%Vargrd(var) = 0.0_wp
         me%Vardir(var) = 0.0_wp
         me%Funvar(var) = 0.0_wp
         me%Senvar(var) = 0.0_wp
      ENDDO

      ! CONSTRAINTS
      me%Numcon = Numcon

      ALLOCATE (me%Conval(me%Numcon+1))
      ALLOCATE (me%Contyp(me%Numcon+1))
      ALLOCATE (me%Conpri(me%Numcon+1))
      ALLOCATE (me%Consca(me%Numcon+1))
      ALLOCATE (me%Constr(me%Numcon+1))
      ALLOCATE (me%Conlen(me%Numcon+1))
      ALLOCATE (me%Conref(me%Numcon+1))
      ALLOCATE (me%Senqua(me%Numcon+1))
      ALLOCATE (me%Sencon(me%Numcon+1))
      ALLOCATE (me%Sendel(me%Numcon+1))
      ALLOCATE (me%Senact(me%Numcon+1))

      DO con = 1 , me%Numcon + 1
         me%Conval(con) = 0.0_wp
         me%Contyp(con) = 0
         me%Conpri(con) = 1
         me%Consca(con) = 1.0_wp
         me%Constr(con) = ""
         me%Conlen(con) = 0
         me%Conref(con) = 0.0_wp
         me%Senqua(con) = 0.0_wp
         me%Sencon(con) = 0.0_wp
         me%Sendel(con) = 0.0_wp
         me%Senact(con) = 0
      ENDDO

      ! CONTROL
      me%Optmet = 2
      me%Maxite = 10
      me%Corite = 10
      me%Optite = 10
      me%Divite = 10
      me%Cnvite = 10
      me%Varmax = 10.0_wp
      me%Varsnd = 1.0_wp

      ! DERIVATIVES
      me%Varder = 1
      ALLOCATE (me%Varper(me%Numvar))
      DO var = 1 , me%Numvar
         me%Varper(var) = 1.0e-03_wp
      ENDDO

      me%Loglun = output_unit     ! LOG FILE
      me%Loglev = 1
      me%Loglup = output_unit     ! PYGMO LOG FILE
      me%Loglev = 0
      me%Matlev = 0               ! MATLAB CONSOLE
      me%Tablun = output_unit     ! TABLE FILE
      me%Tablev = 0

      me%Senopt = 0    ! LINEAR OPTIMISATION MODE

      ! WORKING VECTORS
      ALLOCATE (me%Actcon(me%Numcon+1))
      ALLOCATE (me%Confix(me%Numcon))
      ALLOCATE (me%Conact(me%Numcon+4))
      ALLOCATE (me%Conder(me%Numcon+3,me%Numvar))
      ALLOCATE (me%Conred(me%Numcon+3,me%Numvar))
      ALLOCATE (me%Sender(me%Numcon+3,me%Numvar))
      ! ALLOCATE (me%Conopt(me%Numcon+1))
      me%Numact = 0
      me%Actcon = 0
      me%Conact = 0
      me%Confix = 0
      me%Conder = 0.0_wp
      me%Conred = 0.0_wp
      !me%Conopt = 0

      ! NOTE: should the last element also be set? (not done in original)  ?????
      DO con = 1 , me%Numcon
         me%Sendel(con) = Delcon(con)
      ENDDO

      DO con = 1 , me%Numcon + 1
         me%Conpri(con) = Conpri(con)
         me%Consca(con) = Consca(con)
         me%Constr(con) = Constr(con)
         me%Conlen(con) = min(Conlen(con),maxstr)
         me%Contyp(con) = Contyp(con)
      ENDDO

      me%Varder = Varder

      DO var = 1 , me%Numvar
         me%Varper(var) = Varper(var)
         me%Varsca(var) = Varsca(var)
         me%Varstr(var) = Varstr(var)
         me%Varlen(var) = min(Varlen(var),maxstr)
         me%Vartyp(var) = Vartyp(var)
      ENDDO

      me%Varmax = Varmax
      me%Varsnd = Varsnd

      me%Maxite = Maxite
      me%Corite = Itecor
      me%Optite = Iteopt
      me%Divite = Itediv
      me%Cnvite = Itecnv
      IF ( me%Corite>me%Maxite ) me%Corite = me%Maxite
      IF ( me%Optite>me%Maxite ) me%Optite = me%Maxite
      IF ( me%Divite>me%Corite ) me%Divite = me%Corite
      IF ( me%Cnvite>me%Optite ) me%Cnvite = me%Optite

      me%Loglup = Loglup
      me%Verbos = Verbos
      me%Fevals = 0     ! initialize number of cost fun evaluations
      me%Pygfla = 0     ! pygmo output status flag: 0: continue iterating, 1: final output

      me%Senopt = Senopt

      me%Loglun = Loglun
      me%Loglev = Loglev
      me%Matlev = Matlev
      me%Tablun = Tablun
      me%Tablev = Tablev

      me%Optmet = Optmet

   end subroutine initialize

   SUBROUTINE ogclos(me)
      !! DEALLOCATION OF ARRAYS
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me

      ! VARIABLES
      if (allocated(me%Varval)) DEALLOCATE (me%Varval)
      if (allocated(me%Vartyp)) DEALLOCATE (me%Vartyp)
      if (allocated(me%Varsca)) DEALLOCATE (me%Varsca)
      if (allocated(me%Varstr)) DEALLOCATE (me%Varstr)
      if (allocated(me%Varlen)) DEALLOCATE (me%Varlen)
      if (allocated(me%Varref)) DEALLOCATE (me%Varref)
      if (allocated(me%Vardes)) DEALLOCATE (me%Vardes)
      if (allocated(me%Vargrd)) DEALLOCATE (me%Vargrd)
      if (allocated(me%Vardir)) DEALLOCATE (me%Vardir)
      if (allocated(me%Funvar)) DEALLOCATE (me%Funvar)
      if (allocated(me%Senvar)) DEALLOCATE (me%Senvar)

      ! CONSTRAINTS
      if (allocated(me%Conval)) DEALLOCATE (me%Conval)
      if (allocated(me%Contyp)) DEALLOCATE (me%Contyp)
      if (allocated(me%Conpri)) DEALLOCATE (me%Conpri)
      if (allocated(me%Consca)) DEALLOCATE (me%Consca)
      if (allocated(me%Constr)) DEALLOCATE (me%Constr)
      if (allocated(me%Conlen)) DEALLOCATE (me%Conlen)
      if (allocated(me%Conref)) DEALLOCATE (me%Conref)
      if (allocated(me%Senqua)) DEALLOCATE (me%Senqua)
      if (allocated(me%Sencon)) DEALLOCATE (me%Sencon)
      if (allocated(me%Sendel)) DEALLOCATE (me%Sendel)
      if (allocated(me%Senact)) DEALLOCATE (me%Senact)

      ! DERIVATIVES
      if (allocated(me%Varper)) DEALLOCATE (me%Varper)

      ! WORKING VECTORS
      if (allocated(me%Actcon)) DEALLOCATE (me%Actcon)
      if (allocated(me%Confix)) DEALLOCATE (me%Confix)
      if (allocated(me%Conact)) DEALLOCATE (me%Conact)
      if (allocated(me%Conder)) DEALLOCATE (me%Conder)
      if (allocated(me%Conred)) DEALLOCATE (me%Conred)
      if (allocated(me%Sender)) DEALLOCATE (me%Sender)
      !if (allocated(me%Conopt)) DEALLOCATE (me%Conopt)

   END SUBROUTINE ogclos

   SUBROUTINE ogcorr(me,Varacc,Finish,Toterr,Norerr)

      !! CORRECTION PART
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp) :: Varacc
      integer(ip) :: Finish
      real(wp) :: Toterr
      real(wp) :: Norerr

      integer(ip) :: coritr , numfff , minpri , maxpri , curpri
      real(wp) :: cornor , foldis , cstval , conerr
      real(wp) :: corinv , varvio , conmax , normax
      integer(ip) :: conind , norind , inelop , maxitr
      integer(ip) :: con , var , act , ind , len , cos , stp
      integer(ip) :: typ , cor , pri , vio , fff
      real(wp) :: val , fac , upr , del , co2 , co1 , co0 , de2 , dis
      real(wp) :: eps , err , dlt , sca , dif
      real(wp) :: exc
      CHARACTER(len=256) :: str , nam
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
            ALLOCATE (cosact(me%Numvar))
            ALLOCATE (varvec(me%Numvar))
            ALLOCATE (varsav(me%Numvar))
            ALLOCATE (varcor(me%Numvar))
            ALLOCATE (corvec(me%Numvar))
            ALLOCATE (consav(me%Numcon+1))
            ALLOCATE (conttt(me%Numcon+1))
            ALLOCATE (concor(me%Numcon))
            ALLOCATE (coninc(me%Numcon))
            ALLOCATE (conhit(me%Numcon))
            ALLOCATE (fffcon(me%Numcon))
            ALLOCATE (prisav(me%Numcon))
            ! ======================================================================
            ! CORRECTION PART
            ! ----------------------------------------------------------------------
            coritr = 0
            maxitr = 10
            IF ( me%Senopt/=0 ) maxitr = 1
            ! ======================================================================
            cos = me%Numcon + 1
            vio = me%Numcon + 2
            stp = 0
            Finish = 0
            eps = 1.0e-03_wp
            dlt = 1.0e-06_wp
            varvio = me%Varmax*1.0e+03_wp
            numfff = me%Numact
            fffcon = me%Actcon
            conttt = me%Contyp
            prisav = me%Conpri
            concor = 0
            ! ----------------------------------------------------------------------
            minpri = 1000
            maxpri = -1000
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               minpri = min(minpri,me%Conpri(con))
               maxpri = max(maxpri,me%Conpri(con))
            ENDDO
            ind = 0
            IF ( numfff>0 ) ind = 1
            IF ( me%Senopt>0 ) ind = 1
            minpri = minpri - ind
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3,"PRIORITISE CONSTRAINTS")
            CALL me%ogwrit(3,"")
            IF ( me%Senopt<=0 ) THEN
               DO fff = 1 , numfff
                  con = fffcon(fff)
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  typ = me%Contyp(con)
                  WRITE (str,'("PRI",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
                  CALL me%ogwrit(3,str)
                  me%Conpri(con) = minpri
               ENDDO
            ENDIF
            IF ( me%Senopt>0 ) THEN
               DO con = 1 , me%Numcon
                  IF ( me%Contyp(con)==-2 ) CYCLE
                  IF ( me%Senact(con)<=0 ) CYCLE
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  typ = me%Contyp(con)
                  WRITE (str,'("PRI",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
                  CALL me%ogwrit(3,str)
                  me%Conpri(con) = minpri
               ENDDO
            ENDIF
            CALL me%ogwrit(3,"")
            spag_nextblock_1 = 2
          CASE (2)
            ! ======================================================================
            ! Evaluation loop
            ! ----------------------------------------------------------------------
            WRITE (str,'("CORITR=",I2,1X,I2)') coritr , maxitr
            CALL me%ogwrit(3,str)
            CALL me%ogwrit(3,"")
            ! ======================================================================
            ! Inequality loop
            ! ----------------------------------------------------------------------
            inelop = 2
            IF ( numfff>0 ) inelop = 1
            IF ( me%Senopt>0 ) inelop = 1
            varsav = me%Varval
            consav = me%Conval
            spag_nextblock_1 = 3
          CASE (3)
            ! ----------------------------------------------------------------------
            WRITE (str,'("INELOP=",I2)') inelop
            CALL me%ogwrit(3,str)
            CALL me%ogwrit(3,"")
            ! ----------------------------------------------------------------------
            me%Varval = varsav
            me%Conval = consav
            me%Contyp = conttt
            IF ( inelop==1 ) THEN
               DO fff = 1 , numfff
                  con = fffcon(fff)
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  typ = me%Contyp(con)
                  WRITE (str,'("TYP",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
                  CALL me%ogwrit(3,str)
                  me%Contyp(con) = 0
               ENDDO
            ENDIF
            IF ( inelop==1 .AND.me%Senopt>0 ) THEN
               DO con = 1 , me%Numcon
                  IF ( me%Contyp(con)==-2 ) CYCLE
                  IF ( me%Senact(con)<=0 ) CYCLE
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  typ = me%Contyp(con)
                  WRITE (str,'("FIX",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
                  CALL me%ogwrit(2,str)
                  me%Contyp(con) = 0
               ENDDO
            ENDIF
            me%Numact = 0
            me%Conact = 0
            me%Conred(1:me%Numcon+1,:) = me%Conder(1:me%Numcon+1,:)
            me%Conred(me%Numcon+2,:) = me%Vardir
            ! ----------------------------------------------------------------------
            ! CHECK CONSTRAINTS
            ! ----------------------------------------------------------------------
            conerr = 0.0_wp
            Norerr = 0.0_wp
            conind = 0
            norind = 0
            conmax = 0.0_wp
            normax = 0.0_wp
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               typ = me%Contyp(con)
               val = me%Conval(con)
               IF ( val<-1.0_wp ) THEN
                  err = abs(val)
               ELSEIF ( typ/=0 ) THEN
                  err = 0.0_wp
               ELSEIF ( val>1.0_wp ) THEN
                  err = abs(val)
               ELSE
                  err = 0.0_wp
               ENDIF
               conerr = conerr + err
               IF ( err>conmax ) THEN
                  conind = con
                  conmax = err
               ENDIF
               fac = 0.0_wp
               DO var = 1 , me%Numvar
                  IF ( me%Vartyp(var)==1 ) CYCLE
                  fac = fac + me%Conder(con,var)**2
               ENDDO
               fac = sqrt(fac)
               IF ( err==0.0_wp ) THEN
               ELSEIF ( fac/=0.0_wp ) THEN
                  err = err/fac
               ELSE
                  CALL me%ogwrit(0,"")
                  CALL me%ogwrit(0,"ERROR: CONSTRAINT CAN NOT BE SATISFIED")
                  WRITE (str,'("CON/VAL= ",I5,1X,D13.6)') con , val
                  CALL me%ogwrit(0,str)
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
            CALL me%ogwrit(3,"")
            WRITE (str,'("NUMFFF/TOTERR/NORERR/COSVAL=",I4,3(1X,D13.6))') numfff , Toterr , Norerr , me%Conval(me%Numcon+1)
            CALL me%ogwrit(2,str)
            CALL me%ogwrit(3,"")
            WRITE (str,'("MAXIM TOTAL ERROR.: ",D13.6,I6)') conmax , conind
            CALL me%ogwrit(3,str)
            WRITE (str,'("MAXIM NORM  ERROR.: ",D13.6,I6)') normax , norind
            CALL me%ogwrit(3,str)
            ! ----------------------------------------------------------------------
            IF ( me%Senopt>0 .AND. coritr==maxitr ) THEN
               IF ( cstval==0.0_wp ) THEN
                  Finish = 1
               ELSE
                  Finish = 0
               ENDIF
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( coritr==0 .AND. conerr==0.0_wp ) THEN
               me%Numact = numfff
               me%Actcon = fffcon
               Finish = 1
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( coritr/=0 .AND. conerr==0.0_wp ) THEN
               Finish = 1
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( coritr==maxitr ) THEN
               Finish = 0
               CALL me%ogwrit(3,"")
               WRITE (str,'("CORITR=",I2)') coritr
               CALL me%ogwrit(3,str)
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
            CALL me%ogwrit(3,str)
            CALL me%ogwrit(3,"")
            ! ======================================================================
            ! MINIMUM NORM CORRECTION
            ! ----------------------------------------------------------------------
            CALL me%ogwrit(3,"")
            WRITE (str,'("CORRECTION OF CONSTRAINTS")')
            CALL me%ogwrit(3,str)
            CALL me%ogwrit(3,"")
            WRITE (str,*) "INELOP/CURPRI=" , inelop , curpri
            CALL me%ogwrit(3,str)
            ! ----------------------------------------------------------------------
            conerr = 0.0_wp
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               conhit(con) = 0
               IF ( me%Conval(con)<-eps ) THEN
                  conerr = conerr + abs(me%Conval(con))
               ELSEIF ( me%Contyp(con)/=0 ) THEN
               ELSEIF ( me%Conval(con)>+eps ) THEN
                  conerr = conerr + abs(me%Conval(con))
               ENDIF
            ENDDO
            ! ----------------------------------------------------------------------
            CALL me%ogwrit(3,"")
            WRITE (str,'("LINEAR ERROR.: ",D13.6)') conerr
            CALL me%ogwrit(3,str)
            ! ----------------------------------------------------------------------
            CALL me%ogwrit(3,"")
            WRITE (str,'(" ACT  PAS  MOV",'//' " COST___VAL COST___GRD",'//' " DIST___DEL CONSTRAINT")')
            CALL me%ogwrit(3,str)
            spag_nextblock_1 = 5
          CASE (5)
            ! ======================================================================
            ! Move loop
            ! ----------------------------------------------------------------------
            ! ----------------------------------------------------------------------
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               coninc(con) = 0
               pri = me%Conpri(con)
               val = me%Conval(con)
               typ = me%Contyp(con)
               act = me%Conact(con)
               cor = concor(con)
               IF ( act>0 ) CYCLE
               IF ( pri>curpri ) THEN
                  me%Conact(con) = -1
               ELSEIF ( val<-dlt ) THEN
                  me%Conact(con) = -1
               ELSEIF ( val>+dlt ) THEN
                  me%Conact(con) = -1
               ELSE
                  me%Conact(con) = 0
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
!          IF (ACT /= CONACT(CON) .OR. COR /= CONCOR(CON)) THEN
!              NAM = CONSTR(CON)
!              LEN = CONLEN(CON)
!              WRITE (STR,'(5X,5X,I4,23X,D10.3,1X,A,4I4)')
!     &        CON, CONVAL(CON),NAM(1:LEN),
!     &        CONACT(CON),CONCOR(CON), ACT, COR
!              CALL me%ogwrit (3,STR)
!          ENDIF
            ENDDO
            ! ======================================================================
            ! STEEPEST ASCENT VECTOR
            ! ======================================================================
            ! MERIT VALUE AND DERIVATIVES
            ! ----------------------------------------------------------------------
            cstval = 0.0_wp
            varvec = 0.0_wp
            me%Conred(vio,:) = 0.0_wp
            ! ----------------------------------------------------------------------
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               IF ( me%Conpri(con)>curpri ) CYCLE
               IF ( concor(con)==0 ) CYCLE
               fac = concor(con)
               cstval = cstval - me%Conval(con)*fac
               me%Conred(vio,:) = me%Conred(vio,:) - me%Conred(con,:)*fac
               varvec = varvec - me%Conder(con,:)*fac
            ENDDO
            SPAG_Loop_1_1: DO
               ! ----------------------------------------------------------------------
               ! ----------------------------------------------------------------------
               ! STEEPEST ASCENT VECTOR
               ! ----------------------------------------------------------------------
               corvec = me%Conred(vio,:)
               ! ----------------------------------------------------------------------
               cornor = sqrt(sum(corvec(me%Numact+1:me%Numvar)**2))
               ! ----------------------------------------------------------------------
               ! MERIT PARTIAL W.R.T. CONSTRAINTS
               ! ----------------------------------------------------------------------
               CALL me%ogrigt(corvec,cosact)
               ! ----------------------------------------------------------------------
               ! CONSTRAINT REMOVAL
               ! ----------------------------------------------------------------------
               ind = 0
               exc = 1.0e-12_wp
               upr = exc
               DO act = 1 , me%Numact
                  con = me%Actcon(act)
                  IF ( me%Contyp(con)==0 ) CYCLE
                  val = cosact(act)
                  IF ( val<=exc ) CYCLE
                  IF ( val<upr ) CYCLE
!                 IF (VAL >= UPR .AND. UPR > 0.0_wp) CYCLE
                  upr = val
                  ind = act
               ENDDO
               ! ----------------------------------------------------------------------
               IF ( ind/=0 ) THEN
                  con = me%Actcon(ind)
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  WRITE (str,'(5X,I4,5X,3(1X,D10.3),1X,A)') con , cstval , cornor , upr , nam(1:len)
                  CALL me%ogwrit(3,str)
                  CALL me%ogexcl(ind)
                  IF ( coninc(con)>=5 ) THEN
                     WRITE (str,'("OGCORR-WARNING: CONSTRAINT INCLUDED")')
                     CALL me%ogwrit(1,str)
                     WRITE (str,'("CON/INC/UPR=",2I4,1X,D10.3)') con , coninc(con) , upr
                     CALL me%ogwrit(1,str)
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
               corinv = 1.0_wp/cornor
               corvec(me%Numact+1:me%Numvar) = corvec(me%Numact+1:me%Numvar)*corinv
               ! ----------------------------------------------------------------------
               ! CONSTRAINT INCLUSION
               ! ----------------------------------------------------------------------
               ind = 0
               upr = 0.0_wp
               ! ----------------------------------------------------------------------
               DO con = 1 , me%Numcon
                  IF ( me%Contyp(con)==-2 ) CYCLE
                  IF ( me%Conpri(con)>curpri ) CYCLE
                  IF ( me%Conact(con)/=0 ) CYCLE
                  del = dot_product(me%Conred(con,me%Numact+1:me%Numvar),me%Conred(vio,me%Numact+1:me%Numvar))
                  val = abs(del)*me%Varmax/cornor
                  IF ( val<eps ) CYCLE
                  fac = dot_product(me%Conred(con,1:me%Numvar),me%Conred(con,1:me%Numvar))
                  del = del/sqrt(fac)
                  IF ( del<upr ) THEN
                     upr = del
                     ind = con
                  ENDIF
                  IF ( me%Contyp(con)/=0 ) CYCLE
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
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  WRITE (str,'(I4,5X,5X,3(1X,D10.3),1X,A)') con , cstval , cornor , upr , nam(1:len)
                  CALL me%ogwrit(3,str)
                  CALL me%ogincl(ind)
                  coninc(con) = coninc(con) + 1
                  CYCLE
               ENDIF
               EXIT SPAG_Loop_1_1
            ENDDO SPAG_Loop_1_1
            ! ----------------------------------------------------------------------
            ! ----------------------------------------------------------------------
            DO var = 1 , me%Numvar
               val = varvec(var)
               DO act = 1 , me%Numact
                  con = me%Actcon(act)
                  val = val - me%Conder(con,var)*cosact(act)
               ENDDO
               varvec(var) = val*corinv
            ENDDO
            ! ----------------------------------------------------------------------
            varcor = me%Varval - me%Varref
            co2 = dot_product(varvec,varvec)
            co1 = dot_product(varvec,varcor)*0.5_wp
            co0 = dot_product(varcor,varcor) - me%Varmax**2
            de2 = co1**2 - co2*co0
            IF ( de2>=0.0_wp .AND. co2/=0.0_wp ) THEN
               dis = (sqrt(de2)-co1)/co2
            ELSE
               dis = 0.0_wp
            ENDIF
            ! ----------------------------------------------------------------------
            DO var = 1 , me%Numvar
               fac = varvec(var)
               IF ( fac==0.0_wp ) CYCLE
               dif = me%Varval(var) - me%Varref(var)
               sca = me%Varmax*1.0e-0_wp   ! JW : why is this multiplied by 1 ?
               val = (dif+sca)/fac
               fac = (dif-sca)/fac
               IF ( fac>val ) val = fac
               IF ( val<dis ) dis = val
            ENDDO
            IF ( dis<0.0_wp ) dis = 0.0_wp
            ! ----------------------------------------------------------------------
            foldis = dis
            ! ======================================================================
            ! OPTIMISE DIRETION OF STEPPEST ASCENT
            ! ======================================================================
            IF ( cstval==0.0_wp ) THEN
               ! ----------------------------------------------------------------------
!              WRITE (STR,'("CNV=",3(1X,D10.3))') CSTVAL/CORNOR/VARVIO
!              CALL me%ogwrit (3,STR)
               WRITE (str,'(4X,5X,5X,3(1X,D10.3))') cstval , cornor , foldis
               CALL me%ogwrit(3,str)
               CALL me%ogwrit(3,"")
               IF ( curpri>=maxpri ) THEN
                  WRITE (str,'("MAXPRI=",I3)') maxpri
                  CALL me%ogwrit(3,str)
                  ! ======================================================================
                  ! ======================================================================
                  ! MATCHED INEQUALITY CONSTRAINTS + MINIMUM CORRECTION NORM
                  ! ----------------------------------------------------------------------
                  CALL me%ogwrit(3,"")
                  CALL me%ogwrit(3,"STATUS OF CONSTRAINTS:")
                  CALL me%ogwrit(3,"")
                  CALL me%ogwrit(3," ACT  PAS  NON COST___VAL CONSTRAINT")
                  DO con = 1 , me%Numcon
                     IF ( me%Contyp(con)==-2 ) CYCLE
                     nam = me%Constr(con)
                     len = me%Conlen(con)
                     val = me%Conval(con)
                     IF ( me%Conact(con)>0 ) THEN
                        WRITE (str,'( I4,5X,6X,D10.3,1X,A)') con , val , nam(1:len)
                        CALL me%ogwrit(3,str)
                     ELSEIF ( me%Conact(con)==0 ) THEN
                        WRITE (str,'( 5X,I4,6X,D10.3,1X,A)') con , val , nam(1:len)
                        CALL me%ogwrit(3,str)
                     ELSEIF ( me%Conact(con)<0 ) THEN
                        WRITE (str,'(10X,I4,1X,D10.3,1X,A)') con , val , nam(1:len)
                        CALL me%ogwrit(3,str)
                     ENDIF
                  ENDDO
                  ! ======================================================================
                  ! ======================================================================
                  IF ( me%Senopt<=0 ) CALL me%ogeval(me%Varval,me%Conval,0,me%Conder)
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ELSE
                  WRITE (str,'("CURPRI=",I3)') curpri
                  CALL me%ogwrit(3,str)
                  curpri = curpri + 1
                  spag_nextblock_1 = 4
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               ! ----------------------------------------------------------------------
            ELSEIF ( cstval<-cornor*varvio ) THEN
               ! ----------------------------------------------------------------------
               WRITE (str,'("CNV=",3(1X,D10.3))') cstval/cornor/varvio
               CALL me%ogwrit(2,str)
               CALL me%ogwrit(3,"")
               WRITE (str,'("CSTVAL=",3D10.3)') cstval , cornor , varvio
               CALL me%ogwrit(3,str)
               IF ( inelop==1 ) THEN
                  WRITE (str,'("INELOP=",I3)') inelop
                  CALL me%ogwrit(3,str)
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
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               IF ( me%Conact(con)/=-1 ) CYCLE
               val = dot_product(me%Conred(con,me%Numact+1:me%Numvar),corvec(me%Numact+1:me%Numvar))
               IF ( val==0.0_wp ) CYCLE
               val = -me%Conval(con)/val
               IF ( val<=0.0_wp ) CYCLE
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
            me%Varval = me%Varval + varvec
            ! ----------------------------------------------------------------------
            DO con = 1 , me%Numcon + 1
               val = dot_product(corvec(me%Numact+1:me%Numvar),me%Conred(con,me%Numact+1:me%Numvar))
               me%Conval(con) = me%Conval(con) + val*foldis
            ENDDO
            ! ----------------------------------------------------------------------
            cstval = cstval + foldis*cornor
            ! ======================================================================
            ! MAXIMUM TRAVEL DISTANCE REACHED: NEXT ITERATION
            ! ----------------------------------------------------------------------
            IF ( ind==0 ) THEN
               WRITE (str,'("CNV=",3(1X,D10.3))') cstval/cornor/varvio
               CALL me%ogwrit(3,str)
               WRITE (str,'(4X,5X,5X,3(1X,D10.3))') cstval , cornor , foldis
               CALL me%ogwrit(3,str)
               IF ( inelop==1 ) THEN
                  WRITE (str,'("INELOP=",I3)') inelop
                  CALL me%ogwrit(3,str)
                  inelop = inelop + 1
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ELSE
                  me%Numact = 0
                  spag_nextblock_1 = 6
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDIF
            ! ======================================================================
            ! CONSTRAINT HIT: UPDATE CONSTRAINTS + CORRECT
            ! ----------------------------------------------------------------------
            con = ind
            nam = me%Constr(con)
            len = me%Conlen(con)
            val = me%Conval(con)
            WRITE (str,'(5X,5X,I4,3(1X,D10.3),1X,A)') con , cstval , cornor , foldis , nam(1:len)
            CALL me%ogwrit(3,str)
            IF ( conhit(con)>=20 ) THEN
               WRITE (str,'("OGCORR-WARNING: CONSTRAINT HIT")')
               CALL me%ogwrit(1,str)
               WRITE (str,'("CON/HIT=",2I4)') con , conhit(con)
               CALL me%ogwrit(1,str)
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
            CALL me%ogwrit(3,"")
            WRITE (str,'("CSTVAL:",D13.6)') cstval
            CALL me%ogwrit(3,str)
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3,"STATUS OF CONSTRAINTS:")
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3," ACT  PAS  NON COST___VAL CONSTRAINT")
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               nam = me%Constr(con)
               len = me%Conlen(con)
               val = me%Conval(con)
               IF ( me%Conact(con)>0 ) THEN
                  WRITE (str,'( I4,5X,6X,D10.3,1X,A)') con , val , nam(1:len)
                  CALL me%ogwrit(3,str)
               ELSEIF ( me%Conact(con)==0 ) THEN
                  WRITE (str,'( 5X,I4,6X,D10.3,1X,A)') con , val , nam(1:len)
                  CALL me%ogwrit(3,str)
               ELSEIF ( me%Conact(con)<0 ) THEN
                  WRITE (str,'(10X,I4,1X,D10.3,1X,A)') con , val , nam(1:len)
                  CALL me%ogwrit(3,str)
               ENDIF
            ENDDO
            ! ----------------------------------------------------------------------
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3,"STATUS OF VIOLATED CONSTRAINTS:")
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3," CON COST___VAL CONSTRAINT")
            conerr = 0.0_wp
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               nam = me%Constr(con)
               len = me%Conlen(con)
               val = me%Conval(con)
               IF ( val<-dlt ) THEN
                  conerr = conerr + abs(val)
                  WRITE (str,'( I4,D11.3,1X,A)') con , val , nam(1:len)
                  CALL me%ogwrit(3,str)
               ELSEIF ( me%Contyp(con)/=0 ) THEN
               ELSEIF ( val>dlt ) THEN
                  conerr = conerr + abs(val)
                  WRITE (str,'( I4,D11.3,1X,A)') con , val , nam(1:len)
                  CALL me%ogwrit(3,str)
               ENDIF
            ENDDO
            ! ----------------------------------------------------------------------
            CALL me%ogwrit(3,"")
            WRITE (str,'("LINEAR ERROR.: ",D13.6)') conerr
            CALL me%ogwrit(3,str)
            spag_nextblock_1 = 7
          CASE (7)
            ! ----------------------------------------------------------------------
            ! ----------------------------------------------------------------------
            me%Contyp = conttt
            me%Conpri = prisav
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

   SUBROUTINE ogeval(me,Valvar,Valcon,Varder,Dercon)

      !! COMPUTES SCALED CONTRAINTS+MERIT AND DERIVATIVES
      !! FROM     SCALED VARIABLES
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Valvar(me%Numvar)
      real(wp),intent(out) :: Valcon(me%Numcon+1)
      integer(ip),intent(in) :: Varder !! DERIVATIVES COMPUTATION MODE
                                       !!
                                       !!  * 0: VALUES ONLY
                                       !!  * 1: USER DEFINED
                                       !!  * 2: NUMERIC WITH DOUBLE DIFFERENCING
                                       !!  * 3: NUMERIC WITH SINGLE DIFFERENCING
      real(wp),intent(out) :: Dercon(me%Numcon+1,me%Numvar)

      integer(ip) :: var , con , cod , len , ind , numvio
      real(wp) :: val , sca , fac , per , sav , der , err , conerr , convio
      character(len=3) :: typ
      character(len=3) :: sta
      character(len=maxstr) :: nam
      character(len=256) :: str

      !real(wp) :: ggg(4,4) , bbb(4) , vvv(4)  ! JW : not sure what this was for
      real(wp) :: objval
      real(wp) , DIMENSION(:) , ALLOCATABLE :: varvec
      real(wp) , DIMENSION(:) , ALLOCATABLE :: convec

      ALLOCATE (varvec(me%Numvar))
      ALLOCATE (convec(me%Numcon+1))

      ! ======================================================================
      ! GENERAL
      ! ----------------------------------------------------------------------
      WRITE (str,'()')
      CALL me%ogwrit(3,str)
      select case (Varder)
       case ( 0 );    WRITE (str,'("COMPUTE RESULTS")')
       case ( 1, -1); WRITE (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES USER DEFINED")')
       case ( 2 );    WRITE (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES BY DOUBLE DIFFERENCING")')
       case ( 3 );    WRITE (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES BY SINGLE DIFFERENCING")')
      end select
      CALL me%ogwrit(3,str)

      ! ======================================================================
      ! WRITE VARIABLES
      ! ----------------------------------------------------------------------
      WRITE (str,'()')
      CALL me%ogwrit(3,str)
      WRITE (str,'("VARIABLES NOT SCALED:")')
      CALL me%ogwrit(3,str)
      WRITE (str,'()')
      CALL me%ogwrit(3,str)
      DO var = 1 , me%Numvar
         val = Valvar(var)
         sca = me%Varsca(var)
         cod = me%Vartyp(var)
         select case (cod)
          case ( 0 ); typ = "FRE"
          case ( 1 ); typ = "PAR"
         end select
         nam = me%Varstr(var)
         len = me%Varlen(var)
         WRITE (str,'("VAR/VAL/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') &
            var , val*sca , sca , typ , nam(1:len)
         CALL me%ogwrit(3,str)
      ENDDO

      ! ======================================================================
      ! DE-SCALE VARIABLES
      ! ----------------------------------------------------------------------
      DO var = 1 , me%Numvar
         varvec(var) = Valvar(var)*me%Varsca(var)
      ENDDO

      ! ======================================================================
      ! GET RESULTS
      ! GET DERIVATIVES IF USER DEFINED
      ! ----------------------------------------------------------------------
      select case (Varder)
       case ( 0 );     CALL me%calval(varvec,Valcon,0)
       case ( 1, -1 ); CALL me%calval(varvec,Valcon,1)
       case ( 2 );     CALL me%calval(varvec,Valcon,1)
       case ( 3 );     CALL me%calval(varvec,Valcon,1)
      end select

! ======================================================================
!    IF ( 1==2 ) THEN
!       ggg(1,1) = +1D+01
!       ggg(2,1) = +1D+00
!       ggg(3,1) = +2D+00
!       ggg(4,1) = +3D+00
!       ggg(1,2) = ggg(2,1)
!       ggg(2,2) = +1D+01
!       ggg(3,2) = +4D+00
!       ggg(4,2) = +5D+00
!       ggg(1,3) = ggg(3,1)
!       ggg(2,3) = ggg(3,2)
!       ggg(3,3) = +1D+01
!       ggg(4,3) = +6D+00
!       ggg(1,4) = ggg(4,1)
!       ggg(2,4) = ggg(4,2)
!       ggg(3,4) = ggg(4,3)
!       ggg(4,4) = +1D+01
!       ! ----------------------------------------------------------------------
!       bbb(1) = +1D+01
!       bbb(2) = +1D+01
!       bbb(3) = +1D+01
!       bbb(4) = +1D+01
!       ! ----------------------------------------------------------------------
!       CALL mul2m(ggg,4,1,1,4,Valvar,4,1,1,4,vvv,4,1,1,1)
!       CALL mul2m(Valvar,1,1,1,1,vvv,4,1,1,4,Valcon,1,1,1,1)
!       CALL mul2m(bbb,1,1,1,1,Valvar,4,1,1,4,vvv,1,1,1,1)
!       CALL mulvs(Valcon,0.5_wp,Valcon,1)
!       CALL sum2v(Valcon,vvv,Valcon,1)
!       CALL mulvs(Valcon,-1.0_wp,Valcon,1)
!       WRITE (str,*) "VALCON=" , (Valcon(ind),ind=1,1)
!       CALL me%ogwrit(3,str)
!    ENDIF

      ! ----------------------------------------------------------------------
      ! SCALE RESULTS
      ! ----------------------------------------------------------------------
      DO con = 1 , me%Numcon + 1
         convec(con) = Valcon(con)
         sca = me%Consca(con)
         cod = me%Contyp(con)
         IF ( cod==-1 ) sca = -sca
         Valcon(con) = Valcon(con)/sca
      ENDDO

      ! ======================================================================
      ! WRITE RESULTS
      ! ----------------------------------------------------------------------
      WRITE (str,'()')
      CALL me%ogwrit(3,str)
      WRITE (str,'("RESULTS NOT SCALED:")')
      CALL me%ogwrit(3,str)
      WRITE (str,'()')
      CALL me%ogwrit(3,str)
      conerr = 0.0_wp  ! total constraint error (scaled to constr. threshod)
      convio = 0.0_wp  ! total constaint error norm (unscaled)
      ind = 0          ! index of largest constraint violation
      fac = 0.0_wp     ! value of largest constraint violation
      numvio = 0       ! number of violated constraints
      DO con = 1 , me%Numcon + 1
         val = Valcon(con)
         sca = me%Consca(con)
         cod = me%Contyp(con)
         sta = "   "
         err = 0.0_wp
         IF ( cod==-1 ) sca = -sca
         IF ( cod==-2 ) typ = "DER"
         IF ( cod==-1 .AND. con<=me%Numcon ) typ = "LTE"
         IF ( cod==0 .AND. con<=me%Numcon ) typ = "EQU"
         IF ( cod==1 .AND. con<=me%Numcon ) typ = "GTE"
         IF ( cod==-1 .AND. con>me%Numcon ) typ = "MIN"
         IF ( cod==1 .AND. con>me%Numcon ) typ = "MAX"
         IF ( cod==0 .AND. con<=me%Numcon .AND. abs(val)>1.0_wp ) THEN
            sta = "VIO"
            err = abs(val)
            numvio = numvio + 1
         ENDIF
         IF ( cod/=0 .AND. con<=me%Numcon .AND. cod/=-2 .AND. -val>1.0_wp ) THEN
            sta = "VIO"
            err = abs(val)
            numvio = numvio + 1
         ENDIF
         conerr = conerr + err
         convio = convio + (err*sca)**2
         IF ( err>fac ) ind = con
         IF ( err>fac ) fac = err
         nam = me%Constr(con)
         len = me%Conlen(con)
         WRITE (str,'("CON/VAL/SCA/TYP/STA/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A3,1X,A)') &
            con , val*sca , sca , typ , sta , nam(1:len)
         CALL me%ogwrit(3,str)
      ENDDO
      WRITE (str,'()')
      CALL me%ogwrit(3,str)
      WRITE (str,'("CONSTRAINT ERROR.:",2(1X,D13.6),I6)') conerr , fac , ind
      CALL me%ogwrit(3,str)
      WRITE (str,'()')
      CALL me%ogwrit(3,str)
      ! write pygmo-style log output
      objval = -Valcon(me%Numcon+1)
      convio = sqrt(convio)
      CALL me%ogpwri(objval,numvio,convio)

      ! ======================================================================
      ! NO DERIVATIVES
      ! ----------------------------------------------------------------------
      IF ( Varder==0 ) THEN
         RETURN
      ELSEIF ( Varder==1 .OR. Varder==-1 ) THEN
         CALL me%calder(varvec,convec,Dercon)
      ENDIF

! ----------------------------------------------------------------------
!    IF ( 1==2 ) THEN
!       CALL mul2m(Valvar,1,1,1,1,ggg,4,1,1,4,Dercon,1,1,1,4)
!       CALL sum2v(Dercon,bbb,Dercon,4)
!       CALL mulvs(Dercon,-1.0_wp,Dercon,4)
!       WRITE (str,*) "DERCON=" , (Dercon(1,ind),ind=1,4)
!       CALL me%ogwrit(3,str)
!    ENDIF

      ! ======================================================================
      ! WRITE DERIVATIVES
      ! ----------------------------------------------------------------------
      WRITE (str,'()')
      CALL me%ogwrit(3,str)
      WRITE (str,'("DERIVATIVES SCALED:")')
      CALL me%ogwrit(3,str)
      WRITE (str,'()')
      CALL me%ogwrit(3,str)

      DO var = 1 , me%Numvar

         ! ----------------------------------------------------------------------
         ! WRITE VARIABLE
         ! ----------------------------------------------------------------------
         val = Valvar(var)
         sca = me%Varsca(var)
         cod = me%Vartyp(var)
         IF ( cod==0 ) typ = "FRE"
         IF ( cod==1 ) typ = "PAR"
         nam = me%Varstr(var)
         len = me%Varlen(var)
         WRITE (str,'("VAR/VAL/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') var , val*sca , sca , typ , nam(1:len)
         CALL me%ogwrit(4,str)
         WRITE (str,'()')
         CALL me%ogwrit(4,str)

         select case (Varder)
          case ( 2 )  ! DERIVATIVES BY DOUBLE DIFFERENCING
            per = me%Varper(var)
            sav = varvec(var)
            varvec(var) = sav + per
            CALL me%calval(varvec,Dercon(1:,var),0)  ! JW added : here
            varvec(var) = sav - per
            CALL me%calval(varvec,convec,0)
            fac = 0.5_wp/per
            DO con = 1 , me%Numcon + 1
               Dercon(con,var) = (Dercon(con,var)-convec(con))*fac
            ENDDO
            varvec(var) = sav
          case ( 3 )  ! DERIVATIVES BY SINGLE DIFFERENCING
            per = me%Varper(var)
            sav = varvec(var)
            varvec(var) = sav + per
            CALL me%calval(varvec,Dercon(1:,var),0)  ! JW added : here
            fac = 1.0_wp/per
            DO con = 1 , me%Numcon + 1
               Dercon(con,var) = (Dercon(con,var)-convec(con))*fac
            ENDDO
            varvec(var) = sav
         end select

         ! ----------------------------------------------------------------------
         ! SCALE DERIVATIVES
         ! ----------------------------------------------------------------------
         DO con = 1 , me%Numcon + 1
            fac = me%Varsca(var)/me%Consca(con)
            IF ( me%Contyp(con)==-1 ) fac = -fac
            Dercon(con,var) = Dercon(con,var)*fac
         ENDDO

         ! ----------------------------------------------------------------------
         ! WRITE DERIVATIVES
         ! ======================================================================
         DO con = 1 , me%Numcon + 1
            der = Dercon(con,var)
            IF ( der/=0.0_wp ) THEN
               sca = me%Consca(con)
               cod = me%Contyp(con)
               IF ( cod==-1 ) sca = -sca
               IF ( cod==-2 ) typ = "DER"
               IF ( cod==-1 .AND. con<=me%Numcon ) typ = "LTE"
               IF ( cod==0 .AND. con<=me%Numcon ) typ = "EQU"
               IF ( cod==1 .AND. con<=me%Numcon ) typ = "GTE"
               IF ( cod==-1 .AND. con>me%Numcon ) typ = "MIN"
               IF ( cod==1 .AND. con>me%Numcon ) typ = "MAX"
               nam = me%Constr(con)
               len = me%Conlen(con)
               WRITE (str,'("CON/DER/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') &
                  con , der*sca/me%Varsca(var) , sca , typ , nam(1:len)
               CALL me%ogwrit(4,str)
            ENDIF
         ENDDO
         WRITE (str,'()')
         CALL me%ogwrit(4,str)

      ENDDO

      DEALLOCATE (varvec)
      DEALLOCATE (convec)

   END SUBROUTINE ogeval

   SUBROUTINE ogexcl(me,Exc)

      !! REMOVE CONSTRAINT TO ACTIVE SET AND REDUCES DERIVATIVES
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      integer(ip),intent(in) :: Exc !! CONSTRAINT TO BE REMOVED
                                    !! SEQUENCE NUMBER IN ACTIVE LIST

      real(wp) :: val , bet , gam
      integer(ip) :: row , col , act , con
      CHARACTER(len=256) :: str

      ! ======================================================================
      ! ADJUST LIST OF ACTIVE CONSTRAINTS
      ! ----------------------------------------------------------------------
      con = me%Actcon(Exc)
      me%Conact(con) = 0
      me%Numact = me%Numact - 1
      DO act = Exc , me%Numact
         con = me%Actcon(act+1)
         me%Actcon(act) = con
         me%Conact(con) = me%Conact(con) - 1
      ENDDO
      ! ======================================================================
      ! REDUCE FOR SUBSEQUENT CONSTRAINTS
      ! ----------------------------------------------------------------------
      DO act = Exc , me%Numact
         con = me%Actcon(act)
         val = 0.0_wp
         DO col = act , act + 1
            val = val + me%Conred(con,col)**2
         ENDDO
         val = sqrt(val)
         IF ( me%Conred(con,act)>0.0_wp ) val = -val
         IF ( abs(val)<1.0e-15_wp ) THEN
            WRITE (me%Loglun,*) "OGEXCL-ERROR: CONSTRAINTS SINGULAR"
            CALL me%ogwrit(2,str)
            WRITE (me%Loglun,*) "VAL=" , val
            CALL me%ogwrit(2,str)
            STOP
         ENDIF
         me%Conred(con,act) = me%Conred(con,act) - val
         bet = 1.0_wp/(val*me%Conred(con,act))
         DO row = 1 , me%Numcon + 3
            IF ( me%Conact(row)>act .OR. me%Conact(row)<=0 ) THEN
               gam = 0.0_wp
               DO col = act , act + 1
                  IF ( me%Conred(row,col)/=0.0_wp ) gam = gam + me%Conred(row,col)*me%Conred(con,col)
               ENDDO
               IF ( gam/=0.0_wp ) THEN
                  gam = gam*bet
                  DO col = act , act + 1
                     me%Conred(row,col) = me%Conred(row,col) + me%Conred(con,col)*gam
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         me%Conred(con,act) = val
         DO col = act + 1 , act + 1
            me%Conred(con,col) = 0.0_wp
         ENDDO
      ENDDO
      ! ======================================================================
   END SUBROUTINE ogexcl

   SUBROUTINE ogexec(me,Valvar,Valcon,Finopt,Finite)
      !! Main routine.
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(inout) :: Valvar(me%Numvar) !! VARIABLES VALUE
                                                  !! -> NOT SCALED
      real(wp),intent(out) :: Valcon(me%Numcon+1) !! CONSTRAINTS VALUE (1:NUMCON)
                                                  !! MERIT       VALUE (1+NUMCON)
                                                  !! -> NOT SCALED
      integer(ip),intent(out) :: Finopt !! TERMINATION STATUS
                                        !!
                                        !!  * 1=    MATCHED &     OPTIMAL
                                        !!  * 2=    MATCHED & NOT OPTIMAL
                                        !!  * 3=NOT MATCHED & NOT OPTIMAL
                                        !!  * 4=NOT FEASIBL & NOT OPTIMAL
      integer(ip),intent(out) :: Finite

      integer(ip) :: finish , itecor , iteopt
      integer(ip) :: var , con , typ , len , num , numvio
      real(wp) :: val , sca , red , der , fac , old , convio
      CHARACTER(len=256) :: str , nam

      integer(ip) :: numequ , itediv , itecnv
      real(wp) :: varacc , cosnew , cosold , varsav , meamer
      real(wp) :: conerr , desnor , norerr , meaerr

      real(wp) , DIMENSION(:) , ALLOCATABLE :: varsum
      real(wp) , DIMENSION(:) , ALLOCATABLE :: varcor
      real(wp) , DIMENSION(:) , ALLOCATABLE :: concor
      INTEGER :: spag_nextblock_1

      spag_nextblock_1 = 1
      SPAG_DispatchLoop_1: DO
         SELECT CASE (spag_nextblock_1)
          CASE (1)
            ! ----------------------------------------------------------------------
            ALLOCATE (varsum(me%Numvar))
            ALLOCATE (varcor(me%Numvar))
            ALLOCATE (concor(me%Numcon+1))
            ! ======================================================================
            ! GENERAL
            ! ----------------------------------------------------------------------
            ! LOGLEV = 2
            CALL me%ogwrit(2,"")
            CALL me%ogwrit(2,"OPTGRA START")
            CALL me%ogwrit(2,"")
            Finopt = 3
            itecor = 0
            iteopt = 0
            meaerr = 0.0_wp
            meamer = 0.0_wp
            itediv = 0
            itecnv = 0
            !me%Conopt = 0
            concor = 0.0_wp
            varcor = 0.0_wp
            desnor = 0.0_wp
            ! ----------------------------------------------------------------------
            me%Varstp = me%Varsnd
            me%Numite = 0
            cosnew = 0.0_wp
            DO var = 1 , me%Numvar
               varsum(var) = 0.0_wp
               varcor(var) = 0.0_wp
            ENDDO
            ! ======================================================================
            ! EQUALTIY CONSTRAINTS IN ACTIVE SET
            ! ----------------------------------------------------------------------
            me%Numact = 0
            DO con = 1 , me%Numcon
!              IF (CONTYP(CON) == -2) CYCLE
               nam = me%Constr(con)
               len = me%Conlen(con)
               WRITE (str,*) "CON/PRI=" , con , me%Conpri(con) , " " , nam(1:len)
               CALL me%ogwrit(3,str)
               me%Conact(con) = 0
               IF ( me%Consca(con)>=1.0e+09_wp ) me%Contyp(con) = -2
               IF ( me%Contyp(con)==0 ) THEN
               ELSEIF ( me%Contyp(con)==-2 ) THEN
                  me%Conact(con) = -2
               ENDIF
            ENDDO
            numequ = me%Numact
            me%Conact(me%Numcon+1) = -3
            me%Conact(me%Numcon+2) = -3
            ! ======================================================================
            ! SCALE VARIABLES
            ! ----------------------------------------------------------------------
            DO var = 1 , me%Numvar
               me%Varval(var) = Valvar(var)/me%Varsca(var)
            ENDDO
            ! ======================================================================
            ! HEADER FOR TABLE
            ! ----------------------------------------------------------------------
            IF ( me%Tablev>=1 ) WRITE (me%Tablun,'("ITER",1X,"OPT",1X,*(1X,I10))') (var,var=1,me%Numvar) , (con,con=1,me%Numcon)
            spag_nextblock_1 = 2
          CASE (2)
            SPAG_Loop_1_1: DO
               ! ======================================================================
               !      IF (NUMITE >= 52) MATLEV = 3
               !      IF (NUMITE >= 55) MATLEV = 2
               ! ======================================================================
               ! NEW ITERATION
               ! ----------------------------------------------------------------------
               IF ( me%Numite>=me%Corite .AND. itecor==0 ) THEN
                  Finopt = 3
                  Finite = me%Numite
                  CALL me%ogwrit(1,"")
                  WRITE (str,'("OPTGRA: Converged: not ITERAT=",2I4,2D11.3)') me%Numite , me%Maxite , conerr , desnor
                  CALL me%ogwrit(1,str)

                  ! Final Pygmo output
                  ! TODO: can this final fitness call be avoided (just for output)?
                  me%Pygfla = 3 ! pygmo flag in COMMON: no covergence
                  CALL me%ogeval(me%Varval,me%Conval,me%Varder,me%Conder(1:me%Numcon+1,:))
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ELSEIF ( me%Numite>=me%Maxite .OR. (me%Numite-itecor>=me%Optite-1 .AND. itecor/=0) ) THEN
                  Finopt = 2
                  Finite = iteopt
                  me%Varval = varcor
                  me%Conval = concor
                  CALL me%ogwrit(1,"")
                  WRITE (str,'("OPTGRA: Converged: mat ITERAT=",2I4,2D11.3)') me%Numite , me%Maxite , conerr , desnor
                  CALL me%ogwrit(1,str)
                  CALL me%ogpwri_end(-Valcon(me%Numcon+1),numvio,convio)

                  ! Final Pygmo output
                  ! TODO: can this final fitness call be avoided (just for output)?
                  me%Pygfla = 2 ! pygmo flag in COMMON: constraints matched
                  CALL me%ogeval(me%Varval,me%Conval,me%Varder,me%Conder(1:me%Numcon+1,:))
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               ! ----------------------------------------------------------------------
               me%Numite = me%Numite + 1
               ! ----------------------------------------------------------------------
               CALL me%ogwrit(3,"")
               WRITE (str,'("ITERAT=",I5)')me%Numite
               CALL me%ogwrit(2,str)
               ! ======================================================================
               ! GET VALUES AND GRADIENTS
               ! ======================================================================
               IF ( me%Senopt<=0 ) THEN
                  CALL me%ogeval(me%Varval,me%Conval,me%Varder,me%Conder(1:me%Numcon+1,:))
               ELSEIF ( me%Senopt==+1 .OR. me%Senopt==+3 ) THEN
                  me%Varval = me%Senvar
                  CALL me%ogeval(me%Varval,me%Conval,0,me%Conder(1:me%Numcon+1,:))
               ELSEIF ( me%Senopt==+2 .OR. me%Senopt==+4 ) THEN
                  me%Varval = me%Senvar
                  DO con = 1 , me%Numcon + 1
                     sca = me%Consca(con)
                     IF ( me%Contyp(con)==-1 ) sca = -sca
                     me%Conval(con) = me%Sencon(con) - me%Sendel(con)/sca
                  ENDDO
               ENDIF
               IF ( me%Senopt==-1 ) THEN
                  me%Senvar = me%Varval
                  me%Sencon = me%Conval
               ENDIF
               ! ======================================================================
               IF ( me%Varder==-1 .AND. me%Senopt<=0 ) THEN
                  me%Conred(1:me%Numcon+1,:) = me%Conder(1:me%Numcon+1,:)
                  CALL me%ogeval(me%Varval,me%Conval,2,me%Conder(1:me%Numcon+1,:))
                  WRITE (str,'("GRADIENT CHECK")')
                  CALL me%ogwrit(1,str)
                  DO var = 1 , me%Numvar
                     DO con = 1 , me%Numcon + 1
                        fac = me%Varsca(var)/me%Consca(con)
                        fac = 1.0_wp
                        der = me%Conder(con,var)*fac
                        red = me%Conred(con,var)*fac
                        IF ( abs(der)<1.0e-6_wp .AND. abs(red)<1.0e-6_wp ) CYCLE
                        IF ( abs(der-red)<1.0e-2_wp ) CYCLE
                        IF ( der/=0.0_wp ) THEN
                           fac = red/der
                        ELSE
                           fac = 0.0_wp
                        ENDIF
                        IF ( abs(fac-1.0_wp)<1.0e-2_wp ) CYCLE
                        WRITE (str,'("VAR/CON/ANA/NUM/A2N=",2I4,3(1X,D13.6))') var , con , red , der , fac
                        CALL me%ogwrit(1,str)
                        nam = me%Varstr(var)
                        len = me%Varlen(var)
                        WRITE (str,'("      VAR=",A,1X,D13.6)') nam(1:len) , me%Varval(var)*me%Varsca(var)
                        CALL me%ogwrit(1,str)
                        nam = me%Constr(con)
                        len = me%Conlen(con)
                        WRITE (str,'("      CON=",A,1X,D13.6)') nam(1:len) , me%Conval(con)*me%Consca(con)
                        CALL me%ogwrit(1,str)
                     ENDDO
                  ENDDO
!                 CONDER(1:NUMCON+1,:) = CONRED(1:NUMCON+1,:)
!                 GOTO 9999
               ENDIF
               ! ======================================================================
               me%Sender = me%Conder
               DO var = 1 , me%Numvar
                  IF ( me%Vartyp(var)/=1 ) CYCLE
!                 WRITE (STR,*) "VAR=",VAR,VARVAL(VAR)*VARSCA(VAR)
!                 CALL me%ogwrit (2,STR)
                  me%Conder(1:me%Numcon+1,var) = 0.0_wp
               ENDDO
               ! ======================================================================
               IF ( me%Numite==1 ) THEN
                  me%Vargrd = me%Varval
               ELSE
                  me%Vargrd = me%Varref
               ENDIF
               me%Varref = me%Varval
               me%Conref = me%Conval
               ! ======================================================================
               varacc = 0.0_wp
               ! ======================================================================
               cosold = cosnew
               cosnew = me%Conval(me%Numcon+1)
               CALL me%ogwrit(3,"")
               WRITE (str,'("OPTGRA: VALCOS=",D15.8,1X,D15.8)') cosnew , cosnew - cosold
               CALL me%ogwrit(3,str)
               ! ======================================================================
               ! CORRECTION PART
               ! ----------------------------------------------------------------------
               CALL me%ogcorr(varacc,finish,conerr,norerr)
               ! ----------------------------------------------------------------------
               IF ( me%Tablev>=1 ) WRITE (me%Tablun,'(I4,1X,"COR",1X,*(1X,D10.3))') &
                  me%Numite , (me%Varval(var),var=1,me%Numvar) , (me%Conval(con),con=1,me%Numcon)
               ! ----------------------------------------------------------------------
               IF ( me%Senopt/=0 ) THEN
                  IF ( finish/=0 ) EXIT SPAG_Loop_1_1
                  Finopt = 0
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               ! ----------------------------------------------------------------------
               IF ( finish==0 ) THEN
                  me%Numact = 0
                  old = meaerr
                  itediv = itediv + 1
                  num = min(itediv,me%Divite)
                  meaerr = (meaerr*(num-1)+norerr)/num
!                 WRITE (STR,*) "MEAERR=",MEAERR
!                 CALL me%ogwrit (2,STR)
                  IF ( itediv<me%Divite .OR. meaerr<=old ) CYCLE
                  finish = -1
               ENDIF
               ! ----------------------------------------------------------------------
               IF ( finish==-1 ) THEN
                  Finopt = 4
                  Finite = me%Numite
                  CALL me%ogwrit(1,"")
                  WRITE (str,'("OPTGRA: Converged: unf ITERAT=",2I4,2D11.3)') me%Numite , me%Maxite , conerr , desnor
                  CALL me%ogwrit(1,str)

                  ! Final Pygmo output
                  ! TODO: can this final fitness call be avoided (just for output)?
                  me%Pygfla = 4 ! pygmo flag in COMMON: infeasible
                  CALL me%ogeval(me%Varval,me%Conval,me%Varder,me%Conder(1:me%Numcon+1,:))
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               ! ----------------------------------------------------------------------
               itediv = 0
               iteopt = me%Numite
               IF ( itecor==0 .OR. concor(me%Numcon+1)<me%Conval(me%Numcon+1) ) THEN
                  varcor = me%Varval
                  concor = me%Conval
               ENDIF
               IF ( itecor==0 ) itecor = me%Numite
               ! ----------------------------------------------------------------------
               old = meamer
               itecnv = itecnv + 1
               num = min(itecnv,me%Cnvite)
               meamer = (meamer*(num-1)+concor(me%Numcon+1))/num
!              WRITE (STR,*) "MEAMER=",ITECNV,NUM,MEAMER,OLD,OLD/MEAMER
!              CALL me%ogwrit (-1,STR)
               IF ( itecnv>=me%Cnvite .AND. meamer<old ) THEN
                  Finopt = 2
                  Finite = iteopt
                  me%Varval = varcor
                  me%Conval = concor
                  CALL me%ogwrit(1,"")
                  WRITE (str,'("OPTGRA: Converged: mat ITERAT=",2I4,2D11.3)') me%Numite , me%Maxite , conerr , desnor
                  CALL me%ogwrit(1,str)

                  ! Final Pygmo output
                  ! TODO: can this final fitness call be avoided (just for output)?
                  me%Pygfla = 2 ! pygmo flag in COMMON: matched
                  CALL me%ogeval(me%Varval,me%Conval,me%Varder,me%Conder(1:me%Numcon+1,:))
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               EXIT SPAG_Loop_1_1
            ENDDO SPAG_Loop_1_1
            ! ======================================================================
            ! OPTIMISATION PART
            ! ----------------------------------------------------------------------
            ! ----------------------------------------------------------------------
            IF ( me%Senopt<+3 ) THEN
               varsav = me%Varmax
               me%Varmax = me%Varmax*10.0e-1_wp
               CALL me%ogopti(varacc,numequ,finish,desnor)
               me%Varmax = varsav
            ENDIF
            ! ----------------------------------------------------------------------
            IF ( me%Senopt/=0 ) THEN
               CALL me%ogwrit(1,"")
               IF ( finish==0 ) THEN
                  Finopt = 0
                  CALL me%ogwrit(1,"OPTGRA sensitivity converged: not")
               ELSE
                  Finopt = 1
                  CALL me%ogwrit(1,"OPTGRA sensitivity converged: yes")
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
            IF ( varacc/=0.0_wp ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            ! ======================================================================
            ! CONVERGED
            ! ----------------------------------------------------------------------
            Finopt = 1
            Finite = me%Numite
            CALL me%ogwrit(1,"")
            WRITE (str,'("OPTGRA: Converged: yes ITERAT=",2I4,2D11.3)') me%Numite , me%Maxite , conerr , desnor
            CALL me%ogwrit(1,str)
            CALL me%ogwrit(3,"")

            ! Final Pygmo output
            ! TODO: can this final fitness call be avoided (just for output)?
            me%Pygfla = 1 ! covergence
            CALL me%ogeval(me%Varval,me%Conval,me%Varder,me%Conder(1:me%Numcon+1,:))
            spag_nextblock_1 = 3
          CASE (3)

!           WRITE (STR,*) "DIF=",NORM2(VARVAL-VARREF)
!           CALL me%ogwrit (1,STR)
            ! ======================================================================
            ! ======================================================================
            ! DESCALE VARIABLES
            ! ----------------------------------------------------------------------
!           CALL OGWMAT (3)
!           CALL me%ogeval (VARVAL, VALCON, 0, CONDER)
            DO var = 1 , me%Numvar
               Valvar(var) = me%Varval(var)*me%Varsca(var)
            ENDDO
!           IF (SENOPT /= 0) THEN
!               CALL me%ogeval (VARVAL, VALCON, 0, CONDER)
!           ENDIF
            ! ======================================================================
            ! DESCALE VALUES
            ! ----------------------------------------------------------------------
            DO con = 1 , me%Numcon + 1
               typ = me%Contyp(con)
               sca = me%Consca(con)
               IF ( typ==-1 ) sca = -sca
               Valcon(con) = me%Conval(con)*sca
            ENDDO
            ! ======================================================================
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3,"STATUS OF CONSTRAINTS:")
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3," ACT  PAS  NON COST___VAL CONSTRAINT")
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               nam = me%Constr(con)
               len = me%Conlen(con)
               val = me%Conval(con)
               IF ( me%Conact(con)>0 ) THEN
                  WRITE (str,'( I4,5X,6X,D10.3,1X,A)') con , val , nam(1:len)
                  CALL me%ogwrit(3,str)
               ELSEIF ( me%Conact(con)==0 ) THEN
                  WRITE (str,'( 5X,I4,6X,D10.3,1X,A)') con , val , nam(1:len)
                  CALL me%ogwrit(3,str)
               ELSEIF ( me%Conact(con)<0 ) THEN
                  WRITE (str,'(10X,I4,1X,D10.3,1X,A)') con , val , nam(1:len)
                  CALL me%ogwrit(3,str)
               ENDIF
            ENDDO

            DEALLOCATE (varsum)
            DEALLOCATE (varcor)
            DEALLOCATE (concor)
            ! ----------------------------------------------------------------------
            CALL me%ogwrit(2,"")
            CALL me%ogwrit(2,"OPTGRA END")
            CALL me%ogwrit(2,"")
            EXIT SPAG_DispatchLoop_1
         END SELECT
      ENDDO SPAG_DispatchLoop_1
      ! ======================================================================
   END SUBROUTINE ogexec

   SUBROUTINE oggsst(me,Varsen,Quasen,Consen,Actsen,Dersen,Actsav,Consav,Redsav,Dersav,Actnum)

      !! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
      !!
      !! Function to get sensitivity state data, necessary for serialization.
      !! Do not use this directly except in serialization routines
      !!
      !! 2021/07/19 | M. von Looz | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(out) :: Varsen(me%Numvar) !! STORED VARIABLES VALUE
      real(wp),intent(out) :: Quasen(me%Numcon+1) !! STORED CONSTRAINTS CORRECTION VECTOR
      real(wp),intent(out) :: Consen(me%Numcon+1) !! STORED CONSTRAINTS VALUE
      integer(ip),intent(out) :: Actsen(me%Numcon+1) !! STORED CONSTRAINTS ACTIVE
      real(wp),intent(out) :: Dersen(me%Numcon+1,me%Numvar) !! STORED DERIVATIVE
      integer(ip),intent(out) :: Actsav(me%Numcon+1) !! STORED ACTIVE CONSTRAINTS
      integer(ip),intent(out) :: Consav(me%Numcon+4) !! STORED ACTIVE CONSTRAINTS
      real(wp),intent(out) :: Redsav(me%Numcon+3,me%Numvar) !! STORED DERIVATIVE
      real(wp),intent(out) :: Dersav(me%Numcon+3,me%Numvar) !! STORED DERIVATIVE
      integer(ip),intent(out) :: Actnum !! NUMBER OF ACTIVE CONSTRAINTS

      integer(ip) :: var , con

      ! Variable values saved for sensitivity
      Actnum = me%Numact

      DO var = 1 , me%Numvar
         Varsen(var) = me%Senvar(var)
      ENDDO

      DO con = 1 , me%Numcon + 1
         Quasen(con) = me%Senqua(con)
         Consen(con) = me%Sencon(con)
         Actsen(con) = me%Senact(con)
         DO var = 1 , me%Numvar
            Dersen(con,var) = me%Sender(con,var)
         ENDDO
      ENDDO

      ! Temporary status saved of which constraints are active
      DO con = 1 , me%Numcon + 1
         Actsav(con) = me%Actcon(con)
      ENDDO

      DO con = 1 , me%Numcon + 4
         Consav(con) = me%Conact(con)
      ENDDO

      DO con = 1 , me%Numcon + 3
         DO var = 1 , me%Numvar
            Redsav(con,var) = me%Conred(con,var)
            Dersav(con,var) = me%Conder(con,var)
         ENDDO
      ENDDO

   END SUBROUTINE oggsst

   SUBROUTINE ogincl(me,Inc)

      !! ADDS CONSTRAINT TO ACTIVE SET AND REDUCES DERIVATIVES
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      integer(ip),intent(in) :: Inc !! CONSTRAINT TO BE INCLUDED

      real(wp) :: val , fac , gam , sav , max
      integer(ip) :: row , col , ind , lst
      CHARACTER(len=256) :: str

      ! GENERAL
      me%Numact = me%Numact + 1

      ! PERMUTATION TO GET ZERO DERIVATIVES AT END FOR NEW ACTIVE CONSTRAINT
      lst = me%Numvar
      DO col = me%Numvar , me%Numact , -1
         IF ( me%Conred(Inc,col)==0.0_wp ) THEN
            IF ( col/=lst ) THEN
               DO row = 1 , me%Numcon + 3
                  IF ( me%Conact(row)<=0 ) THEN
                     sav = me%Conred(row,col)
                     me%Conred(row,col) = me%Conred(row,lst)
                     me%Conred(row,lst) = sav
                  ENDIF
               ENDDO
            ENDIF
            lst = lst - 1
         ENDIF
      ENDDO

      ! PERMUTATION TO GET MAXIMUM PIVOT
      ind = me%Numact
      max = abs(me%Conred(Inc,ind))
      DO col = me%Numact + 1 , lst
         val = abs(me%Conred(Inc,col))
         IF ( val>max ) THEN
            ind = col
            max = val
         ENDIF
      ENDDO

      IF ( ind/=me%Numact ) THEN
         DO row = 1 , me%Numcon + 3
            IF ( me%Conact(row)<=0 ) THEN
               sav = me%Conred(row,ind)
               me%Conred(row,ind) = me%Conred(row,me%Numact)
               me%Conred(row,me%Numact) = sav
            ENDIF
         ENDDO
      ENDIF

      ! UPDATE LIST OF ACTIVE CONSTRAINTS
      me%Actcon(me%Numact) = Inc
      me%Conact(Inc) = me%Numact

      ! REDUCE FOR NEW ACTIVE CONSTRAINT
      IF ( abs(me%Conred(Inc,me%Numact))<1.0e-12_wp ) THEN
         WRITE (str,*) "OGINCL-WARNING: CONSTRAINT SINGULAR"
         CALL me%ogwrit(2,str)
         WRITE (str,*) "INC=" , Inc
         CALL me%ogwrit(2,str)
         WRITE (str,*) "PIV=" , me%Conred(Inc,me%Numact)
         CALL me%ogwrit(2,str)
         me%Numact = me%Numact - 1
         me%Conact(Inc) = 0
         RETURN
      ENDIF

      val = sqrt(sum(me%Conred(Inc,me%Numact:lst)**2))
      IF ( me%Conred(Inc,me%Numact)>0.0_wp ) val = -val

      me%Conred(Inc,me%Numact) = me%Conred(Inc,me%Numact) - val

      sav = me%Conred(Inc,me%Numact)
      fac = 1.0_wp/sav
      me%Conred(Inc,me%Numact:lst) = me%Conred(Inc,me%Numact:lst)*fac

      fac = sav/val
      DO row = 1 , me%Numcon + 3
         IF ( me%Conact(row)<=0 ) THEN
            gam = dot_product(me%Conred(row,me%Numact:lst),me%Conred(Inc,me%Numact:lst))
            IF ( gam/=0.0_wp ) THEN
               gam = gam*fac
               me%Conred(row,me%Numact:lst) = me%Conred(row,me%Numact:lst) + me%Conred(Inc,me%Numact:lst)*gam
            ENDIF
         ENDIF
      ENDDO

      me%Conred(Inc,me%Numact) = val
      me%Conred(Inc,me%Numact+1:lst) = 0.0_wp

   END SUBROUTINE ogincl

   SUBROUTINE ogleft(me,Actinp,Actout)

      !! LEFT-MULTIPLIES VECTOR LOWER TRIANGULAR MATRIX OBTAINED BY REDUCTION
      !! AND SUBSEQUENT INVERSION OF DERIVATIVES OF ACTIVE CONSTRAINTS
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Actinp(me%Numcon) !! VECTOR INITAL
      real(wp),intent(out) :: Actout(me%Numcon) !! VECTOR FINAL (MAY BE SAME AS ACTINP)

      integer(ip) :: row , col , act
      real(wp) :: val

      DO act = 1 , me%Numact
         row = me%Actcon(act)
         val = Actinp(act)
         DO col = 1 , act - 1
            val = val - me%Conred(row,col)*Actout(col)
         ENDDO
         Actout(act) = val/me%Conred(row,act)
      ENDDO

   END SUBROUTINE ogleft

   SUBROUTINE ogopti(me,Varacc,Numequ,Finish,Desnor)

      !! OPTIMISATION PART
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(inout) :: Varacc !! ITERATION SCALED DISTANCE ACCUMULATED
      integer(ip),intent(in) :: Numequ !! NUMBER OF EQUALITY CONSTRAINTS
      integer(ip),intent(out) :: Finish !! 0=LIMIT 1=OPTIM

      integer(ip) :: staflg , faccnt , numcor
      integer(ip) :: con , var , cos , act , ind , len , inc
      integer(ip) :: nnn , typ , des , prv , met
      real(wp) :: Desnor , foldis , cosimp , cornor , quacor , refdis
      real(wp) :: co0 , co1 , co2 , nor
      real(wp) :: cosco2 , cosco1
      real(wp) :: maxdis , norprv
      real(wp) :: val , max , det , ccc , dis
      real(wp) :: fac , del , exc , eps , imp
      real(wp) :: bet , tht
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
            ALLOCATE (cosact(me%Numcon))
            ALLOCATE (varvec(me%Numvar))
            ALLOCATE (varwrk(me%Numvar))
            ALLOCATE (corvec(me%Numcon))
            ALLOCATE (desder(me%Numvar))
            ALLOCATE (desprv(me%Numvar))
            ALLOCATE (varprv(me%Numvar))
            ALLOCATE (convec(me%Numcon+1))
            ALLOCATE (conqua(me%Numcon+1))
            ALLOCATE (concor(me%Numcon+1))
            ! ======================================================================
            ! OPTIMISATION PART
            ! ----------------------------------------------------------------------
            cos = me%Numcon + 1
            des = me%Numcon + 2
            prv = me%Numcon + 3
            faccnt = 0
            me%Conred(1:cos,:) = me%Conder(1:cos,:)
            me%Conred(des,:) = me%Vardes
            desder = me%Vardes
            me%Conred(prv,:) = me%Funvar
            me%Conder(prv,:) = me%Funvar
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3,"OPTIMISATION PART")
            ! ----------------------------------------------------------------------
!           WRITE (STR,'("NUMACT = ",I4)') NUMACT
!           CALL me%ogwrit (3,STR)
            numcor = me%Numact
            concor = me%Actcon
            me%Numact = 0
            me%Conact(1:des) = 0
            spag_nextblock_1 = 2
          CASE (2)
   !        DO COR = 1, NUMCOR
   !            CON = CONCOR(COR)
   !            NAM = CONSTR(CON)
   !            LEN = CONLEN(CON)
   !            WRITE (STR,'("ACT = ",I4,5X,1X,A)') CON, NAM(1:LEN)
   !            CALL me%ogwrit (3,STR)
   !            CALL me%OGINCL (CON)
   !        ENDDO
            ! ======================================================================
            ! ======================================================================
            ! VECTOR OF STEEPEST ASCENT
            ! ----------------------------------------------------------------------
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3,"VECTOR OF STEEPEST ASCENT")
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3,"REMOVE AND INCLUDE CONSTRAINTS:")
            CALL me%ogwrit(3,"")
            CALL me%ogwrit(3," REM  INC CONSTRAINT")
            ! ----------------------------------------------------------------------
            ! REMOVE PASSIVE INEQUALITY CONSTRAINTS
            ! ----------------------------------------------------------------------
            DO act = me%Numact , 1 , -1
               con = me%Actcon(act)
               IF ( me%Conval(con)<=1.0_wp ) CYCLE
               nam = me%Constr(con)
               len = me%Conlen(con)
               WRITE (str,'(I4,5X,1X,A)') con , nam(1:len)
               CALL me%ogwrit(2,str)
               CALL me%ogexcl(act)
               me%Conact(con) = -1
            ENDDO
            ! ----------------------------------------------------------------------
            ! INCLUDE VIOLATED INEQUALITY CONSTRAINTS AND SELECT PASSIVE ONES
            ! SELECT PASSIVE INEQUALITY CONSTRAINTS
            ! ----------------------------------------------------------------------
            DO con = 1 , me%Numcon
               IF ( me%Contyp(con)==-2 ) CYCLE
               IF ( me%Conact(con)>0 ) THEN
               ELSEIF ( me%Contyp(con)==0 ) THEN
                  me%Conact(con) = 0
               ELSEIF ( me%Conval(con)<-1.0_wp ) THEN
                  me%Conact(con) = 0
               ELSEIF ( me%Conval(con)<=+1.0_wp ) THEN
                  me%Conact(con) = 0
               ELSE
                  me%Conact(con) = -1
               ENDIF
            ENDDO
            ! ======================================================================
            nnn = 1
            SPAG_Loop_1_1: DO
               nnn = nnn + 1
               IF ( nnn>999 ) THEN
                  Finish = 0
                  WRITE (str,*) "NNN=" , nnn
                  CALL me%ogwrit(2,str)
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               ! ======================================================================
               ! DERIVATIVES OF MERIT W.R.T. ACTIVE CONSTRAINTS
               ! ----------------------------------------------------------------------
               CALL me%ogrigt(-me%Conred(cos,1:me%Numact),cosact)
               Desnor = sqrt(sum(me%Conred(cos,me%Numact+1:me%Numvar)**2))
               ! ----------------------------------------------------------------------
               ! CONSTRAINT REMOVAL
               ! ----------------------------------------------------------------------
               ind = 0
               exc = -1.0e-12_wp
               max = exc
               DO act = 1 , me%Numact
                  con = me%Actcon(act)
                  IF ( me%Contyp(con)==0 ) CYCLE
                  val = cosact(act)
                  fac = dot_product(me%Conred(con,1:me%Numvar),me%Conred(con,1:me%Numvar))
                  fac = sqrt(fac)
                  val = val*fac
                  IF ( val>=exc ) CYCLE
                  IF ( val>max ) CYCLE
                  max = val
                  ind = act
               ENDDO
               ! ----------------------------------------------------------------------
               IF ( ind/=0 ) THEN
                  con = me%Actcon(ind)
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  WRITE (str,'(I4,5X,3(1X,D10.3),1X,A)') con , Desnor , max , me%Varmax , nam(1:len)
                  CALL me%ogwrit(3,str)
                  CALL me%ogexcl(ind)
                  CYCLE
               ENDIF
               ! ----------------------------------------------------------------------
               ! CONSTRAINT INCLUSION
               ! ----------------------------------------------------------------------
               IF ( Desnor/=0.0_wp ) THEN
                  ! ----------------------------------------------------------------------
                  inc = 0
                  eps = 1.0e-03_wp
                  max = -1.0e10_wp
                  max = 0.0_wp
                  ! ----------------------------------------------------------------------
                  DO con = 1 , me%Numcon
                     IF ( me%Contyp(con)==-2 ) CYCLE
                     IF ( me%Conact(con)/=0 ) CYCLE
                     del = dot_product(me%Conred(con,me%Numact+1:me%Numvar),me%Conred(cos,me%Numact+1:me%Numvar))/Desnor
                     val = abs(del)*me%Varmax
                     IF ( val<eps ) CYCLE
                     fac = dot_product(me%Conred(con,1:me%Numvar),me%Conred(con,1:me%Numvar))
                     fac = sqrt(fac)
                     del = del/fac
                     IF ( del<0.0_wp .AND. del<max ) THEN
                        max = del
                        inc = con
                     ENDIF
                     IF ( me%Contyp(con)/=0 ) CYCLE
                     del = -del
                     IF ( del<0.0_wp .AND. del<max ) THEN
                        max = del
                        inc = con
                     ENDIF
                  ENDDO
                  ! ----------------------------------------------------------------------
                  IF ( inc/=0 ) THEN
                     con = inc
                     nam =me%Constr(con)
                     len = me%Conlen(con)
                     WRITE (str,'(5X,I4,3(1X,D10.3),1X,A)') con , Desnor , max , me%Varmax , nam(1:len)
                     CALL me%ogwrit(3,str)
                     CALL me%ogincl(inc)
                     CYCLE
                  ENDIF
               ENDIF
               ! ----------------------------------------------------------------------
               ! ----------------------------------------------------------------------
               ! MATCHED INEQUALITY CONSTRAINTS + STEEPEST ASCENT VECTOR NORM
               ! ----------------------------------------------------------------------
               CALL me%ogwrit(3,"")
               CALL me%ogwrit(3,"STATUS OF MATCHED INEQUALITY CONSTRAINTS:")
               CALL me%ogwrit(3,"")
               CALL me%ogwrit(3," ACT  PAS CONSTRAINT")
               DO con = 1 , me%Numcon
                  IF ( me%Contyp(con)==-2 ) CYCLE
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  IF ( me%Contyp(con)==0 ) THEN
                  ELSEIF ( me%Conact(con)>0 ) THEN
                     WRITE (str,'(I4,5X,1X,A)') con , nam(1:len)
                     CALL me%ogwrit(3,str)
                  ELSEIF ( me%Conact(con)==0 ) THEN
                     WRITE (str,'(5X,I4,1X,A)') con , nam(1:len)
                     CALL me%ogwrit(3,str)
                  ENDIF
               ENDDO
               CALL me%ogwrit(3,"")
               WRITE (str,'("STEEPEST ASCENT NORM: ",D13.6)') Desnor
               CALL me%ogwrit(3,str)
               WRITE (str,'("MAXIMUM DISTANCE....: ",D13.6)') me%Varmax
               CALL me%ogwrit(3,str)
               ! ======================================================================
               Finish = 0
               ! ======================================================================
               ! IF CONVERGENCE
               ! ----------------------------------------------------------------------
               cosimp = Desnor*me%Varmax
               IF ( abs(cosimp)<=1.0_wp ) THEN
                  foldis = 0.0_wp
                  Finish = 1
                  WRITE (str,'("FINAL...............:",1X,D13.6,'//'11X,1(1X,D10.3),1X,D16.9)') &
                     foldis , cosimp , me%Conval(cos) + cosimp
                  CALL me%ogwrit(2,str)
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               ! ======================================================================
               ! IF CONSTRAINT IS HIT
               ! ----------------------------------------------------------------------
               DO var = 1 , me%Numvar
                  val = me%Conder(cos,var)
                  DO act = 1 , me%Numact
                     con = me%Actcon(act)
                     val = val + me%Conder(con,var)*cosact(act)
                  ENDDO
                  me%Vardes(var) = val
               ENDDO
               ! ----------------------------------------------------------------------
               ind = 0
               dis = 1.0e10_wp
               DO con = 1 , me%Numcon
                  IF ( me%Contyp(con)==-2 ) CYCLE
                  IF ( me%Conact(con)/=-1 ) CYCLE
                  val = dot_product(me%Conred(con,me%Numact+1:me%Numvar),me%Conred(cos,me%Numact+1:me%Numvar))
                  IF ( val==0.0_wp ) CYCLE
                  val = -me%Conval(con)/val*Desnor
                  IF ( val<=0.0_wp ) CYCLE
                  IF ( val>=dis ) CYCLE
                  dis = val
                  ind = con
               ENDDO
               ! ----------------------------------------------------------------------
               IF ( ind/=0 ) THEN
                  val = sqrt(sum((me%Varval-me%Varref+me%Vardes*dis/Desnor)**2))
                  IF ( val>me%Varmax ) ind = 0
               ENDIF
               ! ----------------------------------------------------------------------
               IF ( ind/=0 ) THEN
                  IF ( me%Confix(ind)<=0 ) THEN
                     IF ( val>me%Varmax*1.0e-1_wp ) ind = 0
                  ENDIF
               ENDIF
               ! ----------------------------------------------------------------------
               IF ( ind/=0 ) THEN
                  imp = dis*Desnor
                  con = ind
                  tht = 1.0_wp
                  bet = 0.0_wp
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  WRITE (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,11X,1(1X,D10.3),1X,D16.9,22X,1X,I4,1X,A)') &
                     dis , imp , me%Conval(cos) + cosimp , con , nam(1:len)
                  CALL me%ogwrit(2,str)
                  Varacc = Varacc + dis
                  me%Varval = me%Varval + dis*me%Vardes/Desnor
                  DO con = 1 , me%Numcon + 1
                     val = dot_product(me%Conred(con,me%Numact+1:me%Numvar),me%Conred(cos,me%Numact+1:me%Numvar))
                     me%Conval(con) = me%Conval(con) + val*dis/Desnor
                  ENDDO
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               EXIT SPAG_Loop_1_1
            ENDDO SPAG_Loop_1_1
            SPAG_Loop_1_2: DO
               ! ======================================================================
               ! ----------------------------------------------------------------------
               CALL me%ogrigt(-me%Conred(cos,1:me%Numact),cosact)
               DO var = 1 , me%Numvar
                  val = me%Conder(cos,var)
                  DO act = 1 ,me%Numact
                     con = me%Actcon(act)
                     val = val + me%Conder(con,var)*cosact(act)
                  ENDDO
                  desprv(var) = val
               ENDDO
               Desnor = sqrt(sum(desprv**2))
               WRITE (str,'("DESNOR=",D13.6)') Desnor
!              CALL me%ogwrit (2,STR)
               ! ----------------------------------------------------------------------
               CALL me%ogrigt(-me%Conred(prv,1:me%Numact),cosact)
               DO var = 1 , me%Numvar
                  val = me%Conder(prv,var)
                  DO act = 1 , me%Numact
                     con = me%Actcon(act)
                     val = val + me%Conder(con,var)*cosact(act)
                  ENDDO
                  varprv(var) = val
               ENDDO
               norprv = sqrt(sum(varprv**2))
               WRITE (str,'("NORPRV=",D13.6)') norprv
!              CALL me%ogwrit (2,STR)
               ! ----------------------------------------------------------------------
               CALL me%ogrigt(-me%Conred(des,1:me%Numact),cosact)
               DO var = 1 , me%Numvar
                  val = desder(var)
                  DO act = 1 , me%Numact
                     con = me%Actcon(act)
                     val = val + me%Conder(con,var)*cosact(act)
                  ENDDO
                  me%Vardir(var) = val
               ENDDO
               nor = sqrt(sum(me%Vardir**2))
               WRITE (str,'("NOR=",D13.6)') nor
!              CALL me%ogwrit (2,STR)
               met = me%Optmet
               tht = 1.0_wp
               bet = 0.0_wp
               select case (met)
                case ( 0 ) ! STEEPEST DESCENT METHOD
                case ( 1 ) ! MODIFIED SPECTRAL CONJUGATE GRADIENT METHOD
                  varvec = desprv - varprv
                  IF ( norprv**2>1.0e-12_wp ) THEN
                     tht = -dot_product(me%Vardir,varvec)/norprv**2
                     bet = Desnor**2/norprv**2
                  ENDIF
                case ( 2 ) ! SPECTRAL CONJUGATE GRADIENT METHOD
                  varvec = desprv - varprv
                  varwrk = me%Varref - me%Vargrd
                  val = dot_product(varwrk,varvec)
                  fac = dot_product(me%Vardir,varvec)
                  IF ( abs(val)>1.0e-12_wp .AND. abs(fac)>1.0e-12_wp ) THEN
                     tht = -dot_product(varwrk,varwrk)/val
                     varwrk = -varvec*tht - varwrk
                     bet = dot_product(varwrk,desprv)/fac
                  ENDIF
                case ( 3 ) ! CONJUGATE GRADIENT METHOD
                  IF ( norprv/=0.0_wp ) THEN
                     tht = 1.0_wp
                     bet = Desnor**2/norprv**2
                  ENDIF
               end select
!              WRITE (STR,'("THT=",D13.6)') THT
!              CALL me%ogwrit (3,STR)
!              WRITE (STR,'("BET=",D13.6)') BET
!              CALL me%ogwrit (3,STR)
               ! ----------------------------------------------------------------------
               eps = 1.0e-03_wp
!              WRITE (STR,*) "THT/BET=",THT,BET
!              CALL me%ogwrit (2,STR)
               me%Vardes = tht*desprv + bet*me%Vardir
               Desnor = sqrt(sum(me%Vardes**2))
               nor = Desnor
               DO con = 1 , me%Numcon
                  IF ( me%Contyp(con)==-2 ) CYCLE
                  IF ( me%Conact(con)/=0 ) CYCLE
                  del = dot_product(me%Conder(con,1:me%Numvar),me%Vardes(1:me%Numvar))/nor
                  val = abs(del)*me%Varmax
                  IF ( val<eps ) CYCLE
                  fac = dot_product(me%Conder(con,1:me%Numvar),me%Conder(con,1:me%Numvar))
                  del = del/sqrt(fac)
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  typ = me%Contyp(con)
                  IF ( del<0.0_wp ) THEN
!                    CALL OGINCL (CON)
                     act = me%Conact(con)
!                    WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)') &
!                      CON,ACT,+NOR,VAL,DEL,NAM(1:LEN)
!                    CALL me%ogwrit (2,STR)
                     IF ( act/=0 ) CYCLE SPAG_Loop_1_2
                     bet = 0.0_wp
                     tht = 1.0_wp
                  ENDIF
                  IF ( me%Contyp(con)/=0 ) CYCLE
                  del = -del
                  IF ( del<0.0_wp ) THEN
!                    CALL OGINCL (CON)
                     act = me%Conact(con)
!                    WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)') &
!                      CON,ACT,-NOR,VAL,DEL,NAM(1:LEN)
!                    CALL me%ogwrit (2,STR)
                     IF ( act/=0 ) CYCLE SPAG_Loop_1_2
                     bet = 0.0_wp
                     tht = 1.0_wp
                     del = dot_product(me%Conder(con,1:me%Numvar),me%Conder(cos,1:me%Numvar))
                     val = abs(del)*me%Varmax/Desnor
!                    WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)') &
!                       CON,TYP,-DESNOR,VAL,DEL,NAM(1:LEN)
!                    CALL me%ogwrit (2,STR)
                  ENDIF
               ENDDO
               cosco1 = dot_product(me%Vardes,me%Conder(cos,1:me%Numvar))/Desnor
               IF ( cosco1<0.0_wp ) THEN
                  WRITE (str,*) "COSCO1=" , cosco1
                  CALL me%ogwrit(2,str)
                  bet = 0.0_wp
                  tht = 1.0_wp
               ENDIF
               me%Vardes = tht*desprv + bet*me%Vardir
               Desnor = sqrt(sum(me%Vardes**2))
               EXIT SPAG_Loop_1_2
            ENDDO SPAG_Loop_1_2
            SPAG_Loop_1_3: DO
!              WRITE (STR,*) "THT/BET=",THT,BET
!              CALL me%ogwrit (2,STR)
               ! ======================================================================
               ! SECOND ORDER DERIVATIVES BY NUMERICAL DIFFERENCING
               ! ----------------------------------------------------------------------
               ! ----------------------------------------------------------------------
               CALL me%ogwrit(3,"")
               WRITE (str,'("SECOND ORDER CORRECTION")')
               CALL me%ogwrit(3,str)
               me%Vardes = tht*desprv + bet*me%Vardir
               Desnor = sqrt(sum(me%Vardes**2))
               ! ----------------------------------------------------------------------
               ! MAXIMUM TRAVEL DISTANCE
               ! ----------------------------------------------------------------------
               dis = me%Varmax
               varvec = me%Varval - me%Varref
               co0 = dot_product(varvec,varvec) - dis**2
               IF ( co0>=0.0_wp ) THEN
                  dis = 0.0_wp
               ELSE
                  co1 = dot_product(me%Vardes,varvec)
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
               IF ( me%Senopt>=+1 ) THEN
                  DO con = 1 , me%Numcon
                     IF ( me%Contyp(con)==-2 ) CYCLE
                     act = me%Conact(con)
                     ind = me%Senact(con)
                     IF ( act==-1 .AND. ind==-1 ) CYCLE
                     IF ( act==0 .AND. ind==0 ) CYCLE
                     IF ( act>0 .AND. ind>0 ) CYCLE
                     nam = me%Constr(con)
                     len = me%Conlen(con)
                     WRITE (str,'(I4,1X,I4,1X,I4,1X,A)') con , act , ind , nam(1:len)
                     CALL me%ogwrit(2,str)
                     nnn = 1
                  ENDDO
               ENDIF
               IF ( me%Senopt<=0 .OR. nnn==1 ) THEN
                  fac = me%Varstp/Desnor
                  varvec = me%Varref + me%Vardes*me%Varstp/Desnor
                  CALL me%ogeval(varvec,convec,0,me%Conder)
                  conqua = matmul(me%Conder(1:cos,1:me%Numvar),me%Vardes(1:me%Numvar))
                  conqua = 2.0_wp*(convec-me%Conref-conqua*fac)/fac**2
               ENDIF
               IF ( me%Senopt==-1 ) THEN
                  me%Senqua = conqua
               ELSEIF ( me%Senopt>=+1 .AND. nnn==0 ) THEN
                  conqua = me%Senqua
               ENDIF
               ! ======================================================================
               ! COMPUTE CORRECTION VECTOR
               ! ----------------------------------------------------------------------
               DO act = 1 , me%Numact
                  con = me%Actcon(act)
                  corvec(act) = conqua(con)
               ENDDO
               CALL me%ogleft(corvec,corvec)
               ! ----------------------------------------------------------------------
               cornor = sqrt(sum(corvec(1:me%Numact)**2))*0.5_wp/Desnor/Desnor
               CALL me%ogwrit(3,"")
               WRITE (str,'("STEEPEST ASCENT  NORM: ",D13.6)') Desnor
               CALL me%ogwrit(3,str)
               WRITE (str,'("ACCUMULATED  DISTANCE: ",D13.6)') Varacc
               CALL me%ogwrit(3,str)
               ! ======================================================================
               ! GET EXTREMUM
               ! ----------------------------------------------------------------------
               cosco1 = dot_product(me%Vardes,me%Conder(cos,1:me%Numvar))/Desnor
               cosco2 = conqua(cos) - dot_product(me%Conred(cos,1:me%Numact),corvec(1:me%Numact))
               cosco2 = cosco2*0.5_wp/Desnor/Desnor
               WRITE (str,*) "COSCO2/COSCO1=" , cosco2 , cosco1
               CALL me%ogwrit(3,str)
               IF ( cosco1<0.0_wp ) CALL me%ogwrit(2,str)
               ! ----------------------------------------------------------------------
               foldis = 0.0_wp
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               CALL me%ogwrit(3,"")
               WRITE (str,'(    "STEEPEST ASCENT FOLLOW",'//'  5X,"DISTANCE",'//'  1X,"CORRECTION",'//'  2X,"MERIT_DEL",'//'  6X,"MERIT_VALUE")')
               CALL me%ogwrit(3,str)
               WRITE (str,'("INITIAL.............:",1X,D13.6,'//'  2(1X,D10.3),1X,D16.9)') foldis , quacor , cosimp ,me%Conval(cos) + cosimp
               CALL me%ogwrit(3,str)
               ! ======================================================================
               IF ( cosco2<0.0_wp ) THEN
                  foldis = -0.5_wp*cosco1/cosco2
                  quacor = cornor*foldis*foldis
                  cosimp = foldis*(cosco1+foldis*cosco2)
                  WRITE (str,'("MERIT MAXIMUM.......:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp
                  CALL me%ogwrit(3,str)
                  staflg = 1
               ELSEIF ( cosco2>0.0_wp ) THEN
                  foldis = me%Varmax
                  quacor = cornor*foldis*foldis
                  cosimp = foldis*(cosco1+foldis*cosco2)
                  WRITE (str,'("MERIT CONVEX........:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp
                  CALL me%ogwrit(3,str)
                  staflg = 2
               ELSE
                  foldis = me%Varmax
                  quacor = cornor*foldis*foldis
                  cosimp = foldis*(cosco1+foldis*cosco2)
                  WRITE (str,'("MERIT LINEAR........:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp
                  CALL me%ogwrit(3,str)
                  staflg = 2
               ENDIF
               ! ======================================================================
               ! IF MAXIMUM DISTANCE IS HIT
               ! ----------------------------------------------------------------------
               IF ( foldis>me%Varmax ) THEN
                  foldis = me%Varmax
                  quacor = cornor*foldis*foldis
                  cosimp = foldis*(cosco1+foldis*cosco2)
                  WRITE (str,'("MAXIMUM DISTANCE....:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp
                  CALL me%ogwrit(3,str)
                  staflg = 2
               ENDIF
               ! ======================================================================
               ! IF CONVERGENCE
               ! ----------------------------------------------------------------------
               IF ( abs(cosimp)<=1.0_wp ) THEN
                  WRITE (str,'("FINAL...............:",1X,D13.6,'//'  2(1X,D10.3),1X,D16.9,2D11.3)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp , tht , bet
                  CALL me%ogwrit(2,str)
                  IF ( tht/=1.0_wp .OR. bet/=0.0_wp ) THEN
                     tht = 1.0_wp
                     bet = 0.0_wp
                     CYCLE
                  ENDIF
                  foldis = 0.0_wp
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
                  WRITE (str,'("REMAINING DISTANCE..:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp
                  CALL me%ogwrit(3,str)
                  staflg = 2
               ENDIF
               ! ======================================================================
               ! IF CONSTRAINT IS HIT
               ! ----------------------------------------------------------------------
               ind = 0
               DO con = 1 , me%Numcon
                  IF ( me%Contyp(con)==-2 ) CYCLE
                  IF ( me%Conact(con)/=-1 ) CYCLE
                  co2 = conqua(con) - dot_product(me%Conred(con,1:me%Numact),corvec(1:me%Numact))
                  co1 = dot_product(me%Conder(con,1:me%Numvar),me%Vardes(1:me%Numvar))
                  co0 = me%Conval(con)*2.0_wp
                  IF ( co2/=0.0_wp ) THEN
                     det = co1**2 - co2*co0
                     IF ( det<0.0_wp ) CYCLE
                     det = sqrt(det)
                     val = 1.0e10_wp
                     fac = (-co1+det)/co2
                     IF ( fac>0.0_wp .AND. fac<val ) val = fac
                     fac = (-co1-det)/co2
                     IF ( fac>0.0_wp .AND. fac<val ) val = fac
                  ELSEIF ( co1/=0.0_wp ) THEN
                     val = -co0/co1*0.5_wp
                  ELSE
                     CYCLE
                  ENDIF
                  val = val*Desnor
                  IF ( val>0.0_wp .AND. val<foldis ) THEN
                     foldis = val
                     ind = con
                  ENDIF
               ENDDO
               ! ----------------------------------------------------------------------
               IF ( ind/=0 ) THEN
                  quacor = cornor*foldis*foldis
                  cosimp = foldis*(cosco1+foldis*cosco2)
                  con = ind
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  WRITE (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,2(1X,D10.3),1X,D16.9,1X,I4,1X,A)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp , con , nam(1:len)
                  CALL me%ogwrit(3,str)
                  staflg = 3
               ENDIF
               ! ======================================================================
               ! UPDATE
               ! ----------------------------------------------------------------------
               refdis = foldis
               EXIT SPAG_Loop_1_3
            ENDDO SPAG_Loop_1_3

            SPAG_Loop_1_4: DO
               WRITE (str,'("FINAL...............:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                  foldis , quacor , cosimp , me%Conval(cos) + cosimp
               CALL me%ogwrit(3,str)
               ! ----------------------------------------------------------------------
               fac = foldis/Desnor
               ! ----------------------------------------------------------------------
               ! VARIABLE DELTA
               ! ----------------------------------------------------------------------
               CALL me%ogrigt(corvec,cosact)
               dis = 0.0_wp
               DO var = 1 , me%Numvar
                  val = 0.0_wp
                  DO act = 1 , me%Numact
                     ind = me%Actcon(act)
                     val = val - cosact(act)*me%Conder(ind,var)
                  ENDDO
                  varvec(var) = fac*(me%Vardes(var)+(val*fac*0.5_wp))
                  dis = dis + varvec(var)*varvec(var)
               ENDDO
               dis = sqrt(dis)
               ! ----------------------------------------------------------------------
               WRITE (str,*) "REFDIS=" , refdis
               CALL me%ogwrit(3,str)
               WRITE (str,*) "FOLDIS=" , foldis
               CALL me%ogwrit(3,str)
               WRITE (str,*) "DIS=" , dis
               CALL me%ogwrit(3,str)
               IF ( dis>refdis*1.2_wp .AND. me%Senopt>0 ) THEN
                  faccnt = faccnt + 1
                  IF ( faccnt>=10 ) EXIT SPAG_Loop_1_4
                  foldis = foldis*0.5_wp
                  quacor = cornor*foldis*foldis
                  cosimp = foldis*(cosco1+foldis*cosco2)
                  CYCLE
               ENDIF
               ! ----------------------------------------------------------------------
               ! UPDATE VARIABLES
               ! ----------------------------------------------------------------------
               Varacc = Varacc + foldis
               me%Varval = me%Varval + varvec
               ccc = sqrt(sum((me%Varval-me%Varref)**2)) - me%Varmax**2
               IF ( ccc>=0.0_wp ) THEN
                  WRITE (str,*) "CCC > 0" , ccc
                  CALL me%ogwrit(3,str)
                  staflg = 2
               ENDIF
               ! ======================================================================
               ! MAXIMUM REACHED: NEXT ITERATION
               ! ----------------------------------------------------------------------
               IF ( staflg==1 ) THEN
                  WRITE (str,'("MERIT MAXIMUM.......:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9,2D11.3)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp , tht , bet
                  CALL me%ogwrit(2,str)
                  IF ( me%Senopt>0 ) Finish = 1
                  EXIT SPAG_Loop_1_4
               ENDIF
               ! ======================================================================
               ! MAXIMUM TRAVEL DISTANCE REACHED: NEXT ITERATION
               ! ----------------------------------------------------------------------
               IF ( staflg==2 ) THEN
                  WRITE (str,'("REMAINING DISTANCE..:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9,2D11.3)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp , tht , bet
                  CALL me%ogwrit(2,str)
                  EXIT SPAG_Loop_1_4
               ENDIF
               ! ======================================================================
               ! CONSTRAINT HIT: UPDATE CONSTRAINT + CORRECT
               ! ----------------------------------------------------------------------
               IF ( staflg==3 ) THEN
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  WRITE (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,2(1X,D10.3),1X,D16.9,2D11.3,1X,I4,1X,A)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp , tht , bet , con , nam(1:len)
                  CALL me%ogwrit(2,str)
                  convec = conqua - matmul(me%Conred(1:cos,1:me%Numact),corvec(1:me%Numact))
                  convec = convec*fac*0.5_wp
                  convec = convec + matmul(me%Conder(1:cos,1:me%Numvar),me%Vardes(1:me%Numvar))
                  me%Conval = me%Conval + convec*fac
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               EXIT SPAG_Loop_1_4
            ENDDO SPAG_Loop_1_4
            spag_nextblock_1 = 3
          CASE (3)
            ! ======================================================================
            ! ----------------------------------------------------------------------
            me%Funvar = desprv
            me%Confix = me%Conact(1:me%Numcon)
            IF ( me%Senopt==-1 ) me%Senact = me%Conact(1:me%Numcon)
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

   SUBROUTINE ogpwri(me,Objval,Numvio,Convio)

      !! WRITE OPTIMIZATION LOG IN PYGMO FORMAT
      !!
      !! 2023/01/25 | W. MARTENS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Objval !! OBJECTIVE VALUE
      integer(ip),intent(in) :: Numvio !! NUMBER OF VIOLATED CONSTRAINTS
      real(wp),intent(in) :: Convio !! TOTAL CONSTRAINT VIOLATION

      CHARACTER(len=2) :: feas

      IF ( me%Verbos==0 ) RETURN
      ! Print header
      IF ( me%Fevals==0 ) CALL me%ogpwri_start()
      ! Increase counter for cost function evaluations
      me%Fevals = me%Fevals + 1
      ! Every 50 lines print the column names.
      IF ( mod(real(me%Fevals-1.0_wp)/real(me%Verbos),50.0_wp)==0.0_wp ) &
         WRITE (me%Loglup,'(A10,A15,A15,A15,A2)') "objevals:" , "objval:" , "violated:" , "viol. norm:"
      IF ( me%Verbos/=0 .AND. mod(me%Fevals,me%Verbos)==0.0_wp ) THEN
         IF ( Convio>0.0_wp ) THEN
            feas = " i"
         ELSE
            feas = "  "
         ENDIF

         ! Write the log line (different format depending on violation size)
         IF ( Convio==0.0_wp ) THEN
            WRITE (me%Loglup,'(I10,F15.4,I15,I15,A2)') me%Fevals , Objval , Numvio , int(Convio) , feas
         ELSEIF ( Convio>1.0e-3_wp ) THEN
            WRITE (me%Loglup,'(I10,F15.4,I15,F15.6,A2)') me%Fevals , Objval , Numvio , Convio , feas
         ELSE
            WRITE (me%Loglup,'(I10,F15.4,I15,E15.6,A2)') me%Fevals , Objval , Numvio , Convio , feas
         ENDIF
      ENDIF

      ! Write final summary
      IF ( me%Pygfla/=0 ) CALL me%ogpwri_end(Objval,Numvio,Convio)

   END SUBROUTINE ogpwri

   SUBROUTINE ogpwri_end(me,Objval,Numvio,Convio)

      !! WRITE OPTIMIZATION END RESULT IN PYGMO FORMAT
      !!
      !! 2023/01/25 | W. MARTENS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Objval !! OBJECTIVE VALUE
      real(wp),intent(in) :: Convio !! TOTAL CONSTRAINT VIOLATION
      integer(ip),intent(in) :: Numvio !! NUMBER OF VIOLATED CONSTRAINTS

      IF ( me%Pygfla==0 ) RETURN

      ! Write termination message
      WRITE (me%Loglup,'("")')
      WRITE (me%Loglup,'("Final values after iteration        ", I10:)') me%Numite
      WRITE (me%Loglup,'("Final objective value:              ", F10.4)') Objval
      WRITE (me%Loglup,'("Final constraint violation:         ", F10.4)') Convio
      WRITE (me%Loglup,'("Final num. of violated constraints: ",I10)') Numvio

      select case (me%Pygfla)
       case ( 1 ); WRITE (me%Loglup,'("Successful termination: Optimal solution found.")')
       case ( 2 ); WRITE (me%Loglup,'("Successful termination: Constraints matched.")')
       case ( 3 ); WRITE (me%Loglup,'("Not converged.")')
       case ( 4 ); WRITE (me%Loglup,'("Problem appears infeasible.")')
      end select
      WRITE (me%Loglup,'("")')

   END SUBROUTINE ogpwri_end

   SUBROUTINE ogpwri_start(me)

      !! WRITE OPTIMIZATION LOG IN PYGMO FORMAT
      !!
      !! 2023/01/25 | W. MARTENS | NEW

      class(optgra),intent(inout) :: me

      WRITE (me%Loglup,'("OPTGRA plugin for pagmo/pygmo:")')
      select case (me%Varder) ! DERIVATIVES COMPUTATION MODE
       case ( 0 );      WRITE (me%Loglup,'("")') ! VALUES ONLY
       case ( 1, -1 );  WRITE (me%Loglup,'("    User-defined gradients")')
       case ( 2 );      WRITE (me%Loglup,'("    Numerical gradients by double differencing")')
       case ( 3 );      WRITE (me%Loglup,'("    Numerical gradients by single differencing")')
      end select

      select case (me%Optmet)
       case ( 0 ); WRITE (me%Loglup,'("    Steepest descent method")')
       case ( 1 ); WRITE (me%Loglup,'("    Modified spectral conjugate gradient method")')
       case ( 2 ); WRITE (me%Loglup,'("    Spectral conjugate gradient method")')
       case ( 3 ); WRITE (me%Loglup,'("    Conjugate gradient method")')
      end select

      WRITE (me%Loglup,'("")')

   END SUBROUTINE ogpwri_start

   SUBROUTINE ogrigt(me,Actinp,Actout)

      !! RIGHT-MULTIPLIES VECTOR LOWER TRIANGULAR MATRIX OBTAINED BY REDUCTION
      !! AND SUBSEQUENT INVERSION OF DERIVATIVES OF ACTIVE CONSTRAINTS
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Actinp(me%Numcon) !! VECTOR INITAL
      real(wp),intent(inout) :: Actout(me%Numcon) !! VECTOR FINAL (MAY BE SAME AS ACTINP)

      integer(ip) :: row , col , act
      real(wp) :: val

      DO col = me%Numact , 1 , -1
         val = Actinp(col)
         DO act = me%Numact , col + 1 , -1
            row = me%Actcon(act)
            val = val - me%Conred(row,col)*Actout(act)
         ENDDO
         row = me%Actcon(col)
         Actout(col) = val/me%Conred(row,col)
      ENDDO

   END SUBROUTINE ogrigt

   SUBROUTINE ogsens(me,Consta,Concon,Convar,Varcon,Varvar)

      !! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      integer(ip),intent(out) :: Consta(me%Numcon) !! CONSTRAINT STATUS (0=PAS 1=ACT)
      real(wp),intent(out) :: Concon(me%Numcon+1,me%Numcon) !! SENSITIVITY OF CONTRAINTS+MERIT W.R.T. ACTIVE CONSTRAINTS
      real(wp),intent(out) :: Convar(me%Numcon+1,me%Numvar) !! SENSITIVITY OF CONTRAINTS+MERIT W.R.T. PARAMETERS
      real(wp),intent(out) :: Varcon(me%Numvar,me%Numcon) !! SENSITIVITY OF VARIABLES W.R.T. ACTIVE CONSTRAINTS
      real(wp),intent(out) :: Varvar(me%Numvar,me%Numvar) !! SENSITIVITY OF VARIABLES W.R.T. PARAMETERS
      !! -> NOT SCALED

      real(wp) :: val , sca
      integer(ip) :: var , con , act , par , ind , typ

      ! CONVERGED
      Consta = 0
      DO act = 1 , me%Numact
         con = me%Actcon(act)
         Consta(con) = 1
      ENDDO

      ! SENSITIVITY OF CONTRAINTS W.R.T. ACTIVE CONSTRAINTS
      Concon = 0.0_wp
      DO con = 1 , me%Numcon + 1
         IF ( me%Conact(con)>0 ) Concon(con,con) = 1.0_wp
         IF ( me%Conact(con)>0 ) CYCLE
         me%Conref = me%Conred(con,1:me%Numact)
         CALL me%ogrigt(me%Conref,me%Conref)
         DO act = 1 , me%Numact
            ind = me%Actcon(act)
            Concon(con,ind) = -me%Conref(act)
         ENDDO
      ENDDO

      ! SENSITIVITY OF CONSTRAINTS W.R.T. PARAMETERS
      Convar = 0.0_wp
      DO con = 1 , me%Numcon + 1
         IF ( me%Conact(con)>0 ) CYCLE
         DO var = 1 , me%Numvar
            IF ( me%Vartyp(var)==0 ) CYCLE
            val = me%Sender(con,var)
            DO act = 1 , me%Numact
               ind = me%Actcon(act)
               val = val + Concon(con,ind)*me%Sender(ind,var)
            ENDDO
            Convar(con,var) = val
         ENDDO
      ENDDO

      ! SENSITIVITY OF VARIABLES W.R.T. ACTIVE CONSTRAINTS
      Varcon = 0.0_wp
      DO var = 1 , me%Numvar
         IF ( me%Vartyp(var)/=0 ) CYCLE
         DO act = 1 , me%Numact
            con = me%Actcon(act)
            me%Conref(act) = me%Conder(con,var)
         ENDDO
         CALL me%ogleft(me%Conref,me%Conref)
         CALL me%ogrigt(me%Conref,me%Conref)
         DO act = 1 , me%Numact
            con = me%Actcon(act)
            Varcon(var,con) = -me%Conref(act)
         ENDDO
      ENDDO

      ! SENSITIVITY OF VARIABLES W.R.T. PARAMETERS
      Varvar = 0.0_wp
      DO par = 1 , me%Numvar
         Varvar(par,par) = 1.0_wp
         IF ( me%Vartyp(par)/=1 ) CYCLE
         DO var = 1 , me%Numvar
            IF ( me%Vartyp(var)/=0 ) CYCLE
            val = 0.0_wp
            DO act = 1 , me%Numact
               con = me%Actcon(act)
               val = val + Varcon(var,con)*me%Sender(con,par)
            ENDDO
            Varvar(var,par) = val
         ENDDO
      ENDDO

      ! DESCALE SENSITIVITY
      DO con = 1 , me%Numcon + 1
         typ = me%Contyp(con)
         sca = me%Consca(con)
         IF ( typ<0 ) sca = -sca
         Convar(con,1:me%Numvar) = Convar(con,1:me%Numvar)*sca
         Concon(con,1:me%Numcon) = Concon(con,1:me%Numcon)*sca
         IF ( con>me%Numcon ) CYCLE
         Varcon(1:me%Numvar,con) = Varcon(1:me%Numvar,con)/sca
         Concon(1:me%Numcon+1,con) = Concon(1:me%Numcon+1,con)/sca
      ENDDO

      DO var = 1 , me%Numvar
         sca = me%Varsca(var)
         Varcon(var,1:me%Numcon) = Varcon(var,1:me%Numcon)*sca
         Varvar(var,1:me%Numvar) = Varvar(var,1:me%Numvar)*sca
         Convar(1:me%Numcon+1,var) = Convar(1:me%Numcon+1,var)/sca
         Varvar(1:me%Numvar,var) = Varvar(1:me%Numvar,var)/sca
      ENDDO

   END SUBROUTINE ogsens

   SUBROUTINE ogssst(me,Varsen,Quasen,Consen,Actsen,Dersen,Actsav,Consav,Redsav,Dersav,Actnum)

      !! NEAR-LINEAR OPTIMISATION TOOL SENSITIVITY ANALYSIS
      !!
      !! Function to get sensitivity state data, necessary for serialization.
      !! Do not use this directly except in serialization routines
      !!
      !! 2021/07/19 | M. von Looz | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in)    :: Varsen(me%Numvar) !! STORED VARIABLES VALUE
      real(wp),intent(in)    :: Quasen(me%Numcon+1) !! STORED CONSTRAINTS CORRECTION VECTOR
      real(wp),intent(in)    :: Consen(me%Numcon+1) !! STORED CONSTRAINTS VALUE
      integer(ip),intent(in) :: Actsen(me%Numcon+1) !! STORED CONSTRAINTS ACTIVE
      real(wp),intent(in)    :: Dersen(me%Numcon+1,me%Numvar) !! STORED DERIVATIVE
      integer(ip),intent(in) :: Actsav(me%Numcon+1) !! STORED ACTIVE CONSTRAINTS
      integer(ip),intent(in) :: Consav(me%Numcon+4) !! STORED ACTIVE CONSTRAINTS
      real(wp),intent(in)    :: Redsav(me%Numcon+3,me%Numvar) !! STORED DERIVATIVE
      real(wp),intent(in)    :: Dersav(me%Numcon+3,me%Numvar) !! STORED DERIVATIVE

      integer(ip) :: Actnum
      integer(ip) :: var , con

      ! Variable values saved for sensitivity
      me%Numact = Actnum

      DO var = 1 , me%Numvar
         me%Senvar(var) = Varsen(var)
      ENDDO

      DO con = 1 , me%Numcon + 1
         me%Senqua(con) = Quasen(con)
         me%Sencon(con) = Consen(con)
         me%Senact(con) = Actsen(con)
         DO var = 1 , me%Numvar
            me%Sender(con,var) = Dersen(con,var)
         ENDDO
      ENDDO

      ! Temporary status saved of which constraints are active
      DO con = 1 , me%Numcon + 1
         me%Actcon(con) = Actsav(con)
      ENDDO

      DO con = 1 , me%Numcon + 4
         me%Conact(con) = Consav(con)
      ENDDO

      DO con = 1 , me%Numcon + 3
         DO var = 1 , me%Numvar
            me%Conred(con,var) = Redsav(con,var)
         ENDDO
      ENDDO

      DO con = 1 , me%Numcon
         DO var = 1 , me%Numvar
            me%Conder(con,var) = Dersav(con,var)
         ENDDO
      ENDDO
   END SUBROUTINE ogssst

   SUBROUTINE ogwrit(me,Lev,Str)

      !! Write a meessage to the log.
      !!
      !! 2014/07/29 | J. SCHOENMAEKERS | NEW
      !! 2025/02/09 | J. Williams | added trim()

      class(optgra),intent(inout) :: me
      CHARACTER(len=*),intent(in) :: Str !! string to print
      integer(ip),intent(in) :: Lev !! only print if `Loglev` is >= this

      IF ( Lev<=me%Loglev ) THEN
         WRITE (me%Loglun,'(A)') trim(Str)
         FLUSH (me%Loglun)
      ENDIF

   END SUBROUTINE ogwrit

!****************************************************************************************************
end module optgra_module
!****************************************************************************************************
