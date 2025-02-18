!****************************************************************************************************
!>
!  Near-linear optimization tool tailored for s/c trajectory design.
!
!  This is a modernization of the original Fortran 77 code from: [pyoptgra](https://github.com/esa/pyoptgra).
!  It has been extensively refactored. The old common blocks were removed, and all data is now
!  encapsulated in a single [[optgra]] class.
!
!### History
!  * Original code by J. Schoenmaekers et al. (ESA) from [pyoptgra](https://github.com/esa/pyoptgra).
!  * Jacob Williams, Feb 2025 : created new version.
!
!@todo finish removing the `spag_nextblock_1` loop thing which I don't like.

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

   integer,parameter,public :: optgra_wp = wp   !! Working real precision

   integer(ip), parameter :: name_len = 80 !! max string length for var names
   integer(ip), parameter :: str_len = 256 !! max string length for str prints

   type,public :: optgra
      !! Main class for OPTGRA algorithm.
      !!
      !! The main methods are [[optgra:initialize]], [[optgra:solve]], and [[optgra:destroy]].

      private

      integer(ip) :: numvar = 0 !! number of variables
      real(wp),      dimension(:  ), allocatable :: varval
      integer(ip),   dimension(:  ), allocatable :: vartyp
      real(wp),      dimension(:  ), allocatable :: varsca
      character(len=name_len),    dimension(:  ), allocatable :: varstr
      integer(ip),   dimension(:  ), allocatable :: varlen
      real(wp),      dimension(:  ), allocatable :: varref
      real(wp),      dimension(:  ), allocatable :: vardes
      real(wp),      dimension(:  ), allocatable :: vargrd
      real(wp),      dimension(:  ), allocatable :: vardir
      real(wp),      dimension(:  ), allocatable :: funvar
      real(wp),      dimension(:  ), allocatable :: senvar

      integer(ip) :: numcon = 0 !! number of constraints
      real(wp),      dimension(:  ), allocatable :: conval
      integer(ip),   dimension(:  ), allocatable :: contyp
      integer(ip),   dimension(:  ), allocatable :: conpri
      real(wp),      dimension(:  ), allocatable :: consca
      character(len=name_len),    dimension(:  ), allocatable :: constr
      integer(ip),   dimension(:  ), allocatable :: conlen
      real(wp),      dimension(:  ), allocatable :: conref
      real(wp),      dimension(:  ), allocatable :: senqua
      real(wp),      dimension(:  ), allocatable :: sencon
      real(wp),      dimension(:  ), allocatable :: sendel
      integer(ip),   dimension(:  ), allocatable :: senact

      integer(ip)  :: optmet = 2 !! optimization method
      integer(ip)  :: maxite = 10 !! maximum number of iterations
      integer(ip)  :: corite = 10
      integer(ip)  :: optite = 10
      integer(ip)  :: divite = 10
      integer(ip)  :: cnvite = 10
      real(wp)     :: Varmax = 10.0_wp !! maximum distance per iteration
      real(wp)     :: varsnd = 1.0_wp !! perturbation for 2nd order derivatives
      real(wp)     :: varstp = 1.0_wp

      integer(ip) :: varder = 1 !! derivatives computation mode
      real(wp),      dimension(:  ), allocatable :: varper

      integer(ip) :: loglun = output_unit  !! log file unit
      integer(ip) :: loglev = 1  !! log level

      integer(ip) :: loglup = output_unit  !! pygmo log file unit
      integer(ip) :: verbos = 0  !! pygmo verbosity
      integer(ip) :: fevals = 0  !! pygmo: number of const fun evals
      integer(ip) :: pygfla = 0  !! pygmo: flag indicating status of optimization
      integer(ip) :: numite = 0  !! number of iterations

      integer(ip) :: matlev = 0

      integer(ip) :: tablun = output_unit !! logical unit for writing table
      integer(ip) :: tablev = 0 !! level of tab

      integer(ip) :: senopt = 0 !! sensitivity optimization mode

      integer(ip) :: numact = 0
      integer(ip),   dimension(:  ), allocatable :: actcon
      integer(ip),   dimension(:  ), allocatable :: confix
      integer(ip),   dimension(:  ), allocatable :: conact
      real(wp),      dimension(:,:), allocatable :: conder
      real(wp),      dimension(:,:), allocatable :: conred
      real(wp),      dimension(:,:), allocatable :: sender
      !integer(ip) :: CONVER = 0   ! not used ?
      !integer(ip),   DIMENSION(:  ), allocatable :: CONOPT  ! not used ?

      procedure(calval_f),pointer :: calval => null() !! function for values
      procedure(calder_f),pointer :: calder => null() !! function for derivatives

   contains

      private

      procedure,public :: initialize !! set up the problem
      procedure,public :: solve => ogexec !! solve the problem
      procedure,public :: destroy => ogclos !! free memory when finished

      ! are these intented to be user callable?
      procedure,public :: ogsens
      procedure,public :: ogssst

      ! private methods:
      procedure :: ogrigt
      procedure :: ogincl
      procedure :: ogexcl
      procedure :: ogeval
      procedure :: ogcorr
      procedure :: ogopti
      procedure :: ogleft
      procedure :: ogpwri
      procedure :: ogwrit
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

   subroutine initialize(me,Numvar,Numcon,&
                         Calval,Calder,Delcon,Conpri,Consca,Constr,Conlen,&
                         Contyp,Varder,Varper,Varmax,Varsnd,&
                         Maxite,Itecor,Iteopt,Itediv,Itecnv,&
                         Loglup,Verbos,Senopt,Varsca,&
                         Varstr,Varlen,Vartyp,Loglun,Loglev,Matlev,&
                         Tablun,Tablev,Optmet)

      !! Initialize the class. This should be the first routine called.
      !!
      !! Note: this is a combination of the routines from the original code:
      !! oginit, ogcdel, ogcpri, ogcsca, ogcstr, ogctyp, ogderi, ogdist,
      !! ogiter, ogomet, ogplog, ogsopt, ogvsca, ogvstr, ogvtyp, ogwlog,
      !! ogwmat, ogwtab
      !!
      !!@note Could make some of these optional and keep the default values.

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
                                                   !!
                                                   !! lower constraint priorities are fulfilled earlier.
                                                   !! During the initial constraint correction phase,
                                                   !! only constraints with a priority at most k are
                                                   !! considered in iteration k.
      real(wp),intent(in) :: Consca(Numcon+1)  !! * CONSTRAINTS CONVER THRESHOLD (1:NUMCON)
                                               !! * MERIT       CONVER THRESHOLD (1+NUMCON)
      character(len=name_len),intent(in) :: Constr(Numcon+1) !! CONIABLES NAME STRING
      integer(ip),intent(in) :: Conlen(Numcon+1) !! CONIABLES NAME LENGTH
      integer(ip),intent(in) :: Contyp(Numcon+1) !! CONSTRAINTS TYPE (1:NUMCON)
                                                 !!
                                                 !!  *  1 = GTE (positive inequality constraints)
                                                 !!  * -1 = LTE (inequality constraints that should be negative)
                                                 !!  *  0 = EQU (equality constraints)
                                                 !!  * -2 = DERIVED DATA (unenforced constraints)
                                                 !!
                                                 !! MERIT TYPE (1+NUMCON)
                                                 !!
                                                 !!  *  1 = MAX (maximization problems)
                                                 !!  * -1 = MIN (minimization problems)
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
      integer(ip),intent(in) :: Maxite !! maximum number of iterations
      integer(ip),intent(in) :: Itecor !! number of constraint correction iterations
                                       !! in the beginning. If no feasible solution is
                                       !! found within that many iterations, Optgra aborts
      integer(ip),intent(in) :: Iteopt !! some other iter input ?
      integer(ip),intent(in) :: Itediv !! some other iter input ?
      integer(ip),intent(in) :: Itecnv !! some other iter input ?
      integer(ip),intent(in) :: Loglup !! LOGICAL UNIT FOR WRITING PYGMO LOG
      integer(ip),intent(in) :: Verbos !! VERBOSITY LEVEL:
                                       !!
                                       !!  * 0=NO OUTPUT
                                       !!  * 1 OUTPUT EVERY ITERATION
                                       !!  * 2 OUTPUT EVERY 2ND ITERATION
                                       !!  * N OUTPUT EVERY NTH ITERATION
      integer(ip),intent(in) :: Senopt  !! sensitivity optimization mode
                                        !!
                                        !!  *  0: NO
                                        !!  * -1: INITIALIZATION
                                        !!  * +1: WITH CONSTRAINT CALCULATION
                                        !!  * +2: WITH CONSTRAINT BIAS
                                        !!  * +3: WITH CONSTRAINT CALC / NO OPTIM
                                        !!  * +4: WITH CONSTRAINT BIAS / NO OPTIM
      real(wp),intent(in) :: Varsca(Numvar) !! VARIABLES SCALE FACTOR
      character(len=name_len),intent(in) :: Varstr(Numvar) !! VARIABLES NAME STRING
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
      integer(ip),intent(in) :: Tablun !! logical unit for writing table
      integer(ip),intent(in) :: Tablev !! LEVEL OF TAB
                                       !!
                                       !!  * 0=NO OUTPUT
                                       !!  * 1<ALL
      integer(ip),intent(in) :: Optmet !! OPTIMIZATION METHOD:
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

      allocate (me%Varval(me%Numvar))
      allocate (me%Vartyp(me%Numvar))
      allocate (me%Varsca(me%Numvar))
      allocate (me%Varstr(me%Numvar))
      allocate (me%Varlen(me%Numvar))
      allocate (me%Varref(me%Numvar))
      allocate (me%Vardes(me%Numvar))
      allocate (me%Vargrd(me%Numvar))
      allocate (me%Vardir(me%Numvar))
      allocate (me%Funvar(me%Numvar))
      allocate (me%Senvar(me%Numvar))

      do var = 1 , me%Numvar
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
      enddo

      ! CONSTRAINTS
      me%Numcon = Numcon

      allocate (me%Conval(me%Numcon+1))
      allocate (me%Contyp(me%Numcon+1))
      allocate (me%Conpri(me%Numcon+1))
      allocate (me%Consca(me%Numcon+1))
      allocate (me%Constr(me%Numcon+1))
      allocate (me%Conlen(me%Numcon+1))
      allocate (me%Conref(me%Numcon+1))
      allocate (me%Senqua(me%Numcon+1))
      allocate (me%Sencon(me%Numcon+1))
      allocate (me%Sendel(me%Numcon+1))
      allocate (me%Senact(me%Numcon+1))

      do con = 1 , me%Numcon + 1
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
      enddo

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
      allocate (me%Varper(me%Numvar))
      do var = 1 , me%Numvar
         me%Varper(var) = 1.0e-03_wp
      enddo

      me%Loglun = output_unit     ! LOG FILE
      me%Loglev = 1
      me%Loglup = output_unit     ! PYGMO LOG FILE
      me%Loglev = 0
      me%Matlev = 0               ! MATLAB CONSOLE
      me%Tablun = output_unit     ! TABLE FILE
      me%Tablev = 0

      me%Senopt = 0    ! LINEAR OPTIMIZATION MODE

      ! WORKING VECTORS
      allocate (me%Actcon(me%Numcon+1))
      allocate (me%Confix(me%Numcon))
      allocate (me%Conact(me%Numcon+4))
      allocate (me%Conder(me%Numcon+3,me%Numvar))
      allocate (me%Conred(me%Numcon+3,me%Numvar))
      allocate (me%Sender(me%Numcon+3,me%Numvar))
      ! ALLOCATE (me%Conopt(me%Numcon+1))
      me%Numact = 0
      me%Actcon = 0
      me%Conact = 0
      me%Confix = 0
      me%Conder = 0.0_wp
      me%Conred = 0.0_wp
      !me%Conopt = 0

      ! NOTE: should the last element also be set? (not done in original)  ?????
      do con = 1 , me%Numcon
         me%Sendel(con) = Delcon(con)
      enddo

      do con = 1 , me%Numcon + 1
         me%Conpri(con) = Conpri(con)
         me%Consca(con) = Consca(con)
         me%Constr(con) = Constr(con)
         me%Conlen(con) = min(Conlen(con),name_len)
         me%Contyp(con) = Contyp(con)
      enddo

      me%Varder = Varder

      do var = 1 , me%Numvar
         me%Varper(var) = Varper(var)
         me%Varsca(var) = Varsca(var)
         me%Varstr(var) = Varstr(var)
         me%Varlen(var) = min(Varlen(var),name_len)
         me%Vartyp(var) = Vartyp(var)
      enddo

      me%Varmax = Varmax
      me%Varsnd = Varsnd

      me%Maxite = Maxite
      me%Corite = Itecor
      me%Optite = Iteopt
      me%Divite = Itediv
      me%Cnvite = Itecnv
      if ( me%Corite>me%Maxite ) me%Corite = me%Maxite
      if ( me%Optite>me%Maxite ) me%Optite = me%Maxite
      if ( me%Divite>me%Corite ) me%Divite = me%Corite
      if ( me%Cnvite>me%Optite ) me%Cnvite = me%Optite

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

   subroutine ogclos(me)
      !! DEALLOCATION OF ARRAYS
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me

      ! VARIABLES
      if (allocated(me%Varval)) deallocate (me%Varval)
      if (allocated(me%Vartyp)) deallocate (me%Vartyp)
      if (allocated(me%Varsca)) deallocate (me%Varsca)
      if (allocated(me%Varstr)) deallocate (me%Varstr)
      if (allocated(me%Varlen)) deallocate (me%Varlen)
      if (allocated(me%Varref)) deallocate (me%Varref)
      if (allocated(me%Vardes)) deallocate (me%Vardes)
      if (allocated(me%Vargrd)) deallocate (me%Vargrd)
      if (allocated(me%Vardir)) deallocate (me%Vardir)
      if (allocated(me%Funvar)) deallocate (me%Funvar)
      if (allocated(me%Senvar)) deallocate (me%Senvar)

      ! CONSTRAINTS
      if (allocated(me%Conval)) deallocate (me%Conval)
      if (allocated(me%Contyp)) deallocate (me%Contyp)
      if (allocated(me%Conpri)) deallocate (me%Conpri)
      if (allocated(me%Consca)) deallocate (me%Consca)
      if (allocated(me%Constr)) deallocate (me%Constr)
      if (allocated(me%Conlen)) deallocate (me%Conlen)
      if (allocated(me%Conref)) deallocate (me%Conref)
      if (allocated(me%Senqua)) deallocate (me%Senqua)
      if (allocated(me%Sencon)) deallocate (me%Sencon)
      if (allocated(me%Sendel)) deallocate (me%Sendel)
      if (allocated(me%Senact)) deallocate (me%Senact)

      ! DERIVATIVES
      if (allocated(me%Varper)) deallocate (me%Varper)

      ! WORKING VECTORS
      if (allocated(me%Actcon)) deallocate (me%Actcon)
      if (allocated(me%Confix)) deallocate (me%Confix)
      if (allocated(me%Conact)) deallocate (me%Conact)
      if (allocated(me%Conder)) deallocate (me%Conder)
      if (allocated(me%Conred)) deallocate (me%Conred)
      if (allocated(me%Sender)) deallocate (me%Sender)
      !if (allocated(me%Conopt)) DEALLOCATE (me%Conopt)

   end subroutine ogclos

   subroutine ogcorr(me,Varacc,Finish,Toterr,Norerr,error)

      !! CORRECTION PART
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp) :: Varacc
      integer(ip) :: Finish
      real(wp) :: Toterr
      real(wp) :: Norerr
      logical,intent(out) :: error !! if there was a fatal error

      integer(ip) :: coritr , numfff , minpri , maxpri , curpri
      real(wp) :: cornor , foldis , cstval , conerr
      real(wp) :: corinv , varvio , conmax , normax
      integer(ip) :: conind , norind , inelop , maxitr
      integer(ip) :: con , var , act , ind , len , cos , stp
      integer(ip) :: typ , cor , pri , vio , fff
      real(wp) :: val , fac , upr , del , co2 , co1 , co0 , de2 , dis
      real(wp) :: eps , err , dlt , sca , dif
      real(wp) :: exc
      character(len=str_len) :: str
      character(len=name_len) :: nam
      real(wp) , dimension(:) , allocatable :: cosact
      real(wp) , dimension(:) , allocatable :: varvec
      real(wp) , dimension(:) , allocatable :: varsav
      real(wp) , dimension(:) , allocatable :: varcor
      real(wp) , dimension(:) , allocatable :: corvec
      real(wp) , dimension(:) , allocatable :: consav
      integer(ip) , dimension(:) , allocatable :: conttt
      integer(ip) , dimension(:) , allocatable :: concor
      integer(ip) , dimension(:) , allocatable :: coninc
      integer(ip) , dimension(:) , allocatable :: conhit
      integer(ip) , dimension(:) , allocatable :: fffcon
      integer(ip) , dimension(:) , allocatable :: prisav
      integer :: spag_nextblock_1

      error = .false.

      ! initialize   - JW : these could also be in the class.
      !                     that would save the allocate/reallocate
      !                     step every time this routine is called.
      allocate (cosact(me%Numvar))
      allocate (varvec(me%Numvar))
      allocate (varsav(me%Numvar))
      allocate (varcor(me%Numvar))
      allocate (corvec(me%Numvar))
      allocate (consav(me%Numcon+1))
      allocate (conttt(me%Numcon+1))
      allocate (concor(me%Numcon))
      allocate (coninc(me%Numcon))
      allocate (conhit(me%Numcon))
      allocate (fffcon(me%Numcon))
      allocate (prisav(me%Numcon))

      ! ======================================================================
      ! CORRECTION PART
      ! ----------------------------------------------------------------------
      coritr = 0
      maxitr = 10
      if ( me%Senopt/=0 ) maxitr = 1
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
      do con = 1 , me%Numcon
         if ( me%Contyp(con)==-2 ) cycle
         minpri = min(minpri,me%Conpri(con))
         maxpri = max(maxpri,me%Conpri(con))
      enddo
      ind = 0
      if ( numfff>0 ) ind = 1
      if ( me%Senopt>0 ) ind = 1
      minpri = minpri - ind
      call me%ogwrit(3,"")
      call me%ogwrit(3,"PRIORITISE CONSTRAINTS")
      call me%ogwrit(3,"")
      if ( me%Senopt<=0 ) then
         do fff = 1 , numfff
            con = fffcon(fff)
            nam = me%Constr(con)
            len = me%Conlen(con)
            typ = me%Contyp(con)
            write (str,'("PRI",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
            call me%ogwrit(3,str)
            me%Conpri(con) = minpri
         enddo
      endif
      if ( me%Senopt>0 ) then
         do con = 1 , me%Numcon
            if ( me%Contyp(con)==-2 ) cycle
            if ( me%Senact(con)<=0 ) cycle
            nam = me%Constr(con)
            len = me%Conlen(con)
            typ = me%Contyp(con)
            write (str,'("PRI",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
            call me%ogwrit(3,str)
            me%Conpri(con) = minpri
         enddo
      endif
      call me%ogwrit(3,"")
      spag_nextblock_1 = 2

      main: do
         select case (spag_nextblock_1)
          case (2)
            ! ======================================================================
            ! Evaluation loop
            ! ----------------------------------------------------------------------
            write (str,'("CORITR=",I2,1X,I2)') coritr , maxitr
            call me%ogwrit(3,str)
            call me%ogwrit(3,"")
            ! ======================================================================
            ! Inequality loop
            ! ----------------------------------------------------------------------
            inelop = 2
            if ( numfff>0 ) inelop = 1
            if ( me%Senopt>0 ) inelop = 1
            varsav = me%Varval
            consav = me%Conval
            spag_nextblock_1 = 3
          case (3)
            ! ----------------------------------------------------------------------
            write (str,'("INELOP=",I2)') inelop
            call me%ogwrit(3,str)
            call me%ogwrit(3,"")
            ! ----------------------------------------------------------------------
            me%Varval = varsav
            me%Conval = consav
            me%Contyp = conttt
            if ( inelop==1 ) then
               do fff = 1 , numfff
                  con = fffcon(fff)
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  typ = me%Contyp(con)
                  write (str,'("TYP",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
                  call me%ogwrit(3,str)
                  me%Contyp(con) = 0
               enddo
            endif
            if ( inelop==1 .and.me%Senopt>0 ) then
               do con = 1 , me%Numcon
                  if ( me%Contyp(con)==-2 ) cycle
                  if ( me%Senact(con)<=0 ) cycle
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  typ = me%Contyp(con)
                  write (str,'("FIX",1X,I4,1X,I4,1X,A)') con , typ , nam(1:len)
                  call me%ogwrit(2,str)
                  me%Contyp(con) = 0
               enddo
            endif
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
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               typ = me%Contyp(con)
               val = me%Conval(con)
               if ( val<-1.0_wp ) then
                  err = abs(val)
               elseif ( typ/=0 ) then
                  err = 0.0_wp
               elseif ( val>1.0_wp ) then
                  err = abs(val)
               else
                  err = 0.0_wp
               endif
               conerr = conerr + err
               if ( err>conmax ) then
                  conind = con
                  conmax = err
               endif
               fac = 0.0_wp
               do var = 1 , me%Numvar
                  if ( me%Vartyp(var)==1 ) cycle
                  fac = fac + me%Conder(con,var)**2
               enddo
               fac = sqrt(fac)
               if ( err==0.0_wp ) then
               elseif ( fac/=0.0_wp ) then
                  err = err/fac
               else
                  call me%ogwrit(0,"")
                  call me%ogwrit(0,"ERROR: CONSTRAINT CAN NOT BE SATISFIED")
                  write (str,'("CON/VAL= ",I5,1X,D13.6)') con , val
                  call me%ogwrit(0,str)
                  Finish = 0
                  call done()
                  exit main
               endif
               Norerr = Norerr + err
               if ( err>normax ) then
                  norind = con
                  normax = err
               endif
            enddo
            ! ----------------------------------------------------------------------
            Toterr = conerr
            call me%ogwrit(3,"")
            write (str,'("NUMFFF/TOTERR/NORERR/COSVAL=",I4,3(1X,D13.6))') &
                     numfff , Toterr , Norerr , me%Conval(me%Numcon+1)
            call me%ogwrit(2,str)
            call me%ogwrit(3,"")
            write (str,'("MAXIM TOTAL ERROR.: ",D13.6,I6)') conmax , conind
            call me%ogwrit(3,str)
            write (str,'("MAXIM NORM  ERROR.: ",D13.6,I6)') normax , norind
            call me%ogwrit(3,str)
            ! ----------------------------------------------------------------------
            if ( me%Senopt>0 .and. coritr==maxitr ) then
               if ( cstval==0.0_wp ) then
                  Finish = 1
               else
                  Finish = 0
               endif
               exit main
            endif
            if ( coritr==0 .and. conerr==0.0_wp ) then
               me%Numact = numfff
               me%Actcon = fffcon
               Finish = 1
               exit main
            elseif ( coritr/=0 .and. conerr==0.0_wp ) then
               Finish = 1
               exit main
            elseif ( coritr==maxitr ) then
               Finish = 0
               call me%ogwrit(3,"")
               write (str,'("CORITR=",I2)') coritr
               call me%ogwrit(3,str)
               exit main
            else
               Finish = 0
               coritr = coritr + 1
            endif

            ! ======================================================================
            ! Priority loop
            ! ----------------------------------------------------------------------
            curpri = minpri
            spag_nextblock_1 = 4
          case (4)
            ! ----------------------------------------------------------------------
            write (str,'("CURPRI=",I2,1X,I2)') curpri , maxpri
            call me%ogwrit(3,str)
            call me%ogwrit(3,"")
            ! ======================================================================
            ! MINIMUM NORM CORRECTION
            ! ----------------------------------------------------------------------
            call me%ogwrit(3,"")
            write (str,'("CORRECTION OF CONSTRAINTS")')
            call me%ogwrit(3,str)
            call me%ogwrit(3,"")
            write (str,*) "INELOP/CURPRI=" , inelop , curpri
            call me%ogwrit(3,str)
            ! ----------------------------------------------------------------------
            conerr = 0.0_wp
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               conhit(con) = 0
               if ( me%Conval(con)<-eps ) then
                  conerr = conerr + abs(me%Conval(con))
               elseif ( me%Contyp(con)/=0 ) then
               elseif ( me%Conval(con)>+eps ) then
                  conerr = conerr + abs(me%Conval(con))
               endif
            enddo
            ! ----------------------------------------------------------------------
            call me%ogwrit(3,"")
            write (str,'("LINEAR ERROR.: ",D13.6)') conerr
            call me%ogwrit(3,str)
            ! ----------------------------------------------------------------------
            call me%ogwrit(3,"")
            write (str,'(" ACT  PAS  MOV",'//' " COST___VAL COST___GRD",'//' " DIST___DEL CONSTRAINT")')
            call me%ogwrit(3,str)
            spag_nextblock_1 = 5
          case (5)
            ! ======================================================================
            ! Move loop
            ! ----------------------------------------------------------------------
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               coninc(con) = 0
               pri = me%Conpri(con)
               val = me%Conval(con)
               typ = me%Contyp(con)
               act = me%Conact(con)
               cor = concor(con)
               if ( act>0 ) cycle
               if ( pri>curpri ) then
                  me%Conact(con) = -1
               elseif ( val<-dlt ) then
                  me%Conact(con) = -1
               elseif ( val>+dlt ) then
                  me%Conact(con) = -1
               else
                  me%Conact(con) = 0
               endif
               if ( pri>curpri ) then
                  concor(con) = 0
               elseif ( val<-dlt ) then
                  concor(con) = -1
               elseif ( typ/=0 ) then
                  concor(con) = 0
               elseif ( val>+dlt ) then
                  concor(con) = +1
               else
                  concor(con) = 0
               endif
!              IF (ACT /= CONACT(CON) .OR. COR /= CONCOR(CON)) THEN
!                  NAM = CONSTR(CON)
!                  LEN = CONLEN(CON)
!                  WRITE (STR,'(5X,5X,I4,23X,D10.3,1X,A,4I4)') &
!                     CON, CONVAL(CON),NAM(1:LEN), &
!                     CONACT(CON),CONCOR(CON), ACT, COR
!                  CALL me%ogwrit (3,STR)
!              ENDIF
            enddo
            ! ======================================================================
            ! STEEPEST ASCENT VECTOR
            ! ======================================================================
            ! MERIT VALUE AND DERIVATIVES
            ! ----------------------------------------------------------------------
            cstval = 0.0_wp
            varvec = 0.0_wp
            me%Conred(vio,:) = 0.0_wp
            ! ----------------------------------------------------------------------
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               if ( me%Conpri(con)>curpri ) cycle
               if ( concor(con)==0 ) cycle
               fac = concor(con)
               cstval = cstval - me%Conval(con)*fac
               me%Conred(vio,:) = me%Conred(vio,:) - me%Conred(con,:)*fac
               varvec = varvec - me%Conder(con,:)*fac
            enddo
            inner: do
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
               call me%ogrigt(corvec,cosact)
               ! ----------------------------------------------------------------------
               ! CONSTRAINT REMOVAL
               ! ----------------------------------------------------------------------
               ind = 0
               exc = 1.0e-12_wp
               upr = exc
               do act = 1 , me%Numact
                  con = me%Actcon(act)
                  if ( me%Contyp(con)==0 ) cycle
                  val = cosact(act)
                  if ( val<=exc ) cycle
                  if ( val<upr ) cycle
!                 IF (VAL >= UPR .AND. UPR > 0.0_wp) CYCLE
                  upr = val
                  ind = act
               enddo
               ! ----------------------------------------------------------------------
               if ( ind/=0 ) then
                  con = me%Actcon(ind)
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  write (str,'(5X,I4,5X,3(1X,D10.3),1X,A)') con , cstval , cornor , upr , nam(1:len)
                  call me%ogwrit(3,str)
                  call me%ogexcl(ind,error)
                  if (error) return
                  if ( coninc(con)>=5 ) then
                     write (str,'("OGCORR-WARNING: CONSTRAINT INCLUDED")')
                     call me%ogwrit(1,str)
                     write (str,'("CON/INC/UPR=",2I4,1X,D10.3)') con , coninc(con) , upr
                     call me%ogwrit(1,str)
                  endif
                  if ( coninc(con)>=20 ) then
                     Finish = 0
                     call done()
                     exit main
                  endif
                  cycle
               endif
               ! ----------------------------------------------------------------------
               ! NORMALISE STEEPEST ASCEND VECTOR
               ! ----------------------------------------------------------------------
               if ( cstval<-cornor*varvio ) exit inner
               ! ----------------------------------------------------------------------
               corinv = 1.0_wp/cornor
               corvec(me%Numact+1:me%Numvar) = corvec(me%Numact+1:me%Numvar)*corinv
               ! ----------------------------------------------------------------------
               ! CONSTRAINT INCLUSION
               ! ----------------------------------------------------------------------
               ind = 0
               upr = 0.0_wp
               ! ----------------------------------------------------------------------
               do con = 1 , me%Numcon
                  if ( me%Contyp(con)==-2 ) cycle
                  if ( me%Conpri(con)>curpri ) cycle
                  if ( me%Conact(con)/=0 ) cycle
                  del = dot_product(me%Conred(con,me%Numact+1:me%Numvar),me%Conred(vio,me%Numact+1:me%Numvar))
                  val = abs(del)*me%Varmax/cornor
                  if ( val<eps ) cycle
                  fac = dot_product(me%Conred(con,1:me%Numvar),me%Conred(con,1:me%Numvar))
                  del = del/sqrt(fac)
                  if ( del<upr ) then
                     upr = del
                     ind = con
                  endif
                  if ( me%Contyp(con)/=0 ) cycle
                  if ( concor(con)/=0 ) cycle
                  del = -del
                  if ( del<upr ) then
                     upr = del
                     ind = con
                  endif
               enddo
               ! ----------------------------------------------------------------------
               if ( ind/=0 ) then
                  con = ind
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  write (str,'(I4,5X,5X,3(1X,D10.3),1X,A)') con , cstval , cornor , upr , nam(1:len)
                  call me%ogwrit(3,str)
                  call me%ogincl(ind)
                  coninc(con) = coninc(con) + 1
                  cycle inner
               endif
               exit inner
            enddo inner
            ! ----------------------------------------------------------------------
            do var = 1 , me%Numvar
               val = varvec(var)
               do act = 1 , me%Numact
                  con = me%Actcon(act)
                  val = val - me%Conder(con,var)*cosact(act)
               enddo
               varvec(var) = val*corinv
            enddo
            ! ----------------------------------------------------------------------
            varcor = me%Varval - me%Varref
            co2 = dot_product(varvec,varvec)
            co1 = dot_product(varvec,varcor)*0.5_wp
            co0 = dot_product(varcor,varcor) - me%Varmax**2
            de2 = co1**2 - co2*co0
            if ( de2>=0.0_wp .and. co2/=0.0_wp ) then
               dis = (sqrt(de2)-co1)/co2
            else
               dis = 0.0_wp
            endif
            ! ----------------------------------------------------------------------
            do var = 1 , me%Numvar
               fac = varvec(var)
               if ( fac==0.0_wp ) cycle
               dif = me%Varval(var) - me%Varref(var)
               sca = me%Varmax*1.0e-0_wp   ! JW : why is this multiplied by 1 ?
               val = (dif+sca)/fac
               fac = (dif-sca)/fac
               if ( fac>val ) val = fac
               if ( val<dis ) dis = val
            enddo
            dis = max(dis, 0.0_wp)
            ! ----------------------------------------------------------------------
            foldis = dis
            ! ======================================================================
            ! OPTIMISE DIRETION OF STEPPEST ASCENT
            ! ======================================================================
            if ( cstval==0.0_wp ) then
               ! ----------------------------------------------------------------------
!              WRITE (STR,'("CNV=",3(1X,D10.3))') CSTVAL/CORNOR/VARVIO
!              CALL me%ogwrit (3,STR)
               write (str,'(4X,5X,5X,3(1X,D10.3))') cstval , cornor , foldis
               call me%ogwrit(3,str)
               call me%ogwrit(3,"")
               if ( curpri>=maxpri ) then
                  write (str,'("MAXPRI=",I3)') maxpri
                  call me%ogwrit(3,str)
                  ! ======================================================================
                  ! MATCHED INEQUALITY CONSTRAINTS + MINIMUM CORRECTION NORM
                  ! ----------------------------------------------------------------------
                  call me%ogwrit(3,"")
                  call me%ogwrit(3,"STATUS OF CONSTRAINTS:")
                  call me%ogwrit(3,"")
                  call me%ogwrit(3," ACT  PAS  NON COST___VAL CONSTRAINT")
                  do con = 1 , me%Numcon
                     if ( me%Contyp(con)==-2 ) cycle
                     nam = me%Constr(con)
                     len = me%Conlen(con)
                     val = me%Conval(con)
                     if ( me%Conact(con)>0 ) then
                        write (str,'( I4,5X,6X,D10.3,1X,A)') con , val , nam(1:len)
                        call me%ogwrit(3,str)
                     elseif ( me%Conact(con)==0 ) then
                        write (str,'( 5X,I4,6X,D10.3,1X,A)') con , val , nam(1:len)
                        call me%ogwrit(3,str)
                     elseif ( me%Conact(con)<0 ) then
                        write (str,'(10X,I4,1X,D10.3,1X,A)') con , val , nam(1:len)
                        call me%ogwrit(3,str)
                     endif
                  enddo
                  ! ======================================================================
                  if ( me%Senopt<=0 ) call me%ogeval(me%Varval,me%Conval,0,me%Conder)
                  spag_nextblock_1 = 2
                  cycle main
               else
                  write (str,'("CURPRI=",I3)') curpri
                  call me%ogwrit(3,str)
                  curpri = curpri + 1
                  spag_nextblock_1 = 4
                  cycle main
               endif
               ! ----------------------------------------------------------------------
            elseif ( cstval<-cornor*varvio ) then
               ! ----------------------------------------------------------------------
               write (str,'("CNV=",3(1X,D10.3))') cstval/cornor/varvio
               call me%ogwrit(2,str)
               call me%ogwrit(3,"")
               write (str,'("CSTVAL=",3D10.3)') cstval , cornor , varvio
               call me%ogwrit(3,str)
               if ( inelop==1 ) then
                  write (str,'("INELOP=",I3)') inelop
                  call me%ogwrit(3,str)
                  inelop = inelop + 1
                  spag_nextblock_1 = 3
                  cycle main
               else
                  Finish = 0
                  exit main
               endif
               ! ----------------------------------------------------------------------
            endif
            ! ======================================================================
            ! IF CONSTRAINT IS HIT
            ! ----------------------------------------------------------------------
            ind = 0
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               if ( me%Conact(con)/=-1 ) cycle
               val = dot_product(me%Conred(con,me%Numact+1:me%Numvar),corvec(me%Numact+1:me%Numvar))
               if ( val==0.0_wp ) cycle
               val = -me%Conval(con)/val
               if ( val<=0.0_wp ) cycle
               if ( val>=foldis ) cycle
               foldis = val
               ind = con
            enddo
            ! ======================================================================
            ! UPDATE VARIABLES, CONSTRAINTS AND COST FUNCTION
            ! ----------------------------------------------------------------------
            varvec = varvec*foldis
            ! ----------------------------------------------------------------------
            Varacc = Varacc + foldis
            me%Varval = me%Varval + varvec
            ! ----------------------------------------------------------------------
            do con = 1 , me%Numcon + 1
               val = dot_product(corvec(me%Numact+1:me%Numvar),me%Conred(con,me%Numact+1:me%Numvar))
               me%Conval(con) = me%Conval(con) + val*foldis
            enddo
            ! ----------------------------------------------------------------------
            cstval = cstval + foldis*cornor
            ! ======================================================================
            ! MAXIMUM TRAVEL DISTANCE REACHED: NEXT ITERATION
            ! ----------------------------------------------------------------------
            if ( ind==0 ) then
               write (str,'("CNV=",3(1X,D10.3))') cstval/cornor/varvio
               call me%ogwrit(3,str)
               write (str,'(4X,5X,5X,3(1X,D10.3))') cstval , cornor , foldis
               call me%ogwrit(3,str)
               if ( inelop==1 ) then
                  write (str,'("INELOP=",I3)') inelop
                  call me%ogwrit(3,str)
                  inelop = inelop + 1
                  spag_nextblock_1 = 3
                  cycle main
               else
                  me%Numact = 0
                  exit main
               endif
            endif
            ! ======================================================================
            ! CONSTRAINT HIT: UPDATE CONSTRAINTS + CORRECT
            ! ----------------------------------------------------------------------
            con = ind
            nam = me%Constr(con)
            len = me%Conlen(con)
            val = me%Conval(con)
            write (str,'(5X,5X,I4,3(1X,D10.3),1X,A)') con , cstval , cornor , foldis , nam(1:len)
            call me%ogwrit(3,str)
            if ( conhit(con)>=20 ) then
               write (str,'("OGCORR-WARNING: CONSTRAINT HIT")')
               call me%ogwrit(1,str)
               write (str,'("CON/HIT=",2I4)') con , conhit(con)
               call me%ogwrit(1,str)
               Finish = 0
               call done()
               exit main
            endif
            ! ----------------------------------------------------------------------
            conhit(con) = conhit(con) + 1
            spag_nextblock_1 = 5
            cycle main

         end select

      enddo main

      ! ======================================================================
      ! MATCHED INEQUALITY CONSTRAINTS + MINIMUM CORRECTION NORM
      ! ----------------------------------------------------------------------
      call me%ogwrit(3,"")
      write (str,'("CSTVAL:",D13.6)') cstval
      call me%ogwrit(3,str)
      call me%ogwrit(3,"")
      call me%ogwrit(3,"STATUS OF CONSTRAINTS:")
      call me%ogwrit(3,"")
      call me%ogwrit(3," ACT  PAS  NON COST___VAL CONSTRAINT")
      do con = 1 , me%Numcon
         if ( me%Contyp(con)==-2 ) cycle
         nam = me%Constr(con)
         len = me%Conlen(con)
         val = me%Conval(con)
         if ( me%Conact(con)>0 ) then
            write (str,'( I4,5X,6X,D10.3,1X,A)') con , val , nam(1:len)
            call me%ogwrit(3,str)
         elseif ( me%Conact(con)==0 ) then
            write (str,'( 5X,I4,6X,D10.3,1X,A)') con , val , nam(1:len)
            call me%ogwrit(3,str)
         elseif ( me%Conact(con)<0 ) then
            write (str,'(10X,I4,1X,D10.3,1X,A)') con , val , nam(1:len)
            call me%ogwrit(3,str)
         endif
      enddo
      ! ----------------------------------------------------------------------
      call me%ogwrit(3,"")
      call me%ogwrit(3,"STATUS OF VIOLATED CONSTRAINTS:")
      call me%ogwrit(3,"")
      call me%ogwrit(3," CON COST___VAL CONSTRAINT")
      conerr = 0.0_wp
      do con = 1 , me%Numcon
         if ( me%Contyp(con)==-2 ) cycle
         nam = me%Constr(con)
         len = me%Conlen(con)
         val = me%Conval(con)
         if ( val<-dlt ) then
            conerr = conerr + abs(val)
            write (str,'( I4,D11.3,1X,A)') con , val , nam(1:len)
            call me%ogwrit(3,str)
         elseif ( me%Contyp(con)/=0 ) then
         elseif ( val>dlt ) then
            conerr = conerr + abs(val)
            write (str,'( I4,D11.3,1X,A)') con , val , nam(1:len)
            call me%ogwrit(3,str)
         endif
      enddo
      ! ----------------------------------------------------------------------
      call me%ogwrit(3,"")
      write (str,'("LINEAR ERROR.: ",D13.6)') conerr
      call me%ogwrit(3,str)
      call done()

   contains

      subroutine done()   ! 7

         me%Contyp = conttt
         me%Conpri = prisav

         deallocate (cosact)
         deallocate (varvec)
         deallocate (varsav)
         deallocate (varcor)
         deallocate (corvec)
         deallocate (consav)
         deallocate (conttt)
         deallocate (concor)
         deallocate (coninc)
         deallocate (conhit)
         deallocate (fffcon)
         deallocate (prisav)

      end subroutine done

   end subroutine ogcorr

   subroutine ogeval(me,Valvar,Valcon,Varder,Dercon)

      !! COMPUTES SCALED CONTRAINTS+MERIT AND DERIVATIVES
      !! FROM     SCALED VARIABLES
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Valvar(me%Numvar) !! SCALED VARIABLES
      real(wp),intent(out) :: Valcon(me%Numcon+1) !! SCALED CONTRAINTS+MERIT AND DERIVATIVES
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
      character(len=name_len) :: nam
      character(len=str_len) :: str

      !real(wp) :: ggg(4,4) , bbb(4) , vvv(4)  ! JW : not sure what this was for
      real(wp) :: objval
      real(wp) , dimension(:) , allocatable :: varvec
      real(wp) , dimension(:) , allocatable :: convec

      allocate (varvec(me%Numvar))
      allocate (convec(me%Numcon+1))

      ! ======================================================================
      ! GENERAL
      ! ----------------------------------------------------------------------
      write (str,'()')
      call me%ogwrit(3,str)
      select case (Varder)
       case ( 0 );    write (str,'("COMPUTE RESULTS")')
       case ( 1, -1); write (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES USER DEFINED")')
       case ( 2 );    write (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES BY DOUBLE DIFFERENCING")')
       case ( 3 );    write (str,'("COMPUTE RESULTS",'//'   " AND DERIVATIVES BY SINGLE DIFFERENCING")')
      end select
      call me%ogwrit(3,str)

      ! ======================================================================
      ! WRITE VARIABLES
      ! ----------------------------------------------------------------------
      write (str,'()')
      call me%ogwrit(3,str)
      write (str,'("VARIABLES NOT SCALED:")')
      call me%ogwrit(3,str)
      write (str,'()')
      call me%ogwrit(3,str)
      do var = 1 , me%Numvar
         val = Valvar(var)
         sca = me%Varsca(var)
         cod = me%Vartyp(var)
         select case (cod)
          case ( 0 ); typ = "FRE"
          case ( 1 ); typ = "PAR"
         end select
         nam = me%Varstr(var)
         len = me%Varlen(var)
         write (str,'("VAR/VAL/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') &
            var , val*sca , sca , typ , nam(1:len)
         call me%ogwrit(3,str)
      enddo

      ! ======================================================================
      ! DE-SCALE VARIABLES
      ! ----------------------------------------------------------------------
      do var = 1 , me%Numvar
         varvec(var) = Valvar(var)*me%Varsca(var)
      enddo

      ! ======================================================================
      ! GET RESULTS
      ! GET DERIVATIVES IF USER DEFINED
      ! ----------------------------------------------------------------------
      select case (Varder)
       case ( 0 );     call me%calval(varvec,Valcon,0)
       case ( 1, -1 ); call me%calval(varvec,Valcon,1)
       case ( 2 );     call me%calval(varvec,Valcon,1)
       case ( 3 );     call me%calval(varvec,Valcon,1)
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
      do con = 1 , me%Numcon + 1
         convec(con) = Valcon(con)
         sca = me%Consca(con)
         cod = me%Contyp(con)
         if ( cod==-1 ) sca = -sca
         Valcon(con) = Valcon(con)/sca
      enddo

      ! ======================================================================
      ! WRITE RESULTS
      ! ----------------------------------------------------------------------
      write (str,'()')
      call me%ogwrit(3,str)
      write (str,'("RESULTS NOT SCALED:")')
      call me%ogwrit(3,str)
      write (str,'()')
      call me%ogwrit(3,str)
      conerr = 0.0_wp  ! total constraint error (scaled to constr. threshod)
      convio = 0.0_wp  ! total constaint error norm (unscaled)
      ind = 0          ! index of largest constraint violation
      fac = 0.0_wp     ! value of largest constraint violation
      numvio = 0       ! number of violated constraints
      do con = 1 , me%Numcon + 1
         val = Valcon(con)
         sca = me%Consca(con)
         cod = me%Contyp(con)
         sta = "   "
         err = 0.0_wp
         if ( cod==-1 ) sca = -sca
         if ( cod==-2 ) typ = "DER"
         if ( cod==-1 .and. con<=me%Numcon ) typ = "LTE"
         if ( cod==0 .and. con<=me%Numcon ) typ = "EQU"
         if ( cod==1 .and. con<=me%Numcon ) typ = "GTE"
         if ( cod==-1 .and. con>me%Numcon ) typ = "MIN"
         if ( cod==1 .and. con>me%Numcon ) typ = "MAX"
         if ( cod==0 .and. con<=me%Numcon .and. abs(val)>1.0_wp ) then
            sta = "VIO"
            err = abs(val)
            numvio = numvio + 1
         endif
         if ( cod/=0 .and. con<=me%Numcon .and. cod/=-2 .and. -val>1.0_wp ) then
            sta = "VIO"
            err = abs(val)
            numvio = numvio + 1
         endif
         conerr = conerr + err
         convio = convio + (err*sca)**2
         if ( err>fac ) ind = con
         if ( err>fac ) fac = err
         nam = me%Constr(con)
         len = me%Conlen(con)
         write (str,'("CON/VAL/SCA/TYP/STA/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A3,1X,A)') &
            con , val*sca , sca , typ , sta , nam(1:len)
         call me%ogwrit(3,str)
      enddo
      write (str,'()')
      call me%ogwrit(3,str)
      write (str,'("CONSTRAINT ERROR.:",2(1X,D13.6),I6)') conerr , fac , ind
      call me%ogwrit(3,str)
      write (str,'()')
      call me%ogwrit(3,str)
      ! write pygmo-style log output
      objval = -Valcon(me%Numcon+1)
      convio = sqrt(convio)
      call me%ogpwri(objval,numvio,convio)

      ! ======================================================================
      ! NO DERIVATIVES
      ! ----------------------------------------------------------------------
      if ( Varder==0 ) then
         return
      elseif ( Varder==1 .or. Varder==-1 ) then
         call me%calder(varvec,convec,Dercon)
      endif

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
      write (str,'()')
      call me%ogwrit(3,str)
      write (str,'("DERIVATIVES SCALED:")')
      call me%ogwrit(3,str)
      write (str,'()')
      call me%ogwrit(3,str)

      do var = 1 , me%Numvar

         ! ----------------------------------------------------------------------
         ! WRITE VARIABLE
         ! ----------------------------------------------------------------------
         val = Valvar(var)
         sca = me%Varsca(var)
         cod = me%Vartyp(var)
         if ( cod==0 ) typ = "FRE"
         if ( cod==1 ) typ = "PAR"
         nam = me%Varstr(var)
         len = me%Varlen(var)
         write (str,'("VAR/VAL/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') &
                  var , val*sca , sca , typ , nam(1:len)
         call me%ogwrit(4,str)
         write (str,'()')
         call me%ogwrit(4,str)

         select case (Varder)
          case ( 2 )  ! DERIVATIVES BY DOUBLE DIFFERENCING
            per = me%Varper(var)
            sav = varvec(var)
            varvec(var) = sav + per
            call me%calval(varvec,Dercon(1:,var),0)  ! JW added : here
            varvec(var) = sav - per
            call me%calval(varvec,convec,0)
            fac = 0.5_wp/per
            do con = 1 , me%Numcon + 1
               Dercon(con,var) = (Dercon(con,var)-convec(con))*fac
            enddo
            varvec(var) = sav
          case ( 3 )  ! DERIVATIVES BY SINGLE DIFFERENCING
            per = me%Varper(var)
            sav = varvec(var)
            varvec(var) = sav + per
            call me%calval(varvec,Dercon(1:,var),0)  ! JW added : here
            fac = 1.0_wp/per
            do con = 1 , me%Numcon + 1
               Dercon(con,var) = (Dercon(con,var)-convec(con))*fac
            enddo
            varvec(var) = sav
         end select

         ! ----------------------------------------------------------------------
         ! SCALE DERIVATIVES
         ! ----------------------------------------------------------------------
         do con = 1 , me%Numcon + 1
            fac = me%Varsca(var)/me%Consca(con)
            if ( me%Contyp(con)==-1 ) fac = -fac
            Dercon(con,var) = Dercon(con,var)*fac
         enddo

         ! ----------------------------------------------------------------------
         ! WRITE DERIVATIVES
         ! ======================================================================
         do con = 1 , me%Numcon + 1
            der = Dercon(con,var)
            if ( der/=0.0_wp ) then
               sca = me%Consca(con)
               cod = me%Contyp(con)
               if ( cod==-1 ) sca = -sca
               if ( cod==-2 ) typ = "DER"
               if ( cod==-1 .and. con<=me%Numcon ) typ = "LTE"
               if ( cod==0 .and. con<=me%Numcon ) typ = "EQU"
               if ( cod==1 .and. con<=me%Numcon ) typ = "GTE"
               if ( cod==-1 .and. con>me%Numcon ) typ = "MIN"
               if ( cod==1 .and. con>me%Numcon ) typ = "MAX"
               nam = me%Constr(con)
               len = me%Conlen(con)
               write (str,'("CON/DER/SCA/TYP/NAM=",'//'  I5,D14.6,D9.1,1X,A3,1X,A)') &
                  con , der*sca/me%Varsca(var) , sca , typ , nam(1:len)
               call me%ogwrit(4,str)
            endif
         enddo
         write (str,'()')
         call me%ogwrit(4,str)

      enddo

      deallocate (varvec)
      deallocate (convec)

   end subroutine ogeval

   subroutine ogexcl(me,Exc,error)

      !! REMOVE CONSTRAINT TO ACTIVE SET AND REDUCES DERIVATIVES
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      integer(ip),intent(in) :: Exc !! CONSTRAINT TO BE REMOVED
                                    !! SEQUENCE NUMBER IN ACTIVE LIST
      logical,intent(out) :: error !! if there was a fatal error (constraints singular)

      real(wp) :: val , bet , gam
      integer(ip) :: row , col , act , con
      character(len=str_len) :: str

      error = .false.

      ! ======================================================================
      ! ADJUST LIST OF ACTIVE CONSTRAINTS
      ! ----------------------------------------------------------------------
      con = me%Actcon(Exc)
      me%Conact(con) = 0
      me%Numact = me%Numact - 1
      do act = Exc , me%Numact
         con = me%Actcon(act+1)
         me%Actcon(act) = con
         me%Conact(con) = me%Conact(con) - 1
      enddo
      ! ======================================================================
      ! REDUCE FOR SUBSEQUENT CONSTRAINTS
      ! ----------------------------------------------------------------------
      do act = Exc , me%Numact
         con = me%Actcon(act)
         val = 0.0_wp
         do col = act , act + 1
            val = val + me%Conred(con,col)**2
         enddo
         val = sqrt(val)
         if ( me%Conred(con,act)>0.0_wp ) val = -val
         if ( abs(val)<1.0e-15_wp ) then
            write (me%Loglun,*) "OGEXCL-ERROR: CONSTRAINTS SINGULAR"
            call me%ogwrit(2,str)
            write (me%Loglun,*) "VAL=" , val
            call me%ogwrit(2,str)
            error = .true.
            return   ! fatal error
         endif
         me%Conred(con,act) = me%Conred(con,act) - val
         bet = 1.0_wp/(val*me%Conred(con,act))
         do row = 1 , me%Numcon + 3
            if ( me%Conact(row)>act .or. me%Conact(row)<=0 ) then
               gam = 0.0_wp
               do col = act , act + 1
                  if ( me%Conred(row,col)/=0.0_wp ) gam = gam + me%Conred(row,col)*me%Conred(con,col)
               enddo
               if ( gam/=0.0_wp ) then
                  gam = gam*bet
                  do col = act , act + 1
                     me%Conred(row,col) = me%Conred(row,col) + me%Conred(con,col)*gam
                  enddo
               endif
            endif
         enddo
         me%Conred(con,act) = val
         do col = act + 1 , act + 1
            me%Conred(con,col) = 0.0_wp
         enddo
      enddo
      ! ======================================================================
   end subroutine ogexcl

   subroutine ogexec(me,Valvar,Valcon,Finopt,Finite)
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
                                        !!  *  1 =     MATCHED &     OPTIMAL
                                        !!  *  2 =     MATCHED & NOT OPTIMAL
                                        !!  *  3 = NOT MATCHED & NOT OPTIMAL
                                        !!  *  4 = NOT FEASIBL & NOT OPTIMAL
                                        !!  * -1 = Fatal error (constraints singular)
      integer(ip),intent(out) :: Finite !! ?

      integer(ip) :: finish , itecor , iteopt
      integer(ip) :: var , con , typ , len , num , numvio
      real(wp) :: val , sca , red , der , fac , old , convio
      character(len=str_len) :: str
      character(len=name_len) :: nam
      integer(ip) :: numequ , itediv , itecnv
      real(wp) :: varacc , cosnew , cosold , varsav , meamer
      real(wp) :: conerr , desnor , norerr , meaerr
      real(wp) , dimension(:) , allocatable :: varsum
      real(wp) , dimension(:) , allocatable :: varcor
      real(wp) , dimension(:) , allocatable :: concor
      real(wp) , dimension(:,:) , allocatable :: conder_tmp !! JW: created to avoid "array temporary" warning
      logical :: error

      ! initialize:
      allocate (varsum(me%Numvar))
      allocate (varcor(me%Numvar))
      allocate (concor(me%Numcon+1))
      allocate (conder_tmp(me%Numcon+1,me%Numvar))

      ! ======================================================================
      ! GENERAL
      ! ----------------------------------------------------------------------
      ! LOGLEV = 2
      call me%ogwrit(2,"")
      call me%ogwrit(2,"OPTGRA START")
      call me%ogwrit(2,"")
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
      do var = 1 , me%Numvar
         varsum(var) = 0.0_wp
         varcor(var) = 0.0_wp
      enddo
      ! ======================================================================
      ! EQUALTIY CONSTRAINTS IN ACTIVE SET
      ! ----------------------------------------------------------------------
      me%Numact = 0
      do con = 1 , me%Numcon
!        IF (CONTYP(CON) == -2) CYCLE
         nam = me%Constr(con)
         len = me%Conlen(con)
         write (str,*) "CON/PRI=" , con , me%Conpri(con) , " " , nam(1:len)
         call me%ogwrit(3,str)
         me%Conact(con) = 0
         if ( me%Consca(con)>=1.0e+09_wp ) me%Contyp(con) = -2
         if ( me%Contyp(con)==0 ) then
         elseif ( me%Contyp(con)==-2 ) then
            me%Conact(con) = -2
         endif
      enddo
      numequ = me%Numact
      me%Conact(me%Numcon+1) = -3
      me%Conact(me%Numcon+2) = -3
      ! ======================================================================
      ! SCALE VARIABLES
      ! ----------------------------------------------------------------------
      do var = 1 , me%Numvar
         me%Varval(var) = Valvar(var)/me%Varsca(var)
      enddo
      ! ======================================================================
      ! HEADER FOR TABLE
      ! ----------------------------------------------------------------------
      if ( me%Tablev>=1 ) write (me%Tablun,'("ITER",1X,"OPT",1X,*(1X,I10))') (var,var=1,me%Numvar) , (con,con=1,me%Numcon)

      main: do
         inner: do
            ! ======================================================================
            ! IF (NUMITE >= 52) MATLEV = 3
            ! IF (NUMITE >= 55) MATLEV = 2
            ! ======================================================================
            ! NEW ITERATION
            ! ----------------------------------------------------------------------
            if ( me%Numite>=me%Corite .and. itecor==0 ) then
               Finopt = 3
               Finite = me%Numite
               call me%ogwrit(1,"")
               write (str,'("OPTGRA: Converged: not ITERAT=",2I4,2D11.3)') &
                        me%Numite , me%Maxite , conerr , desnor
               call me%ogwrit(1,str)

               ! Final Pygmo output
               ! TODO: can this final fitness call be avoided (just for output)?
               me%Pygfla = 3 ! pygmo flag in COMMON: no covergence
               call evaluation_func_and_der()
               exit main
            elseif ( me%Numite>=me%Maxite .or. (me%Numite-itecor>=me%Optite-1 .and. itecor/=0) ) then
               Finopt = 2
               Finite = iteopt
               me%Varval = varcor
               me%Conval = concor
               call me%ogwrit(1,"")
               write (str,'("OPTGRA: Converged: mat ITERAT=",2I4,2D11.3)') &
                        me%Numite , me%Maxite , conerr , desnor
               call me%ogwrit(1,str)
               call me%ogpwri_end(-Valcon(me%Numcon+1),numvio,convio)

               ! Final Pygmo output
               ! TODO: can this final fitness call be avoided (just for output)?
               me%Pygfla = 2 ! pygmo flag in COMMON: constraints matched
               call evaluation_func_and_der()
               exit main
            endif
            ! ----------------------------------------------------------------------
            me%Numite = me%Numite + 1
            ! ----------------------------------------------------------------------
            call me%ogwrit(3,"")
            write (str,'("ITERAT=",I5)')me%Numite
            call me%ogwrit(2,str)
            ! ======================================================================
            ! GET VALUES AND GRADIENTS
            ! ======================================================================
            if ( me%Senopt<=0 ) then
               call evaluation_func_and_der()
            elseif ( me%Senopt==+1 .or. me%Senopt==+3 ) then
               me%Varval = me%Senvar
               call me%ogeval(me%Varval,me%Conval,0,conder_tmp)
               me%Conder(1:me%Numcon+1,:) = conder_tmp
            elseif ( me%Senopt==+2 .or. me%Senopt==+4 ) then
               me%Varval = me%Senvar
               do con = 1 , me%Numcon + 1
                  sca = me%Consca(con)
                  if ( me%Contyp(con)==-1 ) sca = -sca
                  me%Conval(con) = me%Sencon(con) - me%Sendel(con)/sca
               enddo
            endif
            if ( me%Senopt==-1 ) then
               me%Senvar = me%Varval
               me%Sencon = me%Conval
            endif
            ! ======================================================================
            if ( me%Varder==-1 .and. me%Senopt<=0 ) then
               me%Conred(1:me%Numcon+1,:) = me%Conder(1:me%Numcon+1,:)
               call me%ogeval(me%Varval,me%Conval,2,conder_tmp)
               me%Conder(1:me%Numcon+1,:) = conder_tmp
               write (str,'("GRADIENT CHECK")')
               call me%ogwrit(1,str)
               do var = 1 , me%Numvar
                  do con = 1 , me%Numcon + 1
                     fac = me%Varsca(var)/me%Consca(con)
                     fac = 1.0_wp
                     der = me%Conder(con,var)*fac
                     red = me%Conred(con,var)*fac
                     if ( abs(der)<1.0e-6_wp .and. abs(red)<1.0e-6_wp ) cycle
                     if ( abs(der-red)<1.0e-2_wp ) cycle
                     if ( der/=0.0_wp ) then
                        fac = red/der
                     else
                        fac = 0.0_wp
                     endif
                     if ( abs(fac-1.0_wp)<1.0e-2_wp ) cycle
                     write (str,'("VAR/CON/ANA/NUM/A2N=",2I4,3(1X,D13.6))') var , con , red , der , fac
                     call me%ogwrit(1,str)
                     nam = me%Varstr(var)
                     len = me%Varlen(var)
                     write (str,'("      VAR=",A,1X,D13.6)') nam(1:len) , me%Varval(var)*me%Varsca(var)
                     call me%ogwrit(1,str)
                     nam = me%Constr(con)
                     len = me%Conlen(con)
                     write (str,'("      CON=",A,1X,D13.6)') nam(1:len) , me%Conval(con)*me%Consca(con)
                     call me%ogwrit(1,str)
                  enddo
               enddo
!              CONDER(1:NUMCON+1,:) = CONRED(1:NUMCON+1,:)
!              GOTO 9999
            endif
            ! ======================================================================
            me%Sender = me%Conder
            do var = 1 , me%Numvar
               if ( me%Vartyp(var)/=1 ) cycle
!              WRITE (STR,*) "VAR=",VAR,VARVAL(VAR)*VARSCA(VAR)
!              CALL me%ogwrit (2,STR)
               me%Conder(1:me%Numcon+1,var) = 0.0_wp
            enddo
            ! ======================================================================
            if ( me%Numite==1 ) then
               me%Vargrd = me%Varval
            else
               me%Vargrd = me%Varref
            endif
            me%Varref = me%Varval
            me%Conref = me%Conval
            ! ======================================================================
            varacc = 0.0_wp
            ! ======================================================================
            cosold = cosnew
            cosnew = me%Conval(me%Numcon+1)
            call me%ogwrit(3,"")
            write (str,'("OPTGRA: VALCOS=",D15.8,1X,D15.8)') cosnew , cosnew - cosold
            call me%ogwrit(3,str)
            ! ======================================================================
            ! CORRECTION PART
            ! ----------------------------------------------------------------------
            call me%ogcorr(varacc,finish,conerr,norerr,error)
            if (error) then
               Finopt = -1
               return
            end if
            ! ----------------------------------------------------------------------
            if ( me%Tablev>=1 ) write (me%Tablun,'(I4,1X,"COR",1X,*(1X,D10.3))') &
               me%Numite , (me%Varval(var),var=1,me%Numvar) , (me%Conval(con),con=1,me%Numcon)
            ! ----------------------------------------------------------------------
            if ( me%Senopt/=0 ) then
               if ( finish/=0 ) exit inner
               Finopt = 0
               exit main
            endif
            ! ----------------------------------------------------------------------
            if ( finish==0 ) then
               me%Numact = 0
               old = meaerr
               itediv = itediv + 1
               num = min(itediv,me%Divite)
               meaerr = (meaerr*(num-1)+norerr)/num
!              WRITE (STR,*) "MEAERR=",MEAERR
!              CALL me%ogwrit (2,STR)
               if ( itediv<me%Divite .or. meaerr<=old ) cycle
               finish = -1
            endif
            ! ----------------------------------------------------------------------
            if ( finish==-1 ) then
               Finopt = 4
               Finite = me%Numite
               call me%ogwrit(1,"")
               write (str,'("OPTGRA: Converged: unf ITERAT=",2I4,2D11.3)') &
                        me%Numite , me%Maxite , conerr , desnor
               call me%ogwrit(1,str)

               ! Final Pygmo output
               ! TODO: can this final fitness call be avoided (just for output)?
               me%Pygfla = 4 ! pygmo flag in COMMON: infeasible
               call evaluation_func_and_der()
               exit main
            endif
            ! ----------------------------------------------------------------------
            itediv = 0
            iteopt = me%Numite
            if ( itecor==0 .or. concor(me%Numcon+1)<me%Conval(me%Numcon+1) ) then
               varcor = me%Varval
               concor = me%Conval
            endif
            if ( itecor==0 ) itecor = me%Numite
            ! ----------------------------------------------------------------------
            old = meamer
            itecnv = itecnv + 1
            num = min(itecnv,me%Cnvite)
            meamer = (meamer*(num-1)+concor(me%Numcon+1))/num
!           WRITE (STR,*) "MEAMER=",ITECNV,NUM,MEAMER,OLD,OLD/MEAMER
!           CALL me%ogwrit (-1,STR)
            if ( itecnv>=me%Cnvite .and. meamer<old ) then
               Finopt = 2
               Finite = iteopt
               me%Varval = varcor
               me%Conval = concor
               call me%ogwrit(1,"")
               write (str,'("OPTGRA: Converged: mat ITERAT=",2I4,2D11.3)') &
                        me%Numite , me%Maxite , conerr , desnor
               call me%ogwrit(1,str)

               ! Final Pygmo output
               ! TODO: can this final fitness call be avoided (just for output)?
               me%Pygfla = 2 ! pygmo flag in COMMON: matched
               call evaluation_func_and_der()
               exit main
            endif
            exit inner
         enddo inner

         ! ======================================================================
         ! OPTIMIZATION PART
         ! ----------------------------------------------------------------------
         if ( me%Senopt<+3 ) then
            varsav = me%Varmax
            me%Varmax = me%Varmax*10.0e-1_wp
            call me%ogopti(varacc,numequ,finish,desnor,error)
            if (error) then
               Finopt = -1
               return
            end if
            me%Varmax = varsav
         endif
         ! ----------------------------------------------------------------------
         if ( me%Senopt/=0 ) then
            call me%ogwrit(1,"")
            if ( finish==0 ) then
               Finopt = 0
               call me%ogwrit(1,"OPTGRA sensitivity converged: not")
            else
               Finopt = 1
               call me%ogwrit(1,"OPTGRA sensitivity converged: yes")
            endif
            exit main
         endif
         ! ----------------------------------------------------------------------
         if ( finish==0 ) cycle main

         ! ======================================================================
         ! NOT CONVERGED
         ! ----------------------------------------------------------------------
         if ( varacc/=0.0_wp ) cycle main

         ! ======================================================================
         ! CONVERGED
         ! ----------------------------------------------------------------------
         Finopt = 1
         Finite = me%Numite
         call me%ogwrit(1,"")
         write (str,'("OPTGRA: Converged: yes ITERAT=",2I4,2D11.3)') &
                  me%Numite , me%Maxite , conerr , desnor
         call me%ogwrit(1,str)
         call me%ogwrit(3,"")

         ! Final Pygmo output
         ! TODO: can this final fitness call be avoided (just for output)?
         me%Pygfla = 1 ! covergence
         call evaluation_func_and_der()
         exit main

      enddo main

!     WRITE (STR,*) "DIF=",NORM2(VARVAL-VARREF)
!     CALL me%ogwrit (1,STR)
      ! ======================================================================
      ! DESCALE VARIABLES
      ! ----------------------------------------------------------------------
!     CALL OGWMAT (3)
!     CALL me%ogeval (VARVAL, VALCON, 0, CONDER)
      do var = 1 , me%Numvar
         Valvar(var) = me%Varval(var)*me%Varsca(var)
      enddo
!     IF (SENOPT /= 0) THEN
!         CALL me%ogeval (VARVAL, VALCON, 0, CONDER)
!     ENDIF
      ! ======================================================================
      ! DESCALE VALUES
      ! ----------------------------------------------------------------------
      do con = 1 , me%Numcon + 1
         typ = me%Contyp(con)
         sca = me%Consca(con)
         if ( typ==-1 ) sca = -sca
         Valcon(con) = me%Conval(con)*sca
      enddo
      ! ======================================================================
      call me%ogwrit(3,"")
      call me%ogwrit(3,"STATUS OF CONSTRAINTS:")
      call me%ogwrit(3,"")
      call me%ogwrit(3," ACT  PAS  NON COST___VAL CONSTRAINT")
      do con = 1 , me%Numcon
         if ( me%Contyp(con)==-2 ) cycle
         nam = me%Constr(con)
         len = me%Conlen(con)
         val = me%Conval(con)
         if ( me%Conact(con)>0 ) then
            write (str,'( I4,5X,6X,D10.3,1X,A)') con , val , nam(1:len)
            call me%ogwrit(3,str)
         elseif ( me%Conact(con)==0 ) then
            write (str,'( 5X,I4,6X,D10.3,1X,A)') con , val , nam(1:len)
            call me%ogwrit(3,str)
         elseif ( me%Conact(con)<0 ) then
            write (str,'(10X,I4,1X,D10.3,1X,A)') con , val , nam(1:len)
            call me%ogwrit(3,str)
         endif
      enddo

      deallocate (varsum)
      deallocate (varcor)
      deallocate (concor)
      deallocate (conder_tmp)
      ! ----------------------------------------------------------------------
      call me%ogwrit(2,"")
      call me%ogwrit(2,"OPTGRA END")
      call me%ogwrit(2,"")

      contains

         subroutine evaluation_func_and_der()
            !! evaluate the function and derivatives
            call me%ogeval(me%Varval,me%Conval,me%Varder,conder_tmp)
            me%Conder(1:me%Numcon+1,:) = conder_tmp
         end subroutine evaluation_func_and_der

   end subroutine ogexec

   subroutine oggsst(me,Varsen,Quasen,Consen,Actsen,Dersen,Actsav,Consav,Redsav,Dersav,Actnum)

      !! NEAR-LINEAR OPTIMIZATION TOOL SENSITIVITY ANALYSIS
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

      do var = 1 , me%Numvar
         Varsen(var) = me%Senvar(var)
      enddo

      do con = 1 , me%Numcon + 1
         Quasen(con) = me%Senqua(con)
         Consen(con) = me%Sencon(con)
         Actsen(con) = me%Senact(con)
         do var = 1 , me%Numvar
            Dersen(con,var) = me%Sender(con,var)
         enddo
      enddo

      ! Temporary status saved of which constraints are active
      do con = 1 , me%Numcon + 1
         Actsav(con) = me%Actcon(con)
      enddo

      do con = 1 , me%Numcon + 4
         Consav(con) = me%Conact(con)
      enddo

      do con = 1 , me%Numcon + 3
         do var = 1 , me%Numvar
            Redsav(con,var) = me%Conred(con,var)
            Dersav(con,var) = me%Conder(con,var)
         enddo
      enddo

   end subroutine oggsst

   subroutine ogincl(me,Inc)

      !! ADDS CONSTRAINT TO ACTIVE SET AND REDUCES DERIVATIVES
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      integer(ip),intent(in) :: Inc !! CONSTRAINT TO BE INCLUDED

      real(wp) :: val , fac , gam , sav , max
      integer(ip) :: row , col , ind , lst
      character(len=str_len) :: str

      ! GENERAL
      me%Numact = me%Numact + 1

      ! PERMUTATION TO GET ZERO DERIVATIVES AT END FOR NEW ACTIVE CONSTRAINT
      lst = me%Numvar
      do col = me%Numvar , me%Numact , -1
         if ( me%Conred(Inc,col)==0.0_wp ) then
            if ( col/=lst ) then
               do row = 1 , me%Numcon + 3
                  if ( me%Conact(row)<=0 ) then
                     sav = me%Conred(row,col)
                     me%Conred(row,col) = me%Conred(row,lst)
                     me%Conred(row,lst) = sav
                  endif
               enddo
            endif
            lst = lst - 1
         endif
      enddo

      ! PERMUTATION TO GET MAXIMUM PIVOT
      ind = me%Numact
      max = abs(me%Conred(Inc,ind))
      do col = me%Numact + 1 , lst
         val = abs(me%Conred(Inc,col))
         if ( val>max ) then
            ind = col
            max = val
         endif
      enddo

      if ( ind/=me%Numact ) then
         do row = 1 , me%Numcon + 3
            if ( me%Conact(row)<=0 ) then
               sav = me%Conred(row,ind)
               me%Conred(row,ind) = me%Conred(row,me%Numact)
               me%Conred(row,me%Numact) = sav
            endif
         enddo
      endif

      ! UPDATE LIST OF ACTIVE CONSTRAINTS
      me%Actcon(me%Numact) = Inc
      me%Conact(Inc) = me%Numact

      ! REDUCE FOR NEW ACTIVE CONSTRAINT
      if ( abs(me%Conred(Inc,me%Numact))<1.0e-12_wp ) then
         write (str,*) "OGINCL-WARNING: CONSTRAINT SINGULAR"
         call me%ogwrit(2,str)
         write (str,*) "INC=" , Inc
         call me%ogwrit(2,str)
         write (str,*) "PIV=" , me%Conred(Inc,me%Numact)
         call me%ogwrit(2,str)
         me%Numact = me%Numact - 1
         me%Conact(Inc) = 0
         return
      endif

      val = sqrt(sum(me%Conred(Inc,me%Numact:lst)**2))
      if ( me%Conred(Inc,me%Numact)>0.0_wp ) val = -val

      me%Conred(Inc,me%Numact) = me%Conred(Inc,me%Numact) - val

      sav = me%Conred(Inc,me%Numact)
      fac = 1.0_wp/sav
      me%Conred(Inc,me%Numact:lst) = me%Conred(Inc,me%Numact:lst)*fac

      fac = sav/val
      do row = 1 , me%Numcon + 3
         if ( me%Conact(row)<=0 ) then
            gam = dot_product(me%Conred(row,me%Numact:lst),me%Conred(Inc,me%Numact:lst))
            if ( gam/=0.0_wp ) then
               gam = gam*fac
               me%Conred(row,me%Numact:lst) = me%Conred(row,me%Numact:lst) + &
                                              me%Conred(Inc,me%Numact:lst)*gam
            endif
         endif
      enddo

      me%Conred(Inc,me%Numact) = val
      me%Conred(Inc,me%Numact+1:lst) = 0.0_wp

   end subroutine ogincl

   subroutine ogleft(me,Actinp,Actout)

      !! LEFT-MULTIPLIES VECTOR LOWER TRIANGULAR MATRIX OBTAINED BY REDUCTION
      !! AND SUBSEQUENT INVERSION OF DERIVATIVES OF ACTIVE CONSTRAINTS
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Actinp(me%Numcon) !! VECTOR INITAL
      real(wp),intent(out) :: Actout(me%Numcon) !! VECTOR FINAL (MAY BE SAME AS ACTINP)

      integer(ip) :: row , col , act
      real(wp) :: val

      do act = 1 , me%Numact
         row = me%Actcon(act)
         val = Actinp(act)
         do col = 1 , act - 1
            val = val - me%Conred(row,col)*Actout(col)
         enddo
         Actout(act) = val/me%Conred(row,act)
      enddo

   end subroutine ogleft

   subroutine ogopti(me,Varacc,Numequ,Finish,Desnor,error)

      !! OPTIMIZATION PART
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(inout) :: Varacc !! ITERATION SCALED DISTANCE ACCUMULATED
      integer(ip),intent(in) :: Numequ !! NUMBER OF EQUALITY CONSTRAINTS [not used?]
      integer(ip),intent(out) :: Finish !! 0=LIMIT 1=OPTIM
      logical,intent(out) :: error !! if there was a fatal error

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
      character(len=str_len) :: str
      character(len=name_len) :: nam

      real(wp) , dimension(:) , allocatable :: cosact
      real(wp) , dimension(:) , allocatable :: varvec
      real(wp) , dimension(:) , allocatable :: varwrk
      real(wp) , dimension(:) , allocatable :: corvec
      real(wp) , dimension(:) , allocatable :: desder
      real(wp) , dimension(:) , allocatable :: desprv
      real(wp) , dimension(:) , allocatable :: varprv
      real(wp) , dimension(:) , allocatable :: convec
      real(wp) , dimension(:) , allocatable :: conqua
      integer(ip) , dimension(:) , allocatable :: concor

      error = .false.
      ! initialize:  - JW : these could also be in the class.
      !                     that would save the allocate/reallocate
      !                     step every time this routine is called.
      allocate (cosact(me%Numcon))
      allocate (varvec(me%Numvar))
      allocate (varwrk(me%Numvar))
      allocate (corvec(me%Numcon))
      allocate (desder(me%Numvar))
      allocate (desprv(me%Numvar))
      allocate (varprv(me%Numvar))
      allocate (convec(me%Numcon+1))
      allocate (conqua(me%Numcon+1))
      allocate (concor(me%Numcon+1))

      ! ======================================================================
      ! OPTIMIZATION PART
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
      call me%ogwrit(3,"")
      call me%ogwrit(3,"OPTIMIZATION PART")
      ! ----------------------------------------------------------------------
!     WRITE (STR,'("NUMACT = ",I4)') NUMACT
!     CALL me%ogwrit (3,STR)
      numcor = me%Numact
      concor = me%Actcon
      me%Numact = 0
      me%Conact(1:des) = 0

      main: do
!        DO COR = 1, NUMCOR
!            CON = CONCOR(COR)
!            NAM = CONSTR(CON)
!            LEN = CONLEN(CON)
!            WRITE (STR,'("ACT = ",I4,5X,1X,A)') CON, NAM(1:LEN)
!            CALL me%ogwrit (3,STR)
!            CALL me%OGINCL (CON)
!        ENDDO
         ! ======================================================================
         ! VECTOR OF STEEPEST ASCENT
         ! ----------------------------------------------------------------------
         call me%ogwrit(3,"")
         call me%ogwrit(3,"VECTOR OF STEEPEST ASCENT")
         call me%ogwrit(3,"")
         call me%ogwrit(3,"REMOVE AND INCLUDE CONSTRAINTS:")
         call me%ogwrit(3,"")
         call me%ogwrit(3," REM  INC CONSTRAINT")
         ! ----------------------------------------------------------------------
         ! REMOVE PASSIVE INEQUALITY CONSTRAINTS
         ! ----------------------------------------------------------------------
         do act = me%Numact , 1 , -1
            con = me%Actcon(act)
            if ( me%Conval(con)<=1.0_wp ) cycle
            nam = me%Constr(con)
            len = me%Conlen(con)
            write (str,'(I4,5X,1X,A)') con , nam(1:len)
            call me%ogwrit(2,str)
            call me%ogexcl(act,error)
            if (error) return
            me%Conact(con) = -1
         enddo
         ! ----------------------------------------------------------------------
         ! INCLUDE VIOLATED INEQUALITY CONSTRAINTS AND SELECT PASSIVE ONES
         ! SELECT PASSIVE INEQUALITY CONSTRAINTS
         ! ----------------------------------------------------------------------
         do con = 1 , me%Numcon
            if ( me%Contyp(con)==-2 ) cycle
            if ( me%Conact(con)>0 ) then
            elseif ( me%Contyp(con)==0 ) then
               me%Conact(con) = 0
            elseif ( me%Conval(con)<-1.0_wp ) then
               me%Conact(con) = 0
            elseif ( me%Conval(con)<=+1.0_wp ) then
               me%Conact(con) = 0
            else
               me%Conact(con) = -1
            endif
         enddo
         ! ======================================================================
         nnn = 1
         inner: do
            nnn = nnn + 1
            if ( nnn>999 ) then
               Finish = 0
               write (str,*) "NNN=" , nnn
               call me%ogwrit(2,str)
               exit main
            endif
            ! ======================================================================
            ! DERIVATIVES OF MERIT W.R.T. ACTIVE CONSTRAINTS
            ! ----------------------------------------------------------------------
            call me%ogrigt(-me%Conred(cos,1:me%Numact),cosact)
            Desnor = sqrt(sum(me%Conred(cos,me%Numact+1:me%Numvar)**2))
            ! ----------------------------------------------------------------------
            ! CONSTRAINT REMOVAL
            ! ----------------------------------------------------------------------
            ind = 0
            exc = -1.0e-12_wp
            max = exc
            do act = 1 , me%Numact
               con = me%Actcon(act)
               if ( me%Contyp(con)==0 ) cycle
               val = cosact(act)
               fac = dot_product(me%Conred(con,1:me%Numvar),me%Conred(con,1:me%Numvar))
               fac = sqrt(fac)
               val = val*fac
               if ( val>=exc ) cycle
               if ( val>max ) cycle
               max = val
               ind = act
            enddo
            ! ----------------------------------------------------------------------
            if ( ind/=0 ) then
               con = me%Actcon(ind)
               nam = me%Constr(con)
               len = me%Conlen(con)
               write (str,'(I4,5X,3(1X,D10.3),1X,A)') con , Desnor , max , me%Varmax , nam(1:len)
               call me%ogwrit(3,str)
               call me%ogexcl(ind,error)
               if (error) return
               cycle
            endif
            ! ----------------------------------------------------------------------
            ! CONSTRAINT INCLUSION
            ! ----------------------------------------------------------------------
            if ( Desnor/=0.0_wp ) then
               ! ----------------------------------------------------------------------
               inc = 0
               eps = 1.0e-03_wp
               max = -1.0e10_wp     ! JW : what is going on here ?
               max = 0.0_wp
               ! ----------------------------------------------------------------------
               do con = 1 , me%Numcon
                  if ( me%Contyp(con)==-2 ) cycle
                  if ( me%Conact(con)/=0 ) cycle
                  del = dot_product(me%Conred(con,me%Numact+1:me%Numvar),&
                                    me%Conred(cos,me%Numact+1:me%Numvar))/Desnor
                  val = abs(del)*me%Varmax
                  if ( val<eps ) cycle
                  fac = dot_product(me%Conred(con,1:me%Numvar),me%Conred(con,1:me%Numvar))
                  fac = sqrt(fac)
                  del = del/fac
                  if ( del<0.0_wp .and. del<max ) then
                     max = del
                     inc = con
                  endif
                  if ( me%Contyp(con)/=0 ) cycle
                  del = -del
                  if ( del<0.0_wp .and. del<max ) then
                     max = del
                     inc = con
                  endif
               enddo
               ! ----------------------------------------------------------------------
               if ( inc/=0 ) then
                  con = inc
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  write (str,'(5X,I4,3(1X,D10.3),1X,A)') con , Desnor , max , me%Varmax , nam(1:len)
                  call me%ogwrit(3,str)
                  call me%ogincl(inc)
                  cycle
               endif
            endif
            ! ----------------------------------------------------------------------
            ! MATCHED INEQUALITY CONSTRAINTS + STEEPEST ASCENT VECTOR NORM
            ! ----------------------------------------------------------------------
            call me%ogwrit(3,"")
            call me%ogwrit(3,"STATUS OF MATCHED INEQUALITY CONSTRAINTS:")
            call me%ogwrit(3,"")
            call me%ogwrit(3," ACT  PAS CONSTRAINT")
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               nam = me%Constr(con)
               len = me%Conlen(con)
               if ( me%Contyp(con)==0 ) then
               elseif ( me%Conact(con)>0 ) then
                  write (str,'(I4,5X,1X,A)') con , nam(1:len)
                  call me%ogwrit(3,str)
               elseif ( me%Conact(con)==0 ) then
                  write (str,'(5X,I4,1X,A)') con , nam(1:len)
                  call me%ogwrit(3,str)
               endif
            enddo
            call me%ogwrit(3,"")
            write (str,'("STEEPEST ASCENT NORM: ",D13.6)') Desnor
            call me%ogwrit(3,str)
            write (str,'("MAXIMUM DISTANCE....: ",D13.6)') me%Varmax
            call me%ogwrit(3,str)
            ! ======================================================================
            Finish = 0
            ! ======================================================================
            ! IF CONVERGENCE
            ! ----------------------------------------------------------------------
            cosimp = Desnor*me%Varmax
            if ( abs(cosimp)<=1.0_wp ) then
               foldis = 0.0_wp
               Finish = 1
               write (str,'("FINAL...............:",1X,D13.6,'//'11X,1(1X,D10.3),1X,D16.9)') &
                  foldis , cosimp , me%Conval(cos) + cosimp
               call me%ogwrit(2,str)
               exit main
            endif
            ! ======================================================================
            ! IF CONSTRAINT IS HIT
            ! ----------------------------------------------------------------------
            do var = 1 , me%Numvar
               val = me%Conder(cos,var)
               do act = 1 , me%Numact
                  con = me%Actcon(act)
                  val = val + me%Conder(con,var)*cosact(act)
               enddo
               me%Vardes(var) = val
            enddo
            ! ----------------------------------------------------------------------
            ind = 0
            dis = 1.0e10_wp
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               if ( me%Conact(con)/=-1 ) cycle
               val = dot_product(me%Conred(con,me%Numact+1:me%Numvar),me%Conred(cos,me%Numact+1:me%Numvar))
               if ( val==0.0_wp ) cycle
               val = -me%Conval(con)/val*Desnor
               if ( val<=0.0_wp ) cycle
               if ( val>=dis ) cycle
               dis = val
               ind = con
            enddo
            ! ----------------------------------------------------------------------
            if ( ind/=0 ) then
               val = sqrt(sum((me%Varval-me%Varref+me%Vardes*dis/Desnor)**2))
               if ( val>me%Varmax ) ind = 0
            endif
            ! ----------------------------------------------------------------------
            if ( ind/=0 ) then
               if ( me%Confix(ind)<=0 ) then
                  if ( val>me%Varmax*1.0e-1_wp ) ind = 0
               endif
            endif
            ! ----------------------------------------------------------------------
            if ( ind/=0 ) then
               imp = dis*Desnor
               con = ind
               tht = 1.0_wp
               bet = 0.0_wp
               nam = me%Constr(con)
               len = me%Conlen(con)
               write (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,11X,1(1X,D10.3),1X,D16.9,22X,1X,I4,1X,A)') &
                  dis , imp , me%Conval(cos) + cosimp , con , nam(1:len)
               call me%ogwrit(2,str)
               Varacc = Varacc + dis
               me%Varval = me%Varval + dis*me%Vardes/Desnor
               do con = 1 , me%Numcon + 1
                  val = dot_product(me%Conred(con,me%Numact+1:me%Numvar),me%Conred(cos,me%Numact+1:me%Numvar))
                  me%Conval(con) = me%Conval(con) + val*dis/Desnor
               enddo
               cycle main
            endif
            exit inner
         enddo inner

         inner2: do
            ! ======================================================================
            ! ----------------------------------------------------------------------
            call me%ogrigt(-me%Conred(cos,1:me%Numact),cosact)
            do var = 1 , me%Numvar
               val = me%Conder(cos,var)
               do act = 1 ,me%Numact
                  con = me%Actcon(act)
                  val = val + me%Conder(con,var)*cosact(act)
               enddo
               desprv(var) = val
            enddo
            Desnor = sqrt(sum(desprv**2))
            write (str,'("DESNOR=",D13.6)') Desnor
!              CALL me%ogwrit (2,STR)
            ! ----------------------------------------------------------------------
            call me%ogrigt(-me%Conred(prv,1:me%Numact),cosact)
            do var = 1 , me%Numvar
               val = me%Conder(prv,var)
               do act = 1 , me%Numact
                  con = me%Actcon(act)
                  val = val + me%Conder(con,var)*cosact(act)
               enddo
               varprv(var) = val
            enddo
            norprv = sqrt(sum(varprv**2))
            write (str,'("NORPRV=",D13.6)') norprv
!              CALL me%ogwrit (2,STR)
            ! ----------------------------------------------------------------------
            call me%ogrigt(-me%Conred(des,1:me%Numact),cosact)
            do var = 1 , me%Numvar
               val = desder(var)
               do act = 1 , me%Numact
                  con = me%Actcon(act)
                  val = val + me%Conder(con,var)*cosact(act)
               enddo
               me%Vardir(var) = val
            enddo
            nor = sqrt(sum(me%Vardir**2))
            write (str,'("NOR=",D13.6)') nor
!              CALL me%ogwrit (2,STR)
            met = me%Optmet
            tht = 1.0_wp
            bet = 0.0_wp
            select case (met)
               case ( 0 ) ! STEEPEST DESCENT METHOD
               case ( 1 ) ! MODIFIED SPECTRAL CONJUGATE GRADIENT METHOD
               varvec = desprv - varprv
               if ( norprv**2>1.0e-12_wp ) then
                  tht = -dot_product(me%Vardir,varvec)/norprv**2
                  bet = Desnor**2/norprv**2
               endif
               case ( 2 ) ! SPECTRAL CONJUGATE GRADIENT METHOD
               varvec = desprv - varprv
               varwrk = me%Varref - me%Vargrd
               val = dot_product(varwrk,varvec)
               fac = dot_product(me%Vardir,varvec)
               if ( abs(val)>1.0e-12_wp .and. abs(fac)>1.0e-12_wp ) then
                  tht = -dot_product(varwrk,varwrk)/val
                  varwrk = -varvec*tht - varwrk
                  bet = dot_product(varwrk,desprv)/fac
               endif
               case ( 3 ) ! CONJUGATE GRADIENT METHOD
               if ( norprv/=0.0_wp ) then
                  tht = 1.0_wp
                  bet = Desnor**2/norprv**2
               endif
            end select
!           WRITE (STR,'("THT=",D13.6)') THT
!           CALL me%ogwrit (3,STR)
!           WRITE (STR,'("BET=",D13.6)') BET
!           CALL me%ogwrit (3,STR)
            ! ----------------------------------------------------------------------
            eps = 1.0e-03_wp
!           WRITE (STR,*) "THT/BET=",THT,BET
!           CALL me%ogwrit (2,STR)
            me%Vardes = tht*desprv + bet*me%Vardir
            Desnor = sqrt(sum(me%Vardes**2))
            nor = Desnor
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               if ( me%Conact(con)/=0 ) cycle
               del = dot_product(me%Conder(con,1:me%Numvar),me%Vardes(1:me%Numvar))/nor
               val = abs(del)*me%Varmax
               if ( val<eps ) cycle
               fac = dot_product(me%Conder(con,1:me%Numvar),me%Conder(con,1:me%Numvar))
               del = del/sqrt(fac)
               nam = me%Constr(con)
               len = me%Conlen(con)
               typ = me%Contyp(con)
               if ( del<0.0_wp ) then
!                 CALL OGINCL (CON)
                  act = me%Conact(con)
!                 WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)') &
!                    CON,ACT,+NOR,VAL,DEL,NAM(1:LEN)
!                 CALL me%ogwrit (2,STR)
                  if ( act/=0 ) cycle inner2
                  bet = 0.0_wp
                  tht = 1.0_wp
               endif
               if ( me%Contyp(con)/=0 ) cycle
               del = -del
               if ( del<0.0_wp ) then
!                 CALL OGINCL (CON)
                  act = me%Conact(con)
!                 WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)') &
!                    CON,ACT,-NOR,VAL,DEL,NAM(1:LEN)
!                 CALL me%ogwrit (2,STR)
                  if ( act/=0 ) cycle inner2
                  bet = 0.0_wp
                  tht = 1.0_wp
                  del = dot_product(me%Conder(con,1:me%Numvar),me%Conder(cos,1:me%Numvar))
                  val = abs(del)*me%Varmax/Desnor
!                 WRITE (STR,'(5X,2I4,3(1X,D10.3),1X,A)') &
!                    CON,TYP,-DESNOR,VAL,DEL,NAM(1:LEN)
!                 CALL me%ogwrit (2,STR)
               endif
            enddo
            cosco1 = dot_product(me%Vardes,me%Conder(cos,1:me%Numvar))/Desnor
            if ( cosco1<0.0_wp ) then
               write (str,*) "COSCO1=" , cosco1
               call me%ogwrit(2,str)
               bet = 0.0_wp
               tht = 1.0_wp
            endif
            me%Vardes = tht*desprv + bet*me%Vardir
            Desnor = sqrt(sum(me%Vardes**2))
            exit inner2
         enddo inner2

         inner3: do
!           WRITE (STR,*) "THT/BET=",THT,BET
!           CALL me%ogwrit (2,STR)
            ! ======================================================================
            ! SECOND ORDER DERIVATIVES BY NUMERICAL DIFFERENCING
            ! ----------------------------------------------------------------------
            call me%ogwrit(3,"")
            write (str,'("SECOND ORDER CORRECTION")')
            call me%ogwrit(3,str)
            me%Vardes = tht*desprv + bet*me%Vardir
            Desnor = sqrt(sum(me%Vardes**2))
            ! ----------------------------------------------------------------------
            ! MAXIMUM TRAVEL DISTANCE
            ! ----------------------------------------------------------------------
            dis = me%Varmax
            varvec = me%Varval - me%Varref
            co0 = dot_product(varvec,varvec) - dis**2
            if ( co0>=0.0_wp ) then
               dis = 0.0_wp
            else
               co1 = dot_product(me%Vardes,varvec)
               co2 = Desnor**2
               det = co1**2 - co0*co2
               dis = (sqrt(det)-co1)/co2
            endif
            dis = dis*Desnor
            maxdis = dis
            ! ======================================================================
            ! COMPUTE SECOND ORDER EFFECTS
            ! ----------------------------------------------------------------------
            nnn = 0
            if ( me%Senopt>=+1 ) then
               do con = 1 , me%Numcon
                  if ( me%Contyp(con)==-2 ) cycle
                  act = me%Conact(con)
                  ind = me%Senact(con)
                  if ( act==-1 .and. ind==-1 ) cycle
                  if ( act==0 .and. ind==0 ) cycle
                  if ( act>0 .and. ind>0 ) cycle
                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  write (str,'(I4,1X,I4,1X,I4,1X,A)') con , act , ind , nam(1:len)
                  call me%ogwrit(2,str)
                  nnn = 1
               enddo
            endif
            if ( me%Senopt<=0 .or. nnn==1 ) then
               fac = me%Varstp/Desnor
               varvec = me%Varref + me%Vardes*me%Varstp/Desnor
               call me%ogeval(varvec,convec,0,me%Conder)
               conqua = matmul(me%Conder(1:cos,1:me%Numvar),me%Vardes(1:me%Numvar))
               conqua = 2.0_wp*(convec-me%Conref-conqua*fac)/fac**2
            endif
            if ( me%Senopt==-1 ) then
               me%Senqua = conqua
            elseif ( me%Senopt>=+1 .and. nnn==0 ) then
               conqua = me%Senqua
            endif
            ! ======================================================================
            ! COMPUTE CORRECTION VECTOR
            ! ----------------------------------------------------------------------
            do act = 1 , me%Numact
               con = me%Actcon(act)
               corvec(act) = conqua(con)
            enddo
            call me%ogleft(corvec,corvec)
            ! ----------------------------------------------------------------------
            cornor = sqrt(sum(corvec(1:me%Numact)**2))*0.5_wp/Desnor/Desnor
            call me%ogwrit(3,"")
            write (str,'("STEEPEST ASCENT  NORM: ",D13.6)') Desnor
            call me%ogwrit(3,str)
            write (str,'("ACCUMULATED  DISTANCE: ",D13.6)') Varacc
            call me%ogwrit(3,str)
            ! ======================================================================
            ! GET EXTREMUM
            ! ----------------------------------------------------------------------
            cosco1 = dot_product(me%Vardes,me%Conder(cos,1:me%Numvar))/Desnor
            cosco2 = conqua(cos) - dot_product(me%Conred(cos,1:me%Numact),corvec(1:me%Numact))
            cosco2 = cosco2*0.5_wp/Desnor/Desnor
            write (str,*) "COSCO2/COSCO1=" , cosco2 , cosco1
            call me%ogwrit(3,str)
            if ( cosco1<0.0_wp ) call me%ogwrit(2,str)
            ! ----------------------------------------------------------------------
            foldis = 0.0_wp
            quacor = cornor*foldis*foldis
            cosimp = foldis*(cosco1+foldis*cosco2)
            call me%ogwrit(3,"")
            write (str,'(    "STEEPEST ASCENT FOLLOW",'//'  5X,"DISTANCE",'//&
                       '  1X,"CORRECTION",'//'  2X,"MERIT_DEL",'//'  6X,"MERIT_VALUE")')
            call me%ogwrit(3,str)
            write (str,'("INITIAL.............:",1X,D13.6,'//'  2(1X,D10.3),1X,D16.9)') &
                     foldis , quacor , cosimp ,me%Conval(cos) + cosimp
            call me%ogwrit(3,str)
            ! ======================================================================
            if ( cosco2<0.0_wp ) then
               foldis = -0.5_wp*cosco1/cosco2
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               write (str,'("MERIT MAXIMUM.......:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                  foldis , quacor , cosimp , me%Conval(cos) + cosimp
               call me%ogwrit(3,str)
               staflg = 1
            elseif ( cosco2>0.0_wp ) then
               foldis = me%Varmax
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               write (str,'("MERIT CONVEX........:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                  foldis , quacor , cosimp , me%Conval(cos) + cosimp
               call me%ogwrit(3,str)
               staflg = 2
            else
               foldis = me%Varmax
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               write (str,'("MERIT LINEAR........:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                  foldis , quacor , cosimp , me%Conval(cos) + cosimp
               call me%ogwrit(3,str)
               staflg = 2
            endif
            ! ======================================================================
            ! IF MAXIMUM DISTANCE IS HIT
            ! ----------------------------------------------------------------------
            if ( foldis>me%Varmax ) then
               foldis = me%Varmax
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               write (str,'("MAXIMUM DISTANCE....:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                  foldis , quacor , cosimp , me%Conval(cos) + cosimp
               call me%ogwrit(3,str)
               staflg = 2
            endif
            ! ======================================================================
            ! IF CONVERGENCE
            ! ----------------------------------------------------------------------
            if ( abs(cosimp)<=1.0_wp ) then
               write (str,'("FINAL...............:",1X,D13.6,'//'  2(1X,D10.3),1X,D16.9,2D11.3)') &
                  foldis , quacor , cosimp , me%Conval(cos) + cosimp , tht , bet
               call me%ogwrit(2,str)
               if ( tht/=1.0_wp .or. bet/=0.0_wp ) then
                  tht = 1.0_wp
                  bet = 0.0_wp
                  cycle
               endif
               foldis = 0.0_wp
               Finish = 1
               exit main
            endif
            ! ======================================================================
            ! IF REMAINING DISTANCE IS HIT
            ! ----------------------------------------------------------------------
            if ( foldis>maxdis ) then
               foldis = maxdis
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               write (str,'("REMAINING DISTANCE..:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
                  foldis , quacor , cosimp , me%Conval(cos) + cosimp
               call me%ogwrit(3,str)
               staflg = 2
            endif
            ! ======================================================================
            ! IF CONSTRAINT IS HIT
            ! ----------------------------------------------------------------------
            ind = 0
            do con = 1 , me%Numcon
               if ( me%Contyp(con)==-2 ) cycle
               if ( me%Conact(con)/=-1 ) cycle
               co2 = conqua(con) - dot_product(me%Conred(con,1:me%Numact),corvec(1:me%Numact))
               co1 = dot_product(me%Conder(con,1:me%Numvar),me%Vardes(1:me%Numvar))
               co0 = me%Conval(con)*2.0_wp
               if ( co2/=0.0_wp ) then
                  det = co1**2 - co2*co0
                  if ( det<0.0_wp ) cycle
                  det = sqrt(det)
                  val = 1.0e10_wp
                  fac = (-co1+det)/co2
                  if ( fac>0.0_wp .and. fac<val ) val = fac
                  fac = (-co1-det)/co2
                  if ( fac>0.0_wp .and. fac<val ) val = fac
               elseif ( co1/=0.0_wp ) then
                  val = -co0/co1*0.5_wp
               else
                  cycle
               endif
               val = val*Desnor
               if ( val>0.0_wp .and. val<foldis ) then
                  foldis = val
                  ind = con
               endif
            enddo
            ! ----------------------------------------------------------------------
            if ( ind/=0 ) then
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               con = ind
               nam = me%Constr(con)
               len = me%Conlen(con)
               write (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,2(1X,D10.3),1X,D16.9,1X,I4,1X,A)') &
                  foldis , quacor , cosimp , me%Conval(cos) + cosimp , con , nam(1:len)
               call me%ogwrit(3,str)
               staflg = 3
            endif
            ! ======================================================================
            ! UPDATE
            ! ----------------------------------------------------------------------
            refdis = foldis
            exit inner3
         enddo inner3

         inner4: do
            write (str,'("FINAL...............:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9)') &
               foldis , quacor , cosimp , me%Conval(cos) + cosimp
            call me%ogwrit(3,str)
            ! ----------------------------------------------------------------------
            fac = foldis/Desnor
            ! ----------------------------------------------------------------------
            ! VARIABLE DELTA
            ! ----------------------------------------------------------------------
            call me%ogrigt(corvec,cosact)
            dis = 0.0_wp
            do var = 1 , me%Numvar
               val = 0.0_wp
               do act = 1 , me%Numact
                  ind = me%Actcon(act)
                  val = val - cosact(act)*me%Conder(ind,var)
               enddo
               varvec(var) = fac*(me%Vardes(var)+(val*fac*0.5_wp))
               dis = dis + varvec(var)*varvec(var)
            enddo
            dis = sqrt(dis)
            ! ----------------------------------------------------------------------
            write (str,*) "REFDIS=" , refdis
            call me%ogwrit(3,str)
            write (str,*) "FOLDIS=" , foldis
            call me%ogwrit(3,str)
            write (str,*) "DIS=" , dis
            call me%ogwrit(3,str)
            if ( dis>refdis*1.2_wp .and. me%Senopt>0 ) then
               faccnt = faccnt + 1
               if ( faccnt>=10 ) exit inner4
               foldis = foldis*0.5_wp
               quacor = cornor*foldis*foldis
               cosimp = foldis*(cosco1+foldis*cosco2)
               cycle
            endif
            ! ----------------------------------------------------------------------
            ! UPDATE VARIABLES
            ! ----------------------------------------------------------------------
            Varacc = Varacc + foldis
            me%Varval = me%Varval + varvec
            ccc = sqrt(sum((me%Varval-me%Varref)**2)) - me%Varmax**2
            if ( ccc>=0.0_wp ) then
               write (str,*) "CCC > 0" , ccc
               call me%ogwrit(3,str)
               staflg = 2
            endif

            select case (staflg)

               case ( 1 )  ! MAXIMUM REACHED: NEXT ITERATION

                  write (str,'("MERIT MAXIMUM.......:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9,2D11.3)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp , tht , bet
                  call me%ogwrit(2,str)
                  if ( me%Senopt>0 ) Finish = 1
                  exit inner4

               case ( 2 )   ! MAXIMUM TRAVEL DISTANCE REACHED: NEXT ITERATION

                  write (str,'("REMAINING DISTANCE..:",1X,D13.6,'//'2(1X,D10.3),1X,D16.9,2D11.3)') &
                     foldis , quacor , cosimp , me%Conval(cos) + cosimp , tht , bet
                  call me%ogwrit(2,str)
                  exit inner4

               case ( 3 )   ! CONSTRAINT HIT: UPDATE CONSTRAINT + CORRECT

                  nam = me%Constr(con)
                  len = me%Conlen(con)
                  write (str,'( "CONSTRAINT REACHED..:",'//'1X,D13.6,2(1X,D10.3),1X,D16.9,2D11.3,1X,I4,1X,A)') &
                           foldis , quacor , cosimp , me%Conval(cos) + cosimp , tht , bet , con , nam(1:len)
                  call me%ogwrit(2,str)
                  convec = conqua - matmul(me%Conred(1:cos,1:me%Numact),corvec(1:me%Numact))
                  convec = convec*fac*0.5_wp
                  convec = convec + matmul(me%Conder(1:cos,1:me%Numvar),me%Vardes(1:me%Numvar))
                  me%Conval = me%Conval + convec*fac
                  cycle main

            end select

            exit inner4
         enddo inner4

         exit main

      enddo main

      me%Funvar = desprv
      me%Confix = me%Conact(1:me%Numcon)
      if ( me%Senopt==-1 ) me%Senact = me%Conact(1:me%Numcon)

      deallocate (cosact)
      deallocate (varvec)
      deallocate (varwrk)
      deallocate (corvec)
      deallocate (desder)
      deallocate (desprv)
      deallocate (varprv)
      deallocate (convec)
      deallocate (conqua)
      deallocate (concor)

   end subroutine ogopti

   subroutine ogpwri(me,Objval,Numvio,Convio)

      !! WRITE OPTIMIZATION LOG IN PYGMO FORMAT
      !!
      !! 2023/01/25 | W. MARTENS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Objval !! OBJECTIVE VALUE
      integer(ip),intent(in) :: Numvio !! NUMBER OF VIOLATED CONSTRAINTS
      real(wp),intent(in) :: Convio !! TOTAL CONSTRAINT VIOLATION

      character(len=2) :: feas

      if ( me%Verbos==0 ) return
      ! Print header
      if ( me%Fevals==0 ) call me%ogpwri_start()
      ! Increase counter for cost function evaluations
      me%Fevals = me%Fevals + 1
      ! Every 50 lines print the column names.
      if ( mod(real(me%Fevals-1.0_wp)/real(me%Verbos),50.0_wp)==0.0_wp ) &
         write (me%Loglup,'(A10,A15,A15,A15,A2)') "objevals:" , "objval:" , "violated:" , "viol. norm:"
      if ( me%Verbos/=0 .and. mod(me%Fevals,me%Verbos)==0.0_wp ) then
         if ( Convio>0.0_wp ) then
            feas = " i"
         else
            feas = "  "
         endif

         ! Write the log line (different format depending on violation size)
         if ( Convio==0.0_wp ) then
            write (me%Loglup,'(I10,F15.4,I15,I15,A2)') me%Fevals , Objval , Numvio , int(Convio) , feas
         elseif ( Convio>1.0e-3_wp ) then
            write (me%Loglup,'(I10,F15.4,I15,F15.6,A2)') me%Fevals , Objval , Numvio , Convio , feas
         else
            write (me%Loglup,'(I10,F15.4,I15,E15.6,A2)') me%Fevals , Objval , Numvio , Convio , feas
         endif
      endif

      ! Write final summary
      if ( me%Pygfla/=0 ) call me%ogpwri_end(Objval,Numvio,Convio)

   end subroutine ogpwri

   subroutine ogpwri_end(me,Objval,Numvio,Convio)

      !! WRITE OPTIMIZATION END RESULT IN PYGMO FORMAT
      !!
      !! 2023/01/25 | W. MARTENS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Objval !! OBJECTIVE VALUE
      real(wp),intent(in) :: Convio !! TOTAL CONSTRAINT VIOLATION
      integer(ip),intent(in) :: Numvio !! NUMBER OF VIOLATED CONSTRAINTS

      if ( me%Pygfla==0 ) return

      ! Write termination message
      write (me%Loglup,'("")')
      write (me%Loglup,'("Final values after iteration        ", I10:)')  me%Numite
      write (me%Loglup,'("Final objective value:              ", F10.4)') Objval
      write (me%Loglup,'("Final constraint violation:         ", F10.4)') Convio
      write (me%Loglup,'("Final num. of violated constraints: ", I10)')   Numvio

      select case (me%Pygfla)
       case ( 1 ); write (me%Loglup,'("Successful termination: Optimal solution found.")')
       case ( 2 ); write (me%Loglup,'("Successful termination: Constraints matched.")')
       case ( 3 ); write (me%Loglup,'("Not converged.")')
       case ( 4 ); write (me%Loglup,'("Problem appears infeasible.")')
      end select
      write (me%Loglup,'("")')

   end subroutine ogpwri_end

   subroutine ogpwri_start(me)

      !! WRITE OPTIMIZATION LOG IN PYGMO FORMAT
      !!
      !! 2023/01/25 | W. MARTENS | NEW

      class(optgra),intent(inout) :: me

      write (me%Loglup,'("OPTGRA plugin for pagmo/pygmo:")')
      select case (me%Varder) ! DERIVATIVES COMPUTATION MODE
       case ( 0 );      write (me%Loglup,'("")') ! VALUES ONLY
       case ( 1, -1 );  write (me%Loglup,'("    User-defined gradients")')
       case ( 2 );      write (me%Loglup,'("    Numerical gradients by double differencing")')
       case ( 3 );      write (me%Loglup,'("    Numerical gradients by single differencing")')
      end select

      select case (me%Optmet)
       case ( 0 ); write (me%Loglup,'("    Steepest descent method")')
       case ( 1 ); write (me%Loglup,'("    Modified spectral conjugate gradient method")')
       case ( 2 ); write (me%Loglup,'("    Spectral conjugate gradient method")')
       case ( 3 ); write (me%Loglup,'("    Conjugate gradient method")')
      end select

      write (me%Loglup,'("")')

   end subroutine ogpwri_start

   subroutine ogrigt(me,Actinp,Actout)

      !! RIGHT-MULTIPLIES VECTOR LOWER TRIANGULAR MATRIX OBTAINED BY REDUCTION
      !! AND SUBSEQUENT INVERSION OF DERIVATIVES OF ACTIVE CONSTRAINTS
      !!
      !! 2008/01/16 | J. SCHOENMAEKERS | NEW

      class(optgra),intent(inout) :: me
      real(wp),intent(in) :: Actinp(me%Numcon) !! VECTOR INITAL
      real(wp),intent(inout) :: Actout(me%Numcon) !! VECTOR FINAL (MAY BE SAME AS ACTINP)

      integer(ip) :: row , col , act
      real(wp) :: val

      do col = me%Numact , 1 , -1
         val = Actinp(col)
         do act = me%Numact , col + 1 , -1
            row = me%Actcon(act)
            val = val - me%Conred(row,col)*Actout(act)
         enddo
         row = me%Actcon(col)
         Actout(col) = val/me%Conred(row,col)
      enddo

   end subroutine ogrigt

   subroutine ogsens(me,Consta,Concon,Convar,Varcon,Varvar)

      !! NEAR-LINEAR OPTIMIZATION TOOL SENSITIVITY ANALYSIS
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
      do act = 1 , me%Numact
         con = me%Actcon(act)
         Consta(con) = 1
      enddo

      ! SENSITIVITY OF CONTRAINTS W.R.T. ACTIVE CONSTRAINTS
      Concon = 0.0_wp
      do con = 1 , me%Numcon + 1
         if ( me%Conact(con)>0 ) Concon(con,con) = 1.0_wp
         if ( me%Conact(con)>0 ) cycle
         me%Conref = me%Conred(con,1:me%Numact)
         call me%ogrigt(me%Conref,me%Conref)
         do act = 1 , me%Numact
            ind = me%Actcon(act)
            Concon(con,ind) = -me%Conref(act)
         enddo
      enddo

      ! SENSITIVITY OF CONSTRAINTS W.R.T. PARAMETERS
      Convar = 0.0_wp
      do con = 1 , me%Numcon + 1
         if ( me%Conact(con)>0 ) cycle
         do var = 1 , me%Numvar
            if ( me%Vartyp(var)==0 ) cycle
            val = me%Sender(con,var)
            do act = 1 , me%Numact
               ind = me%Actcon(act)
               val = val + Concon(con,ind)*me%Sender(ind,var)
            enddo
            Convar(con,var) = val
         enddo
      enddo

      ! SENSITIVITY OF VARIABLES W.R.T. ACTIVE CONSTRAINTS
      Varcon = 0.0_wp
      do var = 1 , me%Numvar
         if ( me%Vartyp(var)/=0 ) cycle
         do act = 1 , me%Numact
            con = me%Actcon(act)
            me%Conref(act) = me%Conder(con,var)
         enddo
         call me%ogleft(me%Conref,me%Conref)
         call me%ogrigt(me%Conref,me%Conref)
         do act = 1 , me%Numact
            con = me%Actcon(act)
            Varcon(var,con) = -me%Conref(act)
         enddo
      enddo

      ! SENSITIVITY OF VARIABLES W.R.T. PARAMETERS
      Varvar = 0.0_wp
      do par = 1 , me%Numvar
         Varvar(par,par) = 1.0_wp
         if ( me%Vartyp(par)/=1 ) cycle
         do var = 1 , me%Numvar
            if ( me%Vartyp(var)/=0 ) cycle
            val = 0.0_wp
            do act = 1 , me%Numact
               con = me%Actcon(act)
               val = val + Varcon(var,con)*me%Sender(con,par)
            enddo
            Varvar(var,par) = val
         enddo
      enddo

      ! DESCALE SENSITIVITY
      do con = 1 , me%Numcon + 1
         typ = me%Contyp(con)
         sca = me%Consca(con)
         if ( typ<0 ) sca = -sca
         Convar(con,1:me%Numvar) = Convar(con,1:me%Numvar)*sca
         Concon(con,1:me%Numcon) = Concon(con,1:me%Numcon)*sca
         if ( con>me%Numcon ) cycle
         Varcon(1:me%Numvar,con) = Varcon(1:me%Numvar,con)/sca
         Concon(1:me%Numcon+1,con) = Concon(1:me%Numcon+1,con)/sca
      enddo

      do var = 1 , me%Numvar
         sca = me%Varsca(var)
         Varcon(var,1:me%Numcon) = Varcon(var,1:me%Numcon)*sca
         Varvar(var,1:me%Numvar) = Varvar(var,1:me%Numvar)*sca
         Convar(1:me%Numcon+1,var) = Convar(1:me%Numcon+1,var)/sca
         Varvar(1:me%Numvar,var) = Varvar(1:me%Numvar,var)/sca
      enddo

   end subroutine ogsens

   subroutine ogssst(me,Varsen,Quasen,Consen,Actsen,Dersen,Actsav,Consav,Redsav,Dersav,Actnum)

      !! NEAR-LINEAR OPTIMIZATION TOOL SENSITIVITY ANALYSIS
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

      do var = 1 , me%Numvar
         me%Senvar(var) = Varsen(var)
      enddo

      do con = 1 , me%Numcon + 1
         me%Senqua(con) = Quasen(con)
         me%Sencon(con) = Consen(con)
         me%Senact(con) = Actsen(con)
         do var = 1 , me%Numvar
            me%Sender(con,var) = Dersen(con,var)
         enddo
      enddo

      ! Temporary status saved of which constraints are active
      do con = 1 , me%Numcon + 1
         me%Actcon(con) = Actsav(con)
      enddo

      do con = 1 , me%Numcon + 4
         me%Conact(con) = Consav(con)
      enddo

      do con = 1 , me%Numcon + 3
         do var = 1 , me%Numvar
            me%Conred(con,var) = Redsav(con,var)
         enddo
      enddo

      do con = 1 , me%Numcon
         do var = 1 , me%Numvar
            me%Conder(con,var) = Dersav(con,var)
         enddo
      enddo
   end subroutine ogssst

   subroutine ogwrit(me,Lev,Str)

      !! Write a meessage to the log.
      !!
      !! 2014/07/29 | J. SCHOENMAEKERS | NEW
      !! 2025/02/09 | J. Williams | added trim()

      class(optgra),intent(inout) :: me
      character(len=*),intent(in) :: Str !! string to print
      integer(ip),intent(in) :: Lev !! only print if `Loglev` is >= this

      if ( Lev<=me%Loglev ) then
         write (me%Loglun,'(A)') trim(Str)
         flush (me%Loglun)
      endif

   end subroutine ogwrit

   pure subroutine mul2m(a1,m1,k1,l1,n1,a2,m2,k2,l2,n2,a,m,k,l,n)

      !! Matrix multiply.
      !!
      !! `A(K:K+N1,L:L+N) = A1(K1:K1+N1,L1:L1+N2) * A2(K2:K2+N2,L2:L2+N3)`

      integer,intent(in) :: m1, m2, m, k, k1, k2, l, l1 , l2 , n , n1 , n2
      real(wp),intent(out) :: a(m,*)
      real(wp),intent(in) :: a1(m1,*)
      real(wp),intent(in) :: a2(m2,*)

      real(wp) :: f1 , f2
      integer(ip) :: i , i1 , i2 , ic , ir

      do i1 = k , k + n1 - 1
         do i = l , l + n - 1
            a(i1,i) = 0.0_wp
         enddo
      enddo

      do i1 = 0 , n1 - 1
         do i2 = 0 , n2 - 1
            if ( k1>=0 ) then
               f1 = a1(i1+k1,i2+l1)
            else
               f1 = a1(i2-k1,i1+l1)
            endif
            if ( f1/=0.0_wp ) then
               do i = 0 , n - 1
                  if ( k2>=0 ) then
                     f2 = a2(i2+k2,i+l2)
                  else
                     f2 = a2(i-k2,i2+l2)
                  endif
                  if ( f2/=0.0_wp ) then
                     f2 = f2*f1
                     ic = i1 + k
                     ir = i + l
                     a(ic,ir) = a(ic,ir) + f2
                  endif
               enddo
            endif
         enddo
      enddo

   end subroutine mul2m

   pure subroutine mulvs(x,a,z,Kd)
      !! Scalar Vector multiply.
      !!
      !! `Z (1:KD) = X (1:KD) * A`

      real(wp),intent(in) :: a !! SCALAR
      real(wp),intent(in) :: x(*) !! VECTOR
      real(wp),intent(out) :: z(*) !! VECTOR
      integer(ip),intent(in) :: Kd !! NUMBER OF ELEMENTS TO BE USED
      integer(ip) :: i

      do i = 1 , Kd
         z(i) = x(i)*a
      enddo
   end subroutine mulvs

   pure subroutine sum2v(v1,v2,v,k)
      !! Vector addition.
      !!
      !! `V(1:K) = V1(1:K) + V2(1:K)`

      real(wp),intent(in) :: v1(*) , v2(*)
      real(wp),intent(out) :: v(*)
      integer(ip),intent(in) :: k
      integer(ip) :: i

      do i = 1 , k
         v(i) = v1(i) + v2(i)
      enddo
   end subroutine sum2v

!****************************************************************************************************
end module optgra_module
!****************************************************************************************************
