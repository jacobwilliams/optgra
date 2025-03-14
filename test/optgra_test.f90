program optgra_test

    use optgra_module, wp => optgra_wp

    implicit none

    ! parameters:
    integer,parameter :: varnum = 1     !! number of variables
    integer,parameter :: connum = 1     !! number of constraints

    type(optgra) :: solver

    ! JW : some of these I don't know what to set to ...

    REAL(wp), dimension(connum+1), parameter :: Delcon = [0.0001_wp,0.0001_wp]    !! CONSTRAINTS DELTAS
    integer, dimension(connum+1), parameter :: Pricon = [1,1]               !! CONSTRAINTS PRIORITIES
    REAL(wp), dimension(connum+1), parameter :: Scacon = [1.0_wp,1.0_wp]    !! CONSTRAINTS CONVER THRESHOLD (1:NUMCON)
                                                                            !! MERIT       CONVER THRESHOLD (1+NUMCON)
    character(len=80), dimension(connum+1), parameter :: STRCON = ['f1','j ']    !! CONIABLES NAME STRING
    integer, dimension(connum+1), parameter :: LENCON = [2,1]                !! CONIABLES NAME LENGTH
    INTEGER, dimension(connum+1), parameter :: Typcon = [0,-1]    !! CONSTRAINTS TYPE (1:NUMCON)
                                                                  !! -> 1=GTE -1=LTE 0=EQU -2=DERIVED DATA
                                                                  !! MERIT       TYPE (1+NUMCON)
                                                                  !! -> 1=MAX -1=MIN
    INTEGER, parameter :: Dervar = 1    !! DERIVATIVES COMPUTATION MODE
                                        !! -> 1: USER DEFINED
                                        !! -> 2: NUMERIC WITH DOUBLE DIFFERENCING
                                        !! -> 3: NUMERIC WITH SINGLE DIFFERENCING
    real(wp), dimension(varnum), parameter :: Pervar = [0.00001_wp]  !! VARIABLES PERTURBATION FOR DERIVATIVES

    real(wp), parameter :: MAXVAR = 1.0_wp  !! MAXIMUM DISTANCE PER ITERATION
                                            !! -> SCALED
    real(wp), parameter :: SNDVAR = 1.0_wp !! PERTURBATION FOR 2ND ORDER DERIVATIVES
                                           !! -> SCALED
    integer, parameter :: Itemax = 100 !! MAXIMUM NUMBER OF ITERATIONS
    integer, parameter :: Itecor = 10 !! THESE AREN'T DEFINED IN THE DOCSTRING ...
    integer, parameter :: Iteopt = 10 !! THESE AREN'T DEFINED IN THE DOCSTRING ...
    integer, parameter :: Itediv = 10 !! THESE AREN'T DEFINED IN THE DOCSTRING ...
    integer, parameter :: Itecnv = 10 !! THESE AREN'T DEFINED IN THE DOCSTRING ...
    integer, parameter :: Metopt = 2 !! OPTIMIZATION METHOD
                                     !!
                                     !!  * 3: CONJUGATE GRADIENT METHOD
                                     !!  * 2: SPECTRAL CONJUGATE GRADIENT METHOD
                                     !!  * 1: MODIFIED SPECTRAL CONJUGATE GRADIENT METHOD
                                     !!  * 0: STEEPEST DESCENT METHOD
    integer, parameter :: Bosver = 1 !! VERBOSITY LEVEL FOR WRITING PYGMO LOG
                                     !!
                                     !!  *  0=NO OUTPUT
                                     !!  *  1 OUTPUT EVERY ITERATION
                                     !!  *  2 OUTPUT EVERY 2ND ITERATION
                                     !!  *  N OUTPUT EVERY NTH ITERATION
    integer, parameter :: Optsen = 0  !! SENSITIVITY OPTIMIZATION MODE
                                      !!
                                      !!  *  0: NO
                                      !!  * -1: INITIALIZATION
                                      !!  * +1: WITH CONSTRAINT CALCULATION
                                      !!  * +2: WITH CONSTRAINT BIAS
                                      !!  * +3: WITH CONSTRAINT CALC / NO OPTIM
                                      !!  * +4: WITH CONSTRAINT BIAS / NO OPTIM
    real(wp), dimension(varnum), parameter :: Scavar = [1.0_wp] !! VARIABLES SCALE FACTOR
    character(len=80), dimension(varnum), parameter :: STRVAR = ['x'] !! VARIABLES NAME STRING
    integer, dimension(varnum), parameter :: LENVAR = [1] !! VARIABLES NAME LENGTH
    INTEGER, dimension(varnum), parameter :: Typvar = [0] !! VARIABLES TYPE
                                                           !! -> 0=FREE VARIABLE
                                                           !! -> 1=PARAMETER FOR SENSITIVITY
    integer, parameter :: Levlog = 10 !! LEVEL OF LOG
                                     !! -> 0=NO OUTPUT
                                     !! -> 1<ALL
    integer, parameter :: Levmat = 0 !! Matlab LEVEL OF LOG     <---- JW: don't think this is used
                                     !! -> 0=NO OUTPUT
                                     !! -> 1<ALL
    integer, parameter :: Levtab = 1 !! LEVEL OF TAB
                                      !! -> 0=NO OUTPUT
                                      !! -> 1<ALL

    integer :: Luplog !! LOGICAL UNIT FOR WRITING PYGMO LOG
    integer :: Lunlog !! LOGICAL UNIT FOR WRITING LOG
    integer :: Luntab !! LOGICAL UNIT FOR WRITING TABLE

    ! variables:
    real(wp),dimension(varnum) :: Valvar !! VARIABLES VALUE
                                         !! -> NOT SCALED

    real(wp),dimension(connum+1) :: Valcon !! CONSTRAINTS VALUE (1:NUMCON)
                                            !! MERIT       VALUE (1+NUMCON)
                                            !! -> NOT SCALED
    integer :: finopt !! TERMINATION STATUS
                      !! -> 1=    MATCHED &     OPTIMAL
                      !! -> 2=    MATCHED & NOT OPTIMAL
                      !! -> 3=NOT MATCHED & NOT OPTIMAL
                      !! -> 4=NOT FEASIBL & NOT OPTIMAL
    integer :: Finite !! NOT DOCUMENTED ?

    write(*,*) 'optgra_test'

    open(newunit=Luplog, file = 'pygmo_log.txt', status='REPLACE')
    open(newunit=lunlog, file = 'log.txt', status='REPLACE')
    open(newunit=Luntab, file = 'table_log.txt', status='REPLACE')

    ! setup:
    write(*,*) 'initialize'
    call solver%initialize(Varnum,Connum,Calval,Calder,Delcon,Pricon,Scacon,Strcon,Lencon,&
                            Typcon,Dervar,Pervar,Maxvar,Sndvar,&
                            Itemax,Itecor,Iteopt,Itediv,Itecnv,&
                            Luplog,Bosver,Optsen,Scavar,&
                            Strvar,Lenvar,Typvar,Lunlog,Levlog,Levmat,&
                            Luntab,Levtab,Metopt)
    ! call oginit(Varnum,Connum)  ! i think this may also set defaults so some of below may be optional?

    write(*,*) 'set other params'
    ! call ogcdel(Delcon)
    ! call ogcpri(Pricon)
    ! call ogcsca(Scacon)
    ! call ogcstr(Strcon,Lencon)
    ! call ogctyp(Typcon)
    ! call ogderi(Dervar,Pervar)
    ! call ogdist(Maxvar,Sndvar)
    ! call ogiter(Itemax,Itecor,Iteopt,Itediv,Itecnv)
    ! call ogomet(Metopt)
    ! call ogplog(Luplog,Bosver)
    ! call ogsopt(Optsen)
    ! call ogvsca(Scavar)
    ! call ogvstr(Strvar,Lenvar)
    ! call ogvtyp(Typvar)
    ! call ogwlog(Lunlog,Levlog)
    ! call ogwmat(Levmat)
    ! call ogwtab(Luntab,Levtab)

    ! execute:
    write(*,*) 'execute'
    Valvar = [0.1_wp]
    call solver%solve(Valvar,Valcon,Finopt,Finite)

    write(*,*) ''
    write(*,*) 'solution:'
    write(*,*) ' * Valvar=',Valvar
    write(*,*) ' * Finopt=',Finopt
    write(*,*) ' * Finite=',Finite

    ! finish
    write(*,*) ''
    write(*,*) 'cleanup'
    call solver%destroy()
    close(lunlog)
    close(Luntab)
    close(Luplog)

    contains

    subroutine calval(me,varvec,Valcon,i)

        !! FUNCTION FOR VALUES
        !! INPUT AND OUTPUT NOT SCALED

        ! use ieee_arithmetic, only: ieee_is_finite, ieee_is_nan

        class(optgra),intent(inout) :: me
        real(wp), dimension(:), intent(in) :: varvec !! size is Numvar
        real(wp), dimension(:), intent(out) :: Valcon !! size is Numcon+1
        integer, intent(in) :: i !! JW: THIS IS NOT DOCUMENTED ?

        real(wp) :: x, f, j

        ! write(*,*) 'calval : varvec = ', varvec, 'i=',i

        x = varvec(1)
        f = 10000.0_wp * sin(x)
        j = x

        Valcon = [f,j]

        ! write(*,*) 'calval : Valcon = ', Valcon

    end subroutine calval

    subroutine Calder(me,varvec,convec,Dercon)
        !! FUNCTION FOR VALUES AND DERIVATIVES
        !! INPUT AND OUTPUT NOT SCALED
        class(optgra),intent(inout) :: me
        real(wp), dimension(:), intent(in) :: varvec !! size is Numvar
        real(wp), dimension(:), intent(out) :: convec !! size is Numcon+1
        real(wp), dimension(:,:), intent(out) :: Dercon !! size is Numcon+1,Numvar

        real(wp) :: x

        x = varvec(1)
        call calval(me,varvec,convec,0)
        Dercon(1,1) = 10000.0_wp * cos(x)    ! df/dx
        Dercon(2,1) = 1.0_wp    ! dj/dx

    end subroutine Calder


end program optgra_test