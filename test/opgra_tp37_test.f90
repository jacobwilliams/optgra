program opgra_tp37_test

!! this program executes the example problem from
!! https://mdolab-pyoptsparse.readthedocs-hosted.com/en/latest/quickstart.html
!!
!! this is schittkowski tp37 test problem.
!! see also: https://klaus-schittkowski.de/test_problems.pdf

use optgra_module, wp => optgra_wp
use pyplot_module

implicit none

! parameters:
integer,parameter :: varnum = 3  !! number of variables
integer,parameter :: connum = 2 + 6  !! number of constraints (2 constraints + 4 bounds)

type(optgra) :: solver

! jw : some of these i don't know what to set to ...
!
! not sure what the difference is between delcon and scacon ???
! ... delcon doesn't seem to have must effect on the solution ...

real(wp), dimension(connum+1), parameter :: delcon = [1.0_wp,1.0_wp,&
                                                      1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp,&
                                                      1.0_wp]    !! constraints deltas
integer, dimension(connum+1), parameter :: pricon = [1,1,0,0,0,0,0,0,1]   !! constraints priorities
real(wp), dimension(connum+1), parameter :: scacon = [1.0e-9_wp,1.0e-9_wp,&
                                                      1.0e-9_wp,1.0e-9_wp,1.0e-9_wp,1.0e-9_wp,1.0e-9_wp,1.0e-9_wp,&
                                                      1.0e-4_wp]    !! constraints conver threshold (1:numcon)
                                                                    !! merit       conver threshold (1+numcon)
character(len=80), dimension(connum+1), parameter :: strcon = ['f1 ','f2 ', &
                                                               'x1l', 'x1u', 'x2l', 'x2u', 'x3l', 'x3u', &
                                                               'j  ']    !! coniables name string
integer, dimension(connum+1), parameter :: lencon = [2,2,3,3,3,3,3,3,1]   !! coniables name length
integer, dimension(connum+1), parameter :: typcon = [-1,-1,&
                                                      1,-1,1,-1,1,-1,&
                                                     -1] !!  * constraints type `(1:numcon)`
                                                         !!    *  1 : >= constraint
                                                         !!    * -1 : <= constraint
                                                         !!    *  0 : == constraint
                                                         !!    * -2 : derived data
                                                         !!
                                                         !!  * merit type `(1+numcon)`
                                                         !!    *  1 : max
                                                         !!    * -1 : min
integer, parameter :: dervar = 1    !! derivatives computation mode
                                    !! -> 1: user defined
                                    !! -> 2: numeric with double differencing
                                    !! -> 3: numeric with single differencing
real(wp), dimension(varnum), parameter :: pervar = [0.00001_wp,0.00001_wp,0.00001_wp]  !! variables perturbation for derivatives

real(wp), parameter :: maxvar = 1.0_wp  !! maximum distance per iteration
                                        !! -> scaled
real(wp), parameter :: sndvar = 1.0_wp !! perturbation for 2nd order derivatives
                                       !! -> scaled
integer, parameter :: itemax = 100 !! maximum number of iterations
integer, parameter :: itecor = 10 !! these aren't defined in the docstring ...
integer, parameter :: iteopt = 10 !! these aren't defined in the docstring ...
integer, parameter :: itediv = 10 !! these aren't defined in the docstring ...
integer, parameter :: itecnv = 10 !! these aren't defined in the docstring ...
integer, parameter :: metopt = 2 !! optimization method
                                 !!
                                 !!  * 3: conjugate gradient method
                                 !!  * 2: spectral conjugate gradient method
                                 !!  * 1: modified spectral conjugate gradient method
                                 !!  * 0: steepest descent method
integer, parameter :: bosver = 1 !! verbosity level for writing pygmo log
                                 !!
                                 !!  *  0=no output
                                 !!  *  1 output every iteration
                                 !!  *  2 output every 2nd iteration
                                 !!  *  n output every nth iteration
integer, parameter :: optsen = 0  !! sensitivity optimization mode
                                  !!
                                  !!  *  0: no
                                  !!  * -1: initialization
                                  !!  * +1: with constraint calculation
                                  !!  * +2: with constraint bias
                                  !!  * +3: with constraint calc / no optim
                                  !!  * +4: with constraint bias / no optim
real(wp), dimension(varnum), parameter :: scavar = [20.0_wp,20.0_wp,20.0_wp] !! variables scale factor
character(len=80), dimension(varnum), parameter :: strvar = ['x1','x2','x3'] !! variables name string
integer, dimension(varnum), parameter :: lenvar = [2,2,2] !! variables name length
integer, dimension(varnum), parameter :: typvar = [0,0,0] !! variables type
                                                          !! -> 0=free variable
                                                          !! -> 1=parameter for sensitivity
integer, parameter :: levlog = 10 !! level of log
                                  !! -> 0=no output
                                  !! -> 1<all
integer, parameter :: levmat = 0 !! matlab level of log     <---- jw: don't think this is used
                                 !! -> 0=no output
                                 !! -> 1<all
integer, parameter :: levtab = 1 !! level of tab
                                 !! -> 0=no output
                                 !! -> 1<all

integer :: luplog !! logical unit for writing pygmo log
integer :: lunlog !! logical unit for writing log
integer :: luntab !! logical unit for writing table

! variables:
real(wp),dimension(varnum) :: valvar !! variables value
                                     !! -> not scaled

real(wp),dimension(connum+1) :: valcon  !! constraints value (1:numcon)
                                        !! merit       value (1+numcon)
                                        !! -> not scaled
integer :: finopt !! termination status
                  !! -> 1=    matched &     optimal
                  !! -> 2=    matched & not optimal
                  !! -> 3=not matched & not optimal
                  !! -> 4=not feasibl & not optimal
integer :: finite !! not documented ?

type(pyplot) :: plt

real(wp),dimension(:),allocatable :: iter_hist
real(wp),dimension(:),allocatable :: g1_hist
real(wp),dimension(:),allocatable :: g2_hist
real(wp),dimension(:),allocatable :: obj_hist

! the bounds will be cast as nonlinear constraints
!     vlb(:) = 0.0_wp  ! lower bounds
!     vub(:) = 42.0_wp ! upper bounds
!
! note pyoptgra has a bounds_to_constraints option
! to convert bounds to constraints.

open(newunit=luplog, file = 'pygmo_log.txt', status='replace')
open(newunit=lunlog, file = 'log.txt', status='replace')
open(newunit=luntab, file = 'table_log.txt', status='replace')

! setup:
write(*,*) 'initialize'
call solver%initialize(varnum,connum,calval,calder,delcon,pricon,scacon,strcon,lencon,&
                        typcon,dervar,pervar,maxvar,sndvar,&
                        itemax,itecor,iteopt,itediv,itecnv,&
                        luplog,bosver,optsen,scavar,&
                        strvar,lenvar,typvar,lunlog,levlog,levmat,&
                        luntab,levtab,metopt)

valvar = [10.0_wp, 10.0_wp, 10.0_wp] !! initial values
call solver%solve(valvar,valcon,finopt,finite)

write(*,*) ''
write(*,*) 'solution:'
write(*,*) ' * valvar=',valvar
write(*,*) ' * finopt=',finopt
write(*,*) ' * finite=',finite

! finish
write(*,*) ''
write(*,*) 'cleanup'
call solver%destroy()
close(lunlog)
close(luntab)
close(luplog)

! make a plot of the constraint iteration history
call plt%initialize(grid=.true.,xlabel='iteration',ylabel='constraint violation',&
                    title='tp37 constraints',legend=.true.)
call plt%add_plot(iter_hist,g1_hist,label='g(1)',linestyle='b.-',markersize=5,linewidth=2)
call plt%add_plot(iter_hist,g2_hist,label='g(2)',linestyle='r.-',markersize=5,linewidth=2)
call plt%savefig('iterations.png')

call plt%initialize(grid=.true.,xlabel='iteration',ylabel='objective function value',&
                    title='tp37 objective function',legend=.true.)
call plt%add_plot(iter_hist,obj_hist,label='obj',linestyle='k.-',markersize=5,linewidth=2)
call plt%savefig('obj.png')

contains

    subroutine calval(me,x,valcon,i)

        !! function for values
        !! input and output not scaled

        class(optgra),intent(inout) :: me
        real(wp), dimension(:), intent(in) :: x !! size is numvar
        real(wp), dimension(:), intent(out) :: valcon !! size is numcon+1
        integer, intent(in) :: i !! jw: this is not documented ?

        real(wp) :: f(8), j

        ! constraint values
        f(1) =  x(1) + 2.0_wp * x(2) + 2.0_wp * x(3) - 72.0_wp  ! <= 0
        f(2) = -x(1) - 2.0_wp * x(2) - 2.0_wp * x(3)            ! <= 0

        ! bounds, cast as constraints:
        !   0  <= x1  --> x1      >= 0
        !   x1 <= 42  --> x1 - 42 <= 0
        !   0  <= x2  --> x2      >= 0
        !   x2 <= 42  --> x2 - 42 <= 0
        f(3) = x(1)
        f(4) = x(1) - 42.0_wp
        f(5) = x(2)
        f(6) = x(2) - 42.0_wp
        f(7) = x(3)
        f(8) = x(3) - 42.0_wp

        ! objective function
        j = -x(1) * x(2) * x(3)

        valcon = [f(1:connum),j]

    end subroutine calval

    subroutine calder(me,x,convec,dercon)
        !! function for values and derivatives
        !! input and output not scaled
        class(optgra),intent(inout) :: me
        real(wp), dimension(:), intent(in) :: x !! size is numvar
        real(wp), dimension(:), intent(out) :: convec !! size is numcon+1
        real(wp), dimension(:,:), intent(out) :: dercon !! size is numcon+1,numvar

        ! for now, just call this an iteration.... but not sure if that is always the case ?
        integer,save :: i = 0 !! "iteration" counter
        i = i + 1

        call calval(me,x,convec,0)

        write(*,'(a,i2,a,3(f10.6,1x),a,8(f10.6,1x),f15.6)') 'iter=', i, ' x:', x(1:varnum), ' f:', convec(1:connum+1)

        ! constraint gradient:
        dercon(1,:) = [1.0_wp,   2.0_wp,  2.0_wp]
        dercon(2,:) = [-1.0_wp, -2.0_wp, -2.0_wp]

        dercon(3,:) = [1.0_wp, 0.0_wp, 0.0_wp]  ! for the bounds
        dercon(4,:) = [1.0_wp, 0.0_wp, 0.0_wp]
        dercon(5,:) = [0.0_wp, 1.0_wp, 0.0_wp]
        dercon(6,:) = [0.0_wp, 1.0_wp, 0.0_wp]
        dercon(7,:) = [0.0_wp, 0.0_wp, 1.0_wp]
        dercon(8,:) = [0.0_wp, 0.0_wp, 1.0_wp]

        ! objective gradient:
        dercon(9,:) = [-x(2) * x(3), &
                       -x(1) * x(3), &
                       -x(1) * x(2)]

        ! for the plot:
        if (.not. allocated(iter_hist)) then
            allocate(iter_hist(0))
            allocate(g1_hist(0))
            allocate(g2_hist(0))
            allocate(obj_hist(0))
        end if
        iter_hist = [iter_hist, real(i,wp)]
        g1_hist   = [g1_hist,   convec(1) ]
        g2_hist   = [g2_hist,   convec(2) ]
        obj_hist  = [obj_hist,  convec(connum+1) ]

    end subroutine calder

end program opgra_tp37_test