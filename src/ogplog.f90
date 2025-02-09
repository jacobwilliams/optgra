SUBROUTINE ogplog(Luplog,Bosver)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEFINE WRITING IN PYGMO LOG FORMAT
! ======================================================================
! INP | LUPLOG           | I*4 | LOGICAL UNIT FOR WRITING PYGMO LOG
! ----------------------------------------------------------------------
! INP | BOSVER           | I*4 | VERBOSITY LEVEL
!     |                  |     | -> 0=NO OUTPUT
!     |                  |     | -> 1 OUTPUT EVERY ITERATION
!     |                  |     | -> 2 OUTPUT EVERY 2ND ITERATION
!     |                  |     | -> N OUTPUT EVERY NTH ITERATION
! ======================================================================
! 2023/01/25 | W. MARTENS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
! Yes, we stay true to the original vintage F77 style with our
! variable names just to confuse future developers :P
   INTEGER(4) Luplog
   INTEGER(4) Bosver
! ======================================================================
   Loglup = Luplog
   Verbos = Bosver
   Fevals = 0    ! initialize number of cost fun evaluations
   Pygfla = 0     ! pygmo output status flag: 0: continue iterating, 1: final output
! ======================================================================
END SUBROUTINE ogplog
