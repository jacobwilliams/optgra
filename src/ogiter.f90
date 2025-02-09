SUBROUTINE ogiter(Itemax,Itecor,Iteopt,Itediv,Itecnv)
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEFINE OPTIMISATION CONTROL PARAMETERS
! ======================================================================
! INP | ITEMAX           | I*4 | MAXIMUM NUMBER OF ITERATIONS
! ----------------------------------------------------------------------
! INP | MAXVAR           | R*8 | MAXIMUM DISANCE PER ITERATION
!     |                  |     | -> SCALED
! ----------------------------------------------------------------------
! INP | SNDVAR           | R*8 | PERTURBATION FOR SND ORDER DERIVATIVES
!     |                  |     | -> SCALED
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
   INTEGER(4) Itemax , Itecor , Iteopt , Itediv , Itecnv
! ======================================================================
   Maxite = Itemax
   Corite = Itecor
   Optite = Iteopt
   Divite = Itediv
   Cnvite = Itecnv
   IF ( Corite>Maxite ) Corite = Maxite
   IF ( Optite>Maxite ) Optite = Maxite
   IF ( Divite>Corite ) Divite = Corite
   IF ( Cnvite>Optite ) Cnvite = Optite
! ======================================================================
END SUBROUTINE ogiter
