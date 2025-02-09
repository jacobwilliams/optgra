SUBROUTINE ogclos()
! ======================================================================
! NEAR-LINEAR OPTIMISATION TOOL TAILORED FOR S/C TRAJECTORY DESIGN:
! DEALLOCATION OF ARRAYS
! ======================================================================
! 2008/01/16 | J. SCHOENMAEKERS | NEW
! ======================================================================
   IMPLICIT NONE
! ======================================================================
   INCLUDE "ogdata.inc"
! ======================================================================
! VARIABLES
! ----------------------------------------------------------------------
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
! ======================================================================
! CONSTRAINTS
! ----------------------------------------------------------------------
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
! ======================================================================
! DERIVATIVES
! ----------------------------------------------------------------------
   DEALLOCATE (Varper)
! ======================================================================
! WORKING VECTORS
! ----------------------------------------------------------------------
   DEALLOCATE (Actcon)
   DEALLOCATE (Confix)
   DEALLOCATE (Conact)
   DEALLOCATE (Conder)
   DEALLOCATE (Conred)
   DEALLOCATE (Sender)
   DEALLOCATE (Conopt)
! ======================================================================
END SUBROUTINE ogclos
