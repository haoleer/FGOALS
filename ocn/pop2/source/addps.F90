!  CVS: $Id: addps.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ADDPS
!     ================
!     COMPENSATING THE LOSS OF GROSS MASS
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use domain
      IMPLICIT NONE
      REAL(r8)    :: ERROR,DH00,error0
      integer     :: iblock
      ERROR = 0.0D0
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I), shared(dh00)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JMT
         DO I = 1,IMT
            H0 (I,J,IBLOCK)= (H0 (I,J,IBLOCK) + DH00)* VIT (I,J,1,IBLOCK)
         END DO
      END DO
   END DO
 
      RETURN
      END SUBROUTINE ADDPS
 
