!  CVS: $Id: grids.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      SUBROUTINE GRIDS
!     ================
!     TOPOGRAPHY & GRIDS
!-----------------------------------------------------------------------
!
! Purpose: Set up some constants dependent on model grids.
!
! Author: Yongqiang Yu and Hailong Liu, Dec. 31, 2002
!
!
!-----------------------------------------------------------------------


#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use pmix_mod
use work_mod
use constant_mod
use cdf_mod, only : start1,count1,start2,count2,start3,count3,start4,count4
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
      IMPLICIT NONE
#include <netcdf.inc>

!
      REAL(r8)    :: rpart,efold1,efold2,swarg1,swarg2
      REAL(r8)    :: AJQ,rscl1,rscl2


      if (mytid==0)then
      write(6,*)"Beginning------GRIDS !"
      endif

#if (defined SOLAR)
      rpart = 0.58D0
      efold1 = 0.35D0
      efold2 = 23.0D0
      rscl1 = 1.0D0/ efold1
      rscl2 = 1.0D0/ efold2

      DO k = 1,kmm1
         swarg1 = max (ZKP (k +1)* rscl1, -70.0)
         swarg2 = max (ZKP (k +1)* rscl2, -70.0)
         pen (k) = rpart * exp (swarg1) + (1.0D0- rpart)* exp (swarg2)
      END DO


      DO k = 1,kmm1
         pen (k) = pen (k)* OD0CP
      END DO

#endif

!--------------------------------------------------------------
!     compute boundary of area where vertical mixing coefficients

      if (mytid==0)then
      write(6,*)"END------------GRIDS !"
      endif
      RETURN
      END SUBROUTINE GRIDS


