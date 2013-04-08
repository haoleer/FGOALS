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
use cdf_mod, only : start1,count1,start2,count2,start3,count3,start4,count4
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
      IMPLICIT NONE
#include <netcdf.inc>

!
!     Define Variables.
      integer*4   :: ncid, iret
!
      REAL(r4)    :: ZKP_IN (KMP1)
!      REAL(r8)    :: ZKP (KMP1)
      REAL(r8)    :: OMEGA,DX,ABCD,YU,CURU,EPS
      REAL(r8)    :: rpart,efold1,efold2,swarg1,swarg2
      REAL(r8)    :: AJQ,rscl1,rscl2
      INTEGER :: IREC



      if (mytid==0)then
      write(6,*)"Beginning------GRIDS !"
      endif

!XC
#if (defined TSPAS)
      if (mytid==0)then
      write(6,*)"Use TSPAS advection scheme!"
      endif
#else
      if (mytid==0)then
      write(6,*)"Use CTCS advection scheme!"
      endif
#endif
!XC

!--------------------------------------------------------------
!     SET LOCAL CONSTANTS
!--------------------------------------------------------------
!     AM_TRO=2.0E+3
!     AM_EXT=2.0E+5
!     DLAM  =0.5
      RADIUS = 6371000.0D0
      OMEGA = 0.7292D-4

!--------------------------------------------------------------
!     ZONAL RESOULTION
!--------------------------------------------------------------

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

!      if (mytid==0) print*,pen

      DO k = 1,kmm1
!lhl        pen(k) = pen(k)*solar0*OD0CP
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


