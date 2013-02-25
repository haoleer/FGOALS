!  CVS: $Id: energy.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ENERGY
!     =================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_sys_mod

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0,et0,es0
      REAL (r8) :: EK,EA,EB,ET,ES
      real(r8)::t0,t1,clock0f,nst
      save t0
      integer :: nnn


!---------------------------------------------------------------------
!     TOTAL K.E.
!---------------------------------------------------------------------

      EK = 0.0
      month=(iyfm-1)*12+mon0

      EB = 0.0
      nst=0
!---------------------------------------------------------------------
!     Calculating CPU time for per-day.
!---------------------------------------------------------------------
      if (mytid==0) then
         t1=clock0f()
      endif

      if (mytid==0 )then
      WRITE (6,FMT='(I5,I3,6D25.15)') MONTH,IDAY,EK0,EA0,EB0,ET0,ES0,t1-t0
!      call flush_(6)
      end if

      if (mytid==0) then
         t0=t1
      endif
!
      RETURN
      END SUBROUTINE ENERGY


      SUBROUTINE chk_var3d(var,ch,a)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0
      REAL (r8) :: EK,EA,EB
      real(r8)  :: var(imt,jmt,km)
      integer :: a
      character :: ch*2

      EA = 0.0
      EB = 0.0
!
      if (mytid==0 )then
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
      end if
!
      RETURN
      END SUBROUTINE chk_var3d

      SUBROUTINE chk_var3d1(var,ch,a)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0
      REAL (r8) :: EK,EA,EB
      real(r8)  :: var(imt,jmt,km)
      integer :: a
      character :: ch*2

      EA = 0.0
      EB = 0.0
      if (mytid==0 )then
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
      end if
!
      RETURN
      END SUBROUTINE chk_var3d1

      SUBROUTINE chk_var1(var,ch,a)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0
      REAL (r8) :: EK,EA,EB
      real(r8)  :: var(imt,jmt,km)
      integer :: a
      character :: ch*2

      EB = 0.0

      EK = 0.0
      if (mytid==0 )then
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
      end if

      RETURN
      END SUBROUTINE chk_var1


      SUBROUTINE chk_var2d(var,ch,a)
!     =======================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif

      IMPLICIT NONE

      REAL (r8) :: EK,EA,EB
      REAL (r8) :: EK0,EA0,EB0
      REAL(r8)  :: var(imt,jmt)
      integer :: a
      character :: ch*2

      EA = 0.0
      EB = 0.0
      if (mytid==0)then
      WRITE (6,FMT='(a2,I5,I3,D25.15)') ch,MONTH,IDAY,EK0
      end if

      RETURN
      END SUBROUTINE chk_var2d


