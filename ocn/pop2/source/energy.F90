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
use grid
use constant_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_sys_mod
use domain
use blocks

      IMPLICIT NONE

      real (r8) :: ek0,ea0,eb0,et0,es0
      REAL (r8) :: EK,EA,EB,ET,ES
      real(r8)::t0,t1,clock0f,nst
      save t0
      integer :: nnn,iblock

      type (block) :: this_block

!---------------------------------------------------------------------
!     TOTAL K.E.
!---------------------------------------------------------------------

      EK = 0.0
      month=(iyfm-1)*12+mon0

      EB = 0.0
      nst=0

      EK = 0.0
      month=(iyfm-1)*12+mon0

!!$OMP PARALLEL DO PRIVATE (K,J,I),reduction(+:EK)
      do iblock= 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
               EK = EK + DZP (K)* uarea(i,j,iblock)* VIV (I,J,K,iblock) &
                    * (U (I,J,K,iblock)* U (I,J,K,iblock) + V (I,J,K,iblock)* V (I,J,K,iblock))
            END DO
         END DO
      END DO
      END DO

      call mpi_reduce(ek,ek0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      ek0= ek0*0.5D0* d0

!---------------------------------------------------------------------
!     AVAILABLE P.E.
!---------------------------------------------------------------------

      EA = 0.0
!!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:EA)
      do iblock= 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
            EA = EA + tarea(i,j,iblock)* VIT (I,J,1,iblock)* H0 (I,J,iblock)* H0 (I,J,iblock)
         END DO
      END DO
      END DO

      call mpi_reduce(ea,ea0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      ea0= ea0*0.5D0*d0*g

      EB = 0.0
      nst=0
!!$OMP PARALLEL DO PRIVATE (J,I),reduction(+:EB)
      do iblock= 1, nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
            EB = EB + tarea(i,j,iblock)* VIT (I,J,1,iblock)* AT (I,J,1,1,iblock)
         END DO
      END DO
      end do

      call mpi_reduce(eb,eb0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      eb0= eb0/area_t

!---------------------------------------------------------------------
!     GLOBAL MEAN TEMPERATURE & SALINITY
!---------------------------------------------------------------------

      ET = 0.0D0
      ES = 0.0D0
!!$OMP PARALLEL DO PRIVATE (K,J,I),reduction(+:ET,ES)
      DO iblock= 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         DO J = this_block%jb, this_block%je
            DO I = this_block%ib, this_block%ie
               ET = ET + DZP (K)*tarea(i,j,iblock)* VIT (I,J,K,iblock)* AT(I,J,K,1,iblock)
               ES = ES + DZP (K)*tarea(i,j,iblock)* VIT (I,J,K,iblock)* AT(I,J,K,2,iblock)
            END DO
         END DO
      END DO
      END DO

      call mpi_reduce(et,et0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_reduce(es,es0,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      et0= et0/ volume_t
      es0= es0/ volume_t *1.0D3

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


