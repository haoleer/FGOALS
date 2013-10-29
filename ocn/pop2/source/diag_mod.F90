module diag_mod
#include <def-undef.h>
!-------------------------------------------------------------------------------
!
! Author: Yongqiang YU  ( 1  Dec, 2003)
!
!-------------------------------------------------------------------------------
use precision_mod
use param_mod 
use constant_mod
use pconst_mod 
use output_mod 
use tracer_mod
use work_mod 
use cdf_mod 
use msg_mod
use domain
use grid
use distribution
use gather_scatter
!
      implicit none
!
      contains
!

!     ===================
      SUBROUTINE msf
!     ===================
!
      implicit none
!
      real(r4) :: va(imt,jmt,km,max_blocks_clinic), vvv (imt,jmt,km, max_blocks_clinic)
      integer :: nin,nta(jmt_global), iblock
!
      allocate (work_1(imt,jmt,km,max_blocks_clinic),work_2(imt,jmt,km,max_blocks_clinic))
!
!
      DO  NIN=1,2
!
      IF(NIN.eq.1)then
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
      do iblock = 1, nblocks_clinic
      do k=1,km
      do j=1,jmt
      do i=1,imt
         work_1(i,j,k,iblock)=vit(i,j,k,iblock)
      enddo
      enddo
      enddo
      enddo
!
      ELSE
!$OMP PARALLEL DO PRIVATE (K,J,I)
      do iblock = 1, nblocks_clinic
      do k=1,km
      do j=1,jmt
      do i=1,imt
         if (basin(i,j,iblock)==1.or.basin(i,j,iblock)==2) then
            work_1(i,j,k,iblock)=vit(i,j,k,iblock)
         else
            work_1(i,j,k,iblock)=0.0
         end if
      enddo
      enddo
      enddo
      enddo
      ENDIF
!
      enddo
!
!
      psi = 0.0_r4
!
      DO K=2,km+1
         buffer_r4_local = work_1(:,:,k-1,:) * vsmon(:,:,k-1,:) *dzp(k-1) * dxu/float(imd)
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if (mytid == 0 ) then
            DO J=1,jmt_global
            DO I=1,imt_global
               psi(J,K,NIN)=psi(j,k,nin)+buffer_r4_global(I,J)*1.0E-6
            ENDDO
            ENDDO
            psi(J,K,NIN)=psi(J,K-1,NIN)+psi(J,K,NIN)
         end if
      ENDDO
!
      if (mytid == 0)  then
        where ( abs(psi(:,:,nin)) < 1.0e-25) psi(:,:,nin) = spval
      end if
!
      deallocate (work_1,work_2)
!
      RETURN
      END SUBROUTINE MSF

!
!====================================
      SUBROUTINE BAROSF
!====================================
      implicit none
!
      integer :: iblock
!
      buffer_r4_local= 0.0_r4
      buffer_r4_global= 0.0_r4
!
      do iblock = 1, nblocks_clinic
      do k=1,km
      do j=2,jmt-1
      do i=1,imt
        buffer_r4_local(i,j,iblock)=buffer_r4_local(i,j,iblock)+viv(i,j,k,iblock)*usmon(i,j,k,iblock)* &
                                    1.0e-6*dzp(k)*dyu(i,j,iblock)/float(imd)
      enddo
      enddo
      enddo
      enddo
!
      call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
      do j=2,jmt_global
      do i=1,imt_global
          buffer_r4_global(i,j)= buffer_r4_global(i,j-1) + buffer_r4_global(i,j)
      end do
      end do
!
      return
      END SUBROUTINE BAROSF
!

!========================================
      SUBROUTINE diag_heat_transport(NNN)
!========================================
!
      implicit none
      integer :: NNN,NIN,iblock
      real(r4) ,dimension(imt,jmt,km,max_blocks_clinic) :: ttpp
      allocate (work1_g(imt_global,jmt_global),work2_g(imt_global,jmt_global),work3_g(imt_global,jmt_global))
!
      ttpp=0.0

! vertical averaged

     do nin = 1, 2
!$OMP PARALLEL DO PRIVATE (K,J,I)
      do iblock= 1, nblocks_clinic
      do k=1,km
      do j=2,jmt-1
      do i=2,imt-1
      if (NNN==1) then
          ttpp(i,j,k,iblock)=vsmon(i,j,k,iblock)/nmonth(mon0)*viv(i,j,k,iblock)*0.25*&
                      (tsmon(i  ,j+1,k,iblock)+tsmon(i  ,j,k,iblock)+&
                       tsmon(i-1,j+1,k,iblock)+tsmon(i-1,j,k,iblock))/nmonth(mon0)*dzp(k)*dxu(i,j,iblock)
      else
          if ( basin(j,j,iblock) == 1 .and. basin(i,j,iblock) == 2) then
          ttpp(i,j,k,iblock)=vsmon(i,j,k,iblock)/nmonth(mon0)*viv(i,j,k,iblock)*0.25*&
                      (ssmon(i  ,j+1,k,iblock)+ssmon(i  ,j,k,iblock)+&
                       ssmon(i-1,j+1,k,iblock)+ssmon(i-1,j,k,iblock))/nmonth(mon0)*dzp(k)*dxu(i,j,iblock)
          end if
      endif
      enddo
      enddo
      enddo
      enddo
!
     work_1 = 0.0_r4
     work_2 = 0.0_r4
     work_3 = 0.0_r4
!
      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic) 
         if (mytid == 0 ) then
         do j=1,jmt_global
         do i=1,imt_global
            work1_g(i,j)=work1_g(i,j)+ buffer_r4_global(i,j)
         enddo
         enddo
         end if
      enddo

      do iblock = 1, nblocks_clinic
      do k =1,km
      do j= 1, jmt
      do i= 1, imt
#ifdef ISO
         ttpp(i,j,k,iblock)=ddymon_iso(i,j,k,nnn,iblock)/nmonth(mon0)*dzp(k)*dxu(i,j,iblock)
#else
         ttpp(i,j,k,iblock)=ddymon(i,j,k,nnn,iblock)/nmonth(mon0)*dzp(k)*dxu(i,j,iblock)
#endif
      end do
      end do
      end do
      end do
!
      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic) 
         if ( mytid ==0 ) then
         do j=1,jmt_global
         do i=1,imt_global
            work2_g(i,j)=work2_g(i,j)+ buffer_r4_global(i,j)
         enddo
         enddo
         end if
      enddo

      do iblock = 1, nblocks_clinic
      do k =1,km
      do j= 1, jmt
      do i= 1, imt
#ifdef ISO
         ttpp(i,j,k,iblock)=aymon_iso(i,j,k,nnn,iblock)/nmonth(mon0)*dzp(k)*dxu(i,j,iblock)
#else
         ttpp(i,j,k,iblock)=aymon(i,j,k,nnn,iblock)/nmonth(mon0)*dzp(k)*dxu(i,j,iblock)
#endif
      end do
      end do
      end do
      end do
!
      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic) 
         if ( mytid ==0 ) then
         do j=1,jmt_global
         do i=1,imt_global
            work3_g(i,j)=work3_g(i,j)+ buffer_r4_global(i,j)
         enddo
         enddo
         end if
      enddo

!
      mth_adv = 0.0_r4
      mth_dif = 0.0_r4
#ifdef ISO
      mth_adv_iso = 0.0_r4
#endif
!
      if (mytid == 0 ) then
      do j=1,jmt_global
         do i=1,imt_global
           mth_adv(j,NIN,NNN)=mth_adv(j,NIN,NNN)-work1_g(i,j)
           mth_dif(j,NIN,NNN)=mth_dif(j,NIN,NNN)+work2_g(i,j)
#ifdef ISO
           mth_adv_iso(j,NIN,NNN)=mth_adv_iso(j,NIN,NNN)-work3_g(i,j)
#endif
          enddo
      enddo

      IF (NNN==1) then
      do j=1,jmt_global
         mth_adv(j,NIN,NNN)=mth_adv(j,NIN,NNN)*D0*CP*1.0E-15
         mth_dif(j,NIN,NNN)=mth_dif(j,NIN,NNN)*D0*CP*1.0E-15
#ifdef ISO
         mth_adv_iso(j,NIN,NNN)=mth_adv_iso(j,NIN,NNN)*D0*CP*1.0E-15
#endif
         mth(j,NIN,NNN)=mth_adv(j,NIN,NNN)+mth_dif(j,NIN,NNN)+mth_adv_iso(j,NIN,NNN)
      enddo
      ELSE
      do j=1,jmt_global
         mth_adv(j,NIN,NNN)=(mth_adv(j,NIN,NNN)*1000.+35)/35.
         mth_dif(j,NIN,NNN)=(mth_dif(j,NIN,NNN)*1000.+35)/35.
#ifdef ISO
         mth_adv_iso(j,NIN,NNN)=(mth_adv_iso(j,NIN,NNN)*1000.+35)/35.
#endif
         mth(j,NIN,NNN)=mth_adv(j,NIN,NNN)+mth_dif(j,NIN,NNN)+mth_adv_iso(j,NIN,NNN)
      enddo
      ENDIF
!
      end if
      ENDDO

      deallocate (work_1,work_2,work_3)
      RETURN
      END SUBROUTINE diag_heat_transport

end module diag_mod
