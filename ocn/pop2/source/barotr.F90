      SUBROUTINE BAROTR
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use msg_mod
use domain
use grid
use blocks
use hmix_del2
use hmix_del4
use operators
use smuvh
use POP_GridHorzMod
use POP_HaloMod
use global_reductions
use gather_scatter
use distribution
use constant_mod
      IMPLICIT NONE

      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2, ek0, maxz0, minz0
!     real(r8)  ::  ttt(imt_global,jmt_global)
      integer :: iblock,ii1,jj1,ii2,jj2
      real(r8):: hduk(imt,jmt) , hdvk(imt,jmt), gradx(imt,jmt),grady(imt,jmt), div_out(imt,jmt)
      type(block):: this_block

!      Define the threthold latitute for zonal smoother
       fil_lat1=65.0D0
       fil_lat2=65.0D0
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
      wka=0
      work=0
!---------------------------------------------------------------------
!     EULER BACKWARD SCHEME IS USED FOR THE FIRST STEP OF EVERY MONTH
!     IEB=0: LEAP-FROG SCHEME; IEB=1: EULER BACKWARD SCHEME
!---------------------------------------------------------------------
      IEB = 0 ; IEB_LOOP=0

      IF (ISB == 0)  THEN
         IEB = 1 ; IEB_LOOP=1
      END IF
!
      baro_loop : DO NC = 1,NBB+IEB_LOOP




      if (IEB==1.or.ISB>1) then

!---------------------------------------------------------------------
!     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
!---------------------------------------------------------------------

#if ( defined SMAG1)
#else
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (IBLOCK)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call hdiffu_del4(1, HDUK, HDVK, ubp(:,:,iblock), vbp(:,:,iblock), this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,5,IBLOCK)= hduk(i,j)
               WKA (I,J,6,IBLOCK)= hdvk(i,j)
            END DO
         END DO
     END DO
!
!!!!!!
#else
!$OMP PARALLEL DO PRIVATE (IBLOCK)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call hdiffu_del2(1, HDUK, HDVK, ubp(:,:,iblock), vbp(:,:,iblock), this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,5,IBLOCK)= hduk(i,j)
               WKA (I,J,6,IBLOCK)= hdvk(i,j)
            END DO
         END DO
    END DO

#endif
#endif
!
!       call gather_global(ttt, wka(:,:,5,:), master_task,distrb_clinic)
!       if (mytid == 0) then
!          write(221,*) ((ttt(i,j),i=1,imt_global),j=1,jmt_global)
!          close(221)
!       end if
!
            IF (mod(isb,36)  == 1 ) THEN
!$OMP PARALLEL DO PRIVATE (J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 3, jmt-2
               DO I = 3, imt-2
                  DLUB (I,J,IBLOCK)= DLUB (I,J,IBLOCK) + WKA (I,J,5,IBLOCK)
                  DLVB (I,J,IBLOCK)= DLVB (I,J,IBLOCK) + WKA (I,J,6,IBLOCK)
               END DO
            END DO
     END DO
            END IF
         END IF

!---------------------------------------------------------------------
!     + (g'-1)g*dH/dr
!---------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (IBLOCK)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call grad(1, GRADX, GRADY, H0(:,:,iblock), this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               gstar=(WGP (I,J,IBLOCK) -1.0)*G 
               WKA (I,J,1,IBLOCK) = WKA (I,J,5,IBLOCK) + gstar*GRADX(I,J)
               WKA (I,J,2,IBLOCK) = WKA (I,J,6,IBLOCK) + gstar*GRADY(I,J)
            END DO
         END DO
     END DO

!Yu
!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
        call tgrid_to_ugrid(work(:,:,iblock),h0(:,:,iblock),iblock)
   end do
!
!---------------------------------------------------------------------
!     COMPUTING DU & DV
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,1,IBLOCK)= VIV (I,J,1,IBLOCK)* ( WKA (I,J,1,IBLOCK) + DLUB (I,J,IBLOCK)     &
                              - FCOR(I,J,Iblock)* VBP (I,J,IBLOCK) + &
               PAX (I,J,IBLOCK) + PXB (I,J,IBLOCK) - WORK (I,J,IBLOCK)* WHX (I,J,IBLOCK) )
               WKA (I,J,2,IBLOCK)= VIV (I,J,1,IBLOCK)* ( WKA (I,J,2,IBLOCK) + DLVB (I,J,IBLOCK)     &
                              + FCOR(I,J,Iblock)* UBP (I,J,IBLOCK) + &
               PAY (I,J,IBLOCK) + PYB (I,J,IBLOCK) - WORK (I,J,IBLOCK)* WHY (I,J,IBLOCK) )
            END DO
         END DO
     END DO

!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------

         IF (ISB == 0) THEN
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 3, jmt-2
               DO I = 3, imt-2
                  WKA (I,J,3,IBLOCK)= EBEA(I,J,IBLOCK)* WKA (I,J,1,IBLOCK) - EBEB(I,J,IBLOCK)* WKA (I,J,2,IBLOCK)
                  WKA (I,J,4,IBLOCK)= EBEA(I,J,IBLOCK)* WKA (I,J,2,IBLOCK) + EBEB(I,J,IBLOCK)* WKA (I,J,1,IBLOCK)
               END DO
            END DO
     END DO
         ELSE
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 3,jmt-2
               DO I = 3, imt-2
                  WKA (I,J,3,IBLOCK)= EBLA (I,J,Iblock)* WKA (I,J,1,IBLOCK) - EBLB (I,J,Iblock)* WKA (I,J,2,IBLOCK)
                  WKA (I,J,4,IBLOCK)= EBLA (I,J,Iblock)* WKA (I,J,2,IBLOCK) + EBLB (I,J,Iblock)* WKA (I,J,1,IBLOCK)
               END DO
            END DO
     END DO
         END IF



!---------------------------------------------------------------------
!     COMPUTING DH0
!---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 2, JMT
            DO I = 1,IMT-1
               WKA (I,J,1,IBLOCK)= UB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
               WKA (I,J,2,IBLOCK)= VB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
            END DO
         END DO
     END DO


!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) 
    DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call div(1,DIV_OUT,wka(:,:,1,iblock),wka(:,:,2,iblock),this_block)
         DO J = 2, jmt-2
            DO I = 3,imt-2
               WORK (I,J,IBLOCK)=VIT(I,J,1,IBLOCK)*(-1)*div_out(i,j)
             END DO
          ENDDO
    END DO


!
!---------------------------------------------------------------------
!     PREDICTING VB , UB & H0
!---------------------------------------------------------------------
         call POP_HaloUpdate(work , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
         call POP_HaloUpdate(wka(:,:,3,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
         call POP_HaloUpdate(wka(:,:,4,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
 

!YU  Oct. 24,2005
      if (trim(horiz_grid_opt) == 'lat_lon') then
         CALL SMUV_2D (WKA(:,:,3,:) ,VIV(:,:,1,:),fil_lat1)
         CALL SMUV_2D (WKA(:,:,4,:) ,VIV(:,:,1,:),fil_lat1)
         CALL SMZ0 (WORK,VIT(:,:,1,:),fil_lat1)
      end if
!YU  Oct. 24,2005
!   DO IBLOCK = 1, NBLOCKS_CLINIC
!        DO J = 1, jmt
!           DO I = 1,imt
!     if (mytid ==0 .and. i > 3 .and. i < 20 .and. j >5 .and. j < 8 ) then
!     write(200+isb,*) i,j, wka(i,j,3,1),wka(i,j,4,1),work(i,j,1)
!     end if
!           END DO
!        END DO
!    END DO
!

!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1, jmt
            DO I = 1,imt
               UB (I,J,IBLOCK)= UBP (I,J,IBLOCK) + WKA (I,J,3,IBLOCK)* DTB
               VB (I,J,IBLOCK)= VBP (I,J,IBLOCK) + WKA (I,J,4,IBLOCK)* DTB
               H0 (I,J,IBLOCK)= H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK) * DTB
            END DO
         END DO
     END DO



!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------

         if (trim(horiz_grid_opt) == 'lat_lon') then
            CALL SMUV_2D (UB ,VIV(:,:,1,:),fil_lat1)
            CALL SMUV_2D (VB ,VIV(:,:,1,:),fil_lat1)
            CALL SMZ0 (H0,VIT(:,:,1,:),fil_lat1)
         end if
!        

         IF (IEB == 0) THEN
            ISB = ISB +1
            ubp=ub
            vbp=vb
            h0p=h0
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
               DO I = 1,IMT
                  H0F (I,J,IBLOCK) = H0F (I,J,IBLOCK) + H0 (I,J,IBLOCK)
                  H0BF (I,J,IBLOCK) = H0BF (I,J,IBLOCK) + H0 (I,J,IBLOCK)
               END DO
            END DO
    END DO
            ieb = 1
            cycle baro_loop
         END IF

         IEB = 0


      END DO baro_loop

      deallocate(dlub,dlvb)
  call mpi_barrier(mpi_comm_ocn,ierr)

!     maxz0=0.0D0
!     minz0=0.0D0
!     do iblock =1 ,nblocks_clinic
!     do j = 3, jmt-2
!     do i = 3, imt-2
!        if (maxz0 < h0(i,j,iblock) ) then
!           maxz0 = h0(i,j,iblock)
!           ii1=i
!           jj1=j
!        end if
!        if (minz0 > h0(i,j,iblock) ) then
!           minz0 = h0(i,j,iblock)
!           ii2=i
!           jj2=j
!        end if
!     end do
!     end do
!     end do
!     write(130+mytid,*) mytid,ii1,jj1,maxz0
!     write(130+mytid,*) mytid,ii2,jj2,minz0
!     write(130+mytid,*) (h0(i,1,1),i=1,imt)
!     write(130+mytid,*) (h0(i,2,1),i=1,imt)
!     write(130+mytid,*) (h0(i,3,1),i=1,imt)
!     write(130+mytid,*) (h0(i,4,1),i=1,imt)
!     write(130+mytid,*) (h0(i,5,1),i=1,imt)
      RETURN
      END SUBROUTINE BAROTR


