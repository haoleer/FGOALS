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
use distribution
use constant_mod
      IMPLICIT NONE

      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2
      integer :: iblock
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
      call energy

      if (IEB==1.or.ISB>1) then

!---------------------------------------------------------------------
!     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
!---------------------------------------------------------------------

#if ( defined SMAG1)
#if (defined SMAG_FZ)
#else
#endif
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

      write(120+mytid,*) ((wka(i,j,5,1),i=3,imt-2),j=3,jmt-2)
      close(120+mytid)
      write(140+mytid,*) ((wka(i,j,6,1),i=3,imt-2),j=3,jmt-2)
      close(140+mytid)
      write(160+mytid,*) ((dlub(i,j,1),i=3,imt-2),j=3,jmt-2)
      close(160+mytid)
      stop
!---------------------------------------------------------------------
!     + (g'-1)g*dH/dr
!---------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (IBLOCK)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call grad(1, GRADX, GRADY, H0, this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               gstar=(WGP (I,J,IBLOCK) -1.0)*G *0.5
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
!   write(120+mytid,*) isb, global_maxval(work,distrb_clinic,field_loc_center ), global_minval(work,distrb_clinic,field_loc_center)
!   write(120+mytid,*) global_maxval(wka(:,:,3,:),distrb_clinic,field_loc_center ), global_minval(wka(:,:,3,:),distrb_clinic,field_loc_center)
!   write(120+mytid,*) global_maxval(wka(:,:,4,:),distrb_clinic,field_loc_center ), global_minval(wka(:,:,4,:),distrb_clinic,field_loc_center)
!   write(120+mytid,*) global_maxval(wka(:,:,1,:),distrb_clinic,field_loc_center ), global_minval(wka(:,:,1,:),distrb_clinic,field_loc_center)
!   write(120+mytid,*) global_maxval(wka(:,:,2,:),distrb_clinic,field_loc_center ), global_minval(wka(:,:,2,:),distrb_clinic,field_loc_center)


       if (mytid == 1) then
           write(112,*)"OK-----------1"
           close(112)
       end if
!
!---------------------------------------------------------------------
!     PREDICTING VB , UB & H0
!---------------------------------------------------------------------
         call POP_HaloUpdate(work , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
       if (mytid == 1) then
           write(112,*)"OK-----------2"
           close(112)
       end if
         call POP_HaloUpdate(wka(:,:,3,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
       if (mytid == 1) then
           write(112,*)"OK-----------3"
           close(112)
       end if
         call POP_HaloUpdate(wka(:,:,4,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
 
       if (mytid == 1) then
           write(112,*)"OK-----------4"
           close(112)
       end if

!YU  Oct. 24,2005
         CALL SMUV_2D (WKA(:,:,3,:) ,VIV(:,:,1,:),fil_lat1)
       if (mytid == 1) then
           write(112,*)"OK-----------5"
           close(112)
       end if
         CALL SMUV_2D (WKA(:,:,4,:) ,VIV(:,:,1,:),fil_lat1)
       if (mytid == 1) then
           write(112,*)"OK-----------6"
           close(112)
       end if
         CALL SMZ0 (WORK,VIT(:,:,1,:),fil_lat1)
       if (mytid == 1) then
           write(112,*)"OK-----------7"
           close(112)
       end if
!YU  Oct. 24,2005
!
         IF (ISB < 1) THEN

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

            CALL SMUV_2D (UB ,VIV(:,:,1,:),fil_lat1)
            CALL SMUV_2D (VB ,VIV(:,:,1,:),fil_lat1)
            CALL SMZ0 (H0,VIT(:,:,1,:),fil_lat1)
!        

         IF (IEB == 0) THEN
            ISB = ISB +1
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
               DO I = 1,IMT
                  H0F (I,J,IBLOCK) = H0F (I,J,IBLOCK) + H0 (I,J,IBLOCK)
                  H0BF (I,J,IBLOCK) = H0BF (I,J,IBLOCK) + H0 (I,J,IBLOCK)
               END DO
            END DO
    END DO
            cycle baro_loop
         END IF

         IEB = 0

         cycle baro_loop

      ELSE


!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 3,jmt-2
            DO I = 3,imt-2
               WKA (I,J,1,IBLOCK) = UBP (I,J,IBLOCK) + WKA (I,J,3,IBLOCK)* DTB2
               WKA (I,J,2,IBLOCK) = VBP (I,J,IBLOCK) + WKA (I,J,4,IBLOCK)* DTB2
               WORK(I,J,IBLOCK)   = H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK)* DTB2
            END DO
         END DO
    END DO

!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
!        

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = JST,JET
            DO I = 1,IMT
               UBP (I,J,IBLOCK) = AFB2* UB (I,J,IBLOCK) + AFB1* (UBP (I,J,IBLOCK) + WKA (I,J,1,IBLOCK))
               UB (I,J,IBLOCK) = WKA (I,J,1,IBLOCK)*VIV(I,J,1,IBLOCK)
               VBP (I,J,IBLOCK) = AFB2* VB (I,J,IBLOCK) + AFB1* (VBP (I,J,IBLOCK) + WKA (I,J,2,IBLOCK))
               VB (I,J,IBLOCK) = WKA (I,J,2,IBLOCK)*VIV(I,J,1,IBLOCK)
               H0P (I,J,IBLOCK) = AFB2* H0 (I,J,IBLOCK) + AFB1* (H0P (I,J,IBLOCK) + WORK(I,J,IBLOCK))
               H0 (I,J,IBLOCK) = WORK (I,J,IBLOCK)
            END DO
         END DO
    END DO

!YU  Oct. 24,2005
!lhl0711         IF (MOD(ISB,1200)==0) THEN
         IF (MOD(ISB,1440)==1) THEN
            CALL SMUV_2D (UB ,VIV(:,:,1,:),fil_lat2)
            CALL SMUV_2D (VB ,VIV(:,:,1,:),fil_lat2)
            CALL SMZ0 (H0 ,VIT(:,:,1,:),fil_lat2)
            CALL SMUV_2D (UBP,VIV(:,:,1,:),fil_lat2)
            CALL SMUV_2D (VBP,VIV(:,:,1,:),fil_lat2)
            CALL SMZ0 (H0P,VIT(:,:,1,:),fil_lat2)
         END IF

!YU  Oct. 24,2005

         ISB = ISB +1
      END IF



!$OMP PARALLEL DO PRIVATE (J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
            DO I = 1,IMT
               H0F (I,J,IBLOCK) = H0F (I,J,IBLOCK) + H0 (I,J,IBLOCK)
               H0BF (I,J,IBLOCK) = H0BF (I,J,IBLOCK) + H0 (I,J,IBLOCK)
            END DO
         END DO
    END DO

      END DO baro_loop
      close(120+mytid)

      deallocate(dlub,dlvb)
  call mpi_barrier(mpi_comm_ocn,ierr)

      RETURN
      END SUBROUTINE BAROTR


