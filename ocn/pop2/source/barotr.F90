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
      IMPLICIT NONE

      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2
      integer :: iblock
      real(r8):: hduk(imt,jmt) , hdvk(imt,jmt), gradx(imt,jmt),grady(imt,jmt), div_out(imt,jmt)
      type(block):: this_block
!
!---------------------------------------------------------------------
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
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,11,IBLOCK)= UBP (I,J,IBLOCK)
               WKA (I,J,12,IBLOCK)= VBP (I,J,IBLOCK)
            END DO
         END DO
      END DO
         CALL SMAG2 (1)
#if (defined SMAG_FZ)
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,5,IBLOCK)= VIV (I,J,1,IBLOCK)* (0.5* OUX (J)* (WKA (I +1,J,7,IBLOCK) &
                            - WKA (I -1,J,7,IBLOCK)) &
               - R2E (J)* WKA (I,J +1,8,IBLOCK) + R2F (J)* WKA (I,J -1,8,IBLOCK))
               WKA (I,J,6,IBLOCK)= VIV (I,J,1,IBLOCK)* (0.5* OUX (J)* (WKA (I +1,J,9,IBLOCK) &
                            - WKA (I -1,J,9,IBLOCK)) &
               - R3E (J)* WKA (I,J +1,10,IBLOCK) + R3F (J)* WKA (I,J -1,10,IBLOCK)    &
                            + R4E (J)* WKA (I,J,7,IBLOCK))
            END DO
         END DO
     END DO

!-new
#else
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = JSM,JEM
            DO I = 2,IMM
               WKA (I,J,5,IBLOCK)= VIV (I,J,1,IBLOCK)* (0.5* OUX (J)* (WKA (I +1,J,7,IBLOCK) &
                            - WKA (I -1,J,7,IBLOCK)) &
               - R2E (J)* WKA (I,J +1,8,IBLOCK) + R2F (J)* WKA (I,J -1,8,IBLOCK))
               WKA (I,J,6)= VIV (I,J,1,IBLOCK)* (0.5* OUX (J)* (WKA (I +1,J,9,IBLOCK) &
                            - WKA (I -1,J,9,IBLOCK)) &
               - R3E (J)* WKA (I,J +1,10,IBLOCK) + R3F (J)* WKA (I,J -1,10,IBLOCK)    &
                            + R4E (J)* WKA (I,J,7,IBLOCK))
            END DO
         END DO
     END DO
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

!---------------------------------------------------------------------
!     + (g'-1)g*dH/dr
!---------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (IBLOCK)
     DO IBLOCK = 1, NBLOCKS_CLINIC
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
            DO J = JSM,JEM
               DO I = 2,IMM
                  WKA (I,J,3,IBLOCK)= EBEA(I,J,IBLOCK)* WKA (I,J,1,IBLOCK) - EBEB(I,J,IBLOCK)* WKA (I,J,2,IBLOCK)
                  WKA (I,J,4,IBLOCK)= EBEA(I,J,IBLOCK)* WKA (I,J,2,IBLOCK) + EBEB(I,J,IBLOCK)* WKA (I,J,1,IBLOCK)
               END DO
            END DO
     END DO
         ELSE
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = JSM,JEM
               DO I = 2,IMM
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
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,1,IBLOCK)= UB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
               WKA (I,J,2,IBLOCK)= VB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
            END DO
         END DO
     END DO


!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) 
    DO IBLOCK = 1, NBLOCKS_CLINIC
         call div(1,DIV_OUT,wka(:,:,1,iblock),wka(:,:,2,iblock),this_block)
         DO J = 2, jmt-1
            DO I = 2,imt-1
               WORK (I,J,IBLOCK)=VIT(I,J,1,IBLOCK)*(-1)*div_out(i,j)
             END DO
          ENDDO
    END DO
!
!---------------------------------------------------------------------
!     PREDICTING VB , UB & H0
!---------------------------------------------------------------------
!YU  Oct. 24,2005
         CALL SMUV (WKA(:,:,3,:) ,VIV,1,fil_lat1)
         CALL SMUV (WKA(:,:,4,:) ,VIV,1,fil_lat1)
         CALL SMZ0 (WORK,VIT,fil_lat1)
!YU  Oct. 24,2005
!
         IF (ISB < 1) THEN

!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = JSM,JEM
            DO I = 1,IMT
               UB (I,J,IBLOCK)= UBP (I,J,IBLOCK) + WKA (I,J,3,IBLOCK)* DTB
               VB (I,J,IBLOCK)= VBP (I,J,IBLOCK) + WKA (I,J,4,IBLOCK)* DTB
               H0 (I,J,IBLOCK)= H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK) * DTB
            END DO
         END DO
     END DO

!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------

            CALL SMUV (UB ,VIV,1,fil_lat1)
            CALL SMUV (VB ,VIV,1,fil_lat1)
            CALL SMZ0 (H0,VIT,fil_lat1)
!           call FILTER_TRACER(h0,vit_1d,FIL_LAT1,1)

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
         DO J = JSM,JEM
            DO I = 1,IMT
               WKA (I,J,1,IBLOCK) = UBP (I,J,IBLOCK) + WKA (I,J,3,IBLOCK)* DTB2
               WKA (I,J,2,IBLOCK) = VBP (I,J,IBLOCK) + WKA (I,J,4,IBLOCK)* DTB2
               WORK(I,J,IBLOCK)   = H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK)* DTB2
            END DO
         END DO
    END DO

!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------


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
            CALL SMUV (UB ,VIV,1,fil_lat2)
            CALL SMUV (VB ,VIV,1,fil_lat2)
            CALL SMZ0 (H0 ,VIT,fil_lat2)
            CALL SMUV (UBP,VIV,1,fil_lat2)
            CALL SMUV (VBP,VIV,1,fil_lat2)
            CALL SMZ0 (H0P,VIT,fil_lat2)
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

      deallocate(dlub,dlvb)
      RETURN
      END SUBROUTINE BAROTR


