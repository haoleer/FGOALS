!  CVS: $Id: intfor.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE INTFOR
!     =================
!     INTERPOLATES MONTH MEAN FILEDS
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use forc_mod
use msg_mod
use domain
      IMPLICIT NONE
 
      INTEGER :: IPT1,IPT2, IBLOCK
      REAL(r8):: FACTOR
      INTEGER :: IYEAR,IREC
      CHARACTER (LEN=180) :: FNAME 
      CHARACTER (LEN=4)   :: FHEAD 
      CHARACTER (LEN=3)   :: FTAIL(12)
      REAL(r8),dimension(:,:,:),allocatable :: sx,sy
      REAL(r8) :: epsln
      allocate(sx(imt,jmt,max_blocks_clinic))
      allocate(sy(imt,jmt,max_blocks_clinic))
!
      epsln = 1.0D-25
!
      FTAIL=RESHAPE((/'jan','feb','mar','apr','may','jun', &
                   'jul','aug','sep','oct','nov','dec'/),(/12/)) 
!
! 
      IF ( IDAY <= 15) THEN
         IPT1 = MON0-1
         IF (IPT1 == 0 ) IPT1 = 12
         IPT2 = MON0
         FACTOR = FLOAT (IDAY -15)/ FLOAT (NMONTH (IPT1)) + 1
      ELSE
         IPT1 = MON0
         IPT2 = MOD (MON0,12) +1
         FACTOR = FLOAT (IDAY -15)/ FLOAT (NMONTH (IPT1))
      END IF
 
 
!lhl0711
            PSA3=1000.
!lhl0711

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JET
         DO I = 1,IMT
            SU (I,J,IBLOCK)= ( (SU3 (I,J,IPT2,IBLOCK) - SU3 (I,J,IPT1,IBLOCK))               &
                             * FACTOR + SU3 (I,J,IPT1,IBLOCK))
            SV (I,J,IBLOCK)= ( (SV3 (I,J,IPT2,IBLOCK) - SV3 (I,J,IPT1,IBLOCK))               &
                             * FACTOR + SV3 (I,J,IPT1,IBLOCK))
            PSA (I,J,IBLOCK)= (PSA3 (I,J,IPT2,IBLOCK) - PSA3 (I,J,IPT1,IBLOCK))              &
                             * FACTOR + PSA3 (I,J,IPT1,IBLOCK)
            TSA (I,J,IBLOCK)= (TSA3 (I,J,IPT2,IBLOCK) - TSA3 (I,J,IPT1,IBLOCK))              &
                             * FACTOR + TSA3 (I,J,IPT1,IBLOCK)
            SSS (I,J,IBLOCK)= (SSS3 (I,J,IPT2,IBLOCK) - SSS3 (I,J,IPT1,IBLOCK))              &
                             * FACTOR + SSS3 (I,J,IPT1,IBLOCK)
            SWV (I,J,IBLOCK)= ( (SWV3 (I,J,IPT2,IBLOCK) - SWV3 (I,J,IPT1,IBLOCK))            &
                             * FACTOR + SWV3 (I,J,IPT1,IBLOCK))
            UVA (I,J,IBLOCK)= ( (UVA3 (I,J,IPT2,IBLOCK) - UVA3 (I,J,IPT1,IBLOCK))            &
                             * FACTOR + UVA3 (I,J,IPT1,IBLOCK))
            QAR (I,J,IBLOCK)= ( (QAR3 (I,J,IPT2,IBLOCK) - QAR3 (I,J,IPT1,IBLOCK))            &
                             * FACTOR + QAR3 (I,J,IPT1,IBLOCK))
            CLD (I,J,IBLOCK)= ( (CLD3 (I,J,IPT2,IBLOCK) - CLD3 (I,J,IPT1,IBLOCK))            &
                             * FACTOR + CLD3 (I,J,IPT1,IBLOCK))
            SST (I,J,IBLOCK)= ( (SST3 (I,J,IPT2,IBLOCK) - SST3 (I,J,IPT1,IBLOCK))            &
                             * FACTOR + SST3 (I,J,IPT1,IBLOCK))
!lhl
           NSWV (I,J,IBLOCK)= ( (NSWV3(I,J,IPT2,IBLOCK) - NSWV3(I,J,IPT1,IBLOCK))            &
                             * FACTOR + NSWV3 (I,J,IPT1,IBLOCK))
           DQDT (I,J,IBLOCK)= ( (DQDT3(I,J,IPT2,IBLOCK) - DQDT3(I,J,IPT1,IBLOCK))            &
                             * FACTOR + DQDT3(I,J,IPT1,IBLOCK))
!lhl
         END DO
      END DO
   END DO
!lhl1204
!
! Calculate the friction velocity at T-grid
!
!$OMP PARALLEL DO PRIVATE (IBLOCK)
      DO IBLOCK = 1, NBLOCKS_CLINIC
          call ugrid_to_tgrid(sx(:,:,iblock),su(:,:,iblock),iblock,1)
          call ugrid_to_tgrid(sy(:,:,iblock),sv(:,:,iblock),iblock,1)
      END DO

!M
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JET
         DO I = 1,IMT
            USTAR(I,J,iblock)=sqrt(sqrt(sx(I,J,iblock)*sx(I,J,iblock)+ &
                                        sy(I,J,iblock)*sy(I,J,iblock))*OD0)*vit(i,j,1,iblock)
         END DO
      END DO
      END DO

#ifdef DEBUG
      call chk_var2d(sx,"sx",1)
      call chk_var2d(sy,"sy",1)
      call chk_var2d(ustar,"us",1)
#endif
!
!!lhl1204
   

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J = JST,JET
         DO I = 1,IMT 
           seaice(I,J,iblock)= ( (seaice3 (I,J,IPT2,IBLOCK) - seaice3 (I,J,IPT1,IBLOCK))     &
                    * FACTOR + seaice3 (I,J,IPT1,IBLOCK))
           END DO
        END DO
      END DO

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J = JST,JET
         DO I = 1,IMT
           runoff(I,J,iblock)= runoff3(i,j,1,iblock)
           END DO
        END DO
      END DO

#if (defined SOLARCHLORO) 
!M
!$OMP PARALLEL DO PRIVATE (J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J = JST,JET
         DO I = 1,IMT
           chloro(I,J,iblock)= ( (chloro3 (I,J,IPT2,IBLOCK) - chloro3 (I,J,IPT1,IBLOCK))     &
                    * FACTOR + chloro3 (I,J,IPT1,IBLOCK))
           END DO
        END DO
      END DO
#endif
 
      deallocate(sx)
      deallocate(sy)
      RETURN
      END SUBROUTINE INTFOR
 
 
