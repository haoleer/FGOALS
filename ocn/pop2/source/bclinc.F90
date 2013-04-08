!  CVS: $Id: bclinc.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE BCLINC
!     =================
!     INTEGRATION OF MOMENTUM EQUATION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use forc_mod, only: psa,su,sv
use domain
use grid
use blocks
use smuvh
use operators

      IMPLICIT NONE
      REAL(r8)    :: AA,GGU,WK1,WK2,fil_lat1,fil_lat2
      integer     :: iblock
      REAL(r8)    :: gradx(imt,jmt), grady(imt,jmt)
      type (block) :: this_block

!lhl0711
#ifdef CANUTO
      REAL(r8)    :: AIDIF
      REAL(r8)    :: SBCX(imt,jmt,max_blocks_clinic),BBCX(imt,jmt,max_blocks_clinic)
      REAL(r8)    :: SBCY(imt,jmt,max_blocks_clinic),BBCY(imt,jmt,max_blocks_clinic)


!M
!$OMP PARALLEL DO PRIVATE (J,I) 
      do iblock = 1, nblocks_clinic
         DO J = JSM,JEM
            DO I = 2,IMM
               SBCX(I,J,IBLOCK) = SU (I,J,IBLOCK)* OD0
               SBCY(I,J,IBLOCK) = SV (I,J,IBLOCK)* OD0
               BBCX(I,J,IBLOCK)= C0F*SQRT(UP(I,J,KM,IBLOCK)*UP(I,J,KM,IBLOCK)+VP(I,J,KM,IBLOCK)*VP(I,J,KM,IBLOCK))&
                          *(UP(I,J,KM,IBLOCK)*CAG+SNLAT(J)*VP (I,J,KM,IBLOCK)*SAG)
               BBCY(I,J,IBLOCK)= C0F*SQRT(UP(I,J,KM,IBLOCK)*UP(I,J,KM,IBLOCK)+VP(I,J,KM,IBLOCK)*VP(I,J,KM,IBLOCK))&
                          *(-SNLAT(J)*UP(I,J,KM,IBLOCK)*SAG+VP(I,J,KM,IBLOCK)*CAG)
            ENDDO
         ENDDO
     ENDDO

      AIDIF=0.0  
!      if (ISC/=0)  AIDIF=0.5
 
#endif
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=63.0D0
       fil_lat2=63.0D0
!lhl0711
!---------------------------------------------------------------------
!     ADVECTION + DIFFUSION + CORIOLIS
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) - FCOR(I,J,IBLOCK) * VP (I,J,K,IBLOCK)
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + FCOR(I,J,IBLOCK) * UP (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
 
!---------------------------------------------------------------------
!     PRESSURE GRADIENT FORCES
!---------------------------------------------------------------------
 
!!@@@@ DP'/DX
 
      AA = 0.0
      IF (ISC /= 0) AA = 0.5
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JET
         DO I = 1,IMT
            H0BF (I,J,IBLOCK)= H0BF (I,J,IBLOCK)* ONBB
         END DO
      END DO
  END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JET
         DO I = 1,IMT
            WORK (I,J,IBLOCK)= AA * H0BF (I,J,IBLOCK) + (1.0- AA)* H0BL (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,WKK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JET
         DO I = 1,IMT
            WKK (1)= (PSA (I,J,IBLOCK)* OD0+ WORK (I,J,IBLOCK)* G)* VIT (I,J,1,IBLOCK)
            DO K = 1,KM
               WKK (K +1)= WKK (K) - GG (I,J,K,IBLOCK)* DZP (K)* VIT (I,J,K,IBLOCK)
            END DO
 
            DO K = 1,KM
               WKA (I,J,K,IBLOCK)= 0.25* (WKK (K) + WKK (K +1))
            END DO
         END DO
      END DO
  END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call grad(k, GRADX, GRADY, wka, this_block)
         DO J = JSM,JEM
            DO I = 2,IMM
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) - grady(i,j)
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) - gradx(i,j)
            END DO
         END DO
      END DO
  END DO
 
!!@@@@ G'DH/DX
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (1.0+ OHBT (I,J,IBLOCK)* ZKT (K))* WORK (I,J,IBLOCK)      &
                            * VIT (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,GGU)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         call grad(k, GRADX, GRADY, wka, this_block)
         DO J = JSM,JEM
            DO I = 2,IMM
               GGU = 0.125* (GG (I,J,K,IBLOCK) + GG (I -1,J,K,IBLOCK) + GG (I,J +1,K,IBLOCK) &
                     + GG (I -1,J +1,K,IBLOCK))
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + GGU * grady(i,j)
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) + GGU * gradx(i,j)
            END DO
         END DO
      END DO
   END DO

 
!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------
 
      IF (ISC == 0) THEN
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,WK1,WK2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
                  WK1 = EPEA(I,J,IBLOCK)* DLV(I,J,K,IBLOCK) + EPEB(I,J,IBLOCK)* DLU (I,J,K,IBLOCK)
                  WK2 = EPEA(I,J,IBLOCK)* DLU(I,J,K,IBLOCK) - EPEB(I,J,IBLOCK)* DLV (I,J,K,IBLOCK)
                  DLV (I,J,K,IBLOCK)= WK1* VIV (I,J,K,IBLOCK)
                  DLU (I,J,K,IBLOCK)= WK2* VIV (I,J,K,IBLOCK)
               END DO
            END DO
         END DO
  END DO
      ELSE
!$OMP PARALLEL DO PRIVATE (K,J,I,WK1,WK2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
                  WK1 = EPLA(I,J,IBLOCK)* DLV(I,J,K,IBLOCK) + EPLB(I,J,IBLOCK)*DLU(I,J,K,IBLOCK)
                  WK2 = EPLA(I,J,IBLOCK)* DLU(I,J,K,IBLOCK) - EPLB(I,J,IBLOCK)*DLV(I,J,K,IBLOCK)
                  DLV (I,J,K,IBLOCK)= WK1* VIV(I,J,K,IBLOCK)
                  DLU (I,J,K,IBLOCK)= WK2* VIV(I,J,K,IBLOCK)
               END DO
            END DO
         END DO
  END DO
      END IF
 
 
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
      CALL SMUV (DLU,VIV,KM,fil_lat1)
      CALL SMUV (DLV,VIV,KM,fil_lat1)
!YU 

 
!---------------------------------------------------------------------
!     PREDICTING VC & UC
!---------------------------------------------------------------------
 
      IF (ISC < 1) THEN
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               V (I,J,K,IBLOCK)= VP (I,J,K,IBLOCK) + DLV (I,J,K,IBLOCK)* DTC
               U (I,J,K,IBLOCK)= UP (I,J,K,IBLOCK) + DLU (I,J,K,IBLOCK)* DTC
            END DO
         END DO
      END DO
  END DO
 
!---------------------------------------------------------------------
!@@@  INTERACTION BETWEEN BAROTROPIC AND BAROCLINIC MODES
!---------------------------------------------------------------------
      CALL VINTEG (U,WORK)
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               U (I,J,K,IBLOCK)= (U (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + UB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
!
      CALL VINTEG (V,WORK)
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               V (I,J,K,IBLOCK)= (V (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + VB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO

 
      ISC = ISC +1
    
      ELSE 

 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= VP (I,J,K,IBLOCK) + DLV (I,J,K,IBLOCK)* DTC2
            END DO
         END DO
      END DO
   END DO

 
      CALL VINTEG (WKA,WORK)
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (WKA (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + VB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
 
!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               VP (I,J,K,IBLOCK) = AFC2* V (I,J,K,IBLOCK) + AFC1* (VP (I,J,K,IBLOCK) + WKA (I,J,K,IBLOCK))
               V (I,J,K,IBLOCK) = WKA (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= UP (I,J,K,IBLOCK) + DLU (I,J,K,IBLOCK)* DTC2
            END DO
         END DO
      END DO
  END DO
 
 
      CALL VINTEG (WKA,WORK)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (WKA (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + UB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UP (I,J,K,IBLOCK) = AFC2* U (I,J,K,IBLOCK) + AFC1* (UP (I,J,K,IBLOCK) + WKA (I,J,K,IBLOCK))
               U (I,J,K,IBLOCK) = WKA (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
!YU  Oct. 24, 2005
!lhl0711     IF (MOD(ISC,60)==0) THEN
     IF (MOD(ISC,160)==1) THEN
        CALL SMUV (U,VIV,KM,fil_lat2)
        CALL SMUV (V,VIV,KM,fil_lat2)
        CALL SMUV (UP,VIV,KM,fil_lat2)
        CALL SMUV (VP,VIV,KM,fil_lat2)
     END IF
!YU 
 
      ISC = ISC +1
      END IF

 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) 
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK) + U (I,J,K,IBLOCK)
               VTF (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK) + V (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
 
!
      RETURN
      END SUBROUTINE BCLINC
