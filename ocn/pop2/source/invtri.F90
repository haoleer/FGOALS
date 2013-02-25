!  CVS: $Id: invtri.F90,v 1.5 2003/08/12 09:06:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE INVTRI (WK,TOPBC,DCB,AIDIF,C2DTTS)
!     =================
use param_mod
use pconst_mod
use domain
      IMPLICIT NONE
      REAL    :: WK (IMT,JMT,KM,max_blocks_clinic),TOPBC (IMT,JMT,max_blocks_clinic),DCB (IMT,JMT,KM,max_blocks_clinic)
      REAL    :: A8 (KM),B8 (KM),C8 (KM),D8 (KM),E8 (0:KM),F8 (0:KM)
      REAL    :: AIDIF,C2DTTS,G0
      INTEGER :: KZ, IBLOCK
 
!$OMP PARALLEL DO PRIVATE ( iblock, j, i)
  do iblock = 1, nblocks_clinic
      JJJ: DO J = 2,JMM
         III: DO I = 2,IMM
            IF (ITNU (I,J,iblock) == 0) cycle III
            DO K = 2,KM  !lyc
               A8 (K) = DCB (I,J,K -1,iblock)* ODZT (K )* ODZP (K)* C2DTTS * AIDIF
               D8 (K) = WK (I,J,K,iblock)
            ENDDO
            DO K=2,KM-1
               C8 (K) = DCB (I,J,K,iblock)* ODZT (K+1 )* ODZP (K)* C2DTTS * AIDIF
               B8 (K) = 1.0+ A8 (K) + C8 (K)
               E8 (K -1) = 0.0
               F8 (K -1) = 0.0
            END DO
 
!     B. C. AT TOP
            K = 1
            A8 (K) = ODZP (K)* C2DTTS * AIDIF
            C8 (K) = DCB (I,J,K,iblock)* ODZT (K+1)* ODZP (K)* C2DTTS * AIDIF
            B8 (K) = 1.0+ C8 (K)
            D8 (K) = WK (I,J,K,iblock)
            E8 (K -1) = 0.0
            F8 (K -1) = 0.0
!     B. C. AT BOTTOM
            KZ = ITNU (I,J,iblock)
            IF (KZ /= 0) THEN
               B8 (KZ) = 1.0+ A8 (KZ)
               C8 (KZ) = ODZP (KZ)* C2DTTS * AIDIF
               E8 (KZ) = 0.0
               F8 (KZ) = 0.0
            END IF
 
!     NOW INVERT
            DO K = KM,1, -1
               IF (K <= ITNU (I,J,iblock)) THEN
                  G0 = 1.0/ (B8 (K) - C8 (K)* E8 (K))
                  E8 (K -1) = A8 (K)* G0
                  F8 (K -1) = (D8 (K) + C8 (K)* F8 (K))* G0
               END IF
 
            END DO
 
!     B.C. AT SURFACE
            WK (I,J,1,iblock) = (E8 (0)* TOPBC (I,J,iblock) + F8 (0))* VIT (I,J,1,iblock)
            DO K = 2,KM
               WK (I,J,K,iblock)= (E8 (K -1)* WK (I,J,K -1,iblock) + F8 (K -1))* VIT (I,J,K,iblock)
            END DO
 
         END DO III
      END DO JJJ
  end do
!
      RETURN
      END SUBROUTINE INVTRI
 
 
#else
      SUBROUTINE INVTRI ()
      RETURN
      END SUBROUTINE INVTRI
#endif 
