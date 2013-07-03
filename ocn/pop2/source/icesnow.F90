!  CVS: $Id: icesnow.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ==================
      SUBROUTINE ICESNOW
!     ==================
!     Sea Ice Model
#include <def-undef.h> 
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use domain
use grid, only :  kmt
      IMPLICIT NONE
!
      real(r8) :: sal_ocn, sal_ice, tdiff, heat_ice_fusion
      real(r8) ::  t_mix, s_mix
      integer :: iblock
!
      heat_ice_fusion = 3.337d+5        !J/Kg
      sal_ocn=35.0D0                     !Reference salinity for ocean
      sal_ice= 0.0D0                     !Reference salinity for sea ice
!
!------------------------------------------------------------
!  if SST exceed -1.8C to restore it to -1.8C
!------------------------------------------------------------
 
      DO IBLOCK = 1, NBLOCKS_CLINIC
 
      KKK : DO K = 1, 1
      JJJ : DO J = 1, JMT
         III : DO I = 1,IMT
 
            IF (KMT(I,J,IBLOCK) == 0) CYCLE III
 
            IF (AT (I,J,K,1,IBLOCK) < TBICE) THEN
#ifdef  COUP
               tdiff        = TBICE- AT(i,j,k,1,IBLOCK)
               licomqice (I,J  ,IBLOCK) = licomqice(i,j,IBLOCK)+tdiff*dzp(k)/dzp(1)
               at(i,j,1,2,IBLOCK)=at(i,j,1,2,IBLOCK)+tdiff*(sal_ocn-sal_ice)  &
                          *CP/heat_ice_fusion*0.001D0
#endif
               AT (I,J,k,1,iblock) = TBICE
            END IF
 
         END DO  III
      END DO JJJ
      END DO KKK
!
#ifdef  COUP
      DO J=1,JMT
      DO I=1,IMT
         IF (licomqice(I,J,IBLOCK) > 0.0 .and. at(i,j,1,1,IBLOCK) > TBICE) THEN
              tdiff=min((at(i,j,1,1,iblock)-tbice),licomqice(i,j,iblock))
              licomqice(i,j,IBLOCK)=licomqice(i,j,IBLOCK)-tdiff
              at(i,j,1,1,IBLOCK)=at(i,j,1,1,IBLOCK)-tdiff
              at(i,j,1,2,IBLOCK)=at(i,j,1,2,IBLOCK)-tdiff*(sal_ocn-sal_ice)   &
                         *CP/heat_ice_fusion*0.001D0
          END IF
      END DO
      END DO 
#endif
!
      DO K=2,KM
      DO J=1,JMT
      DO I=1,IMT
         IF (AT(I,J,K,1,IBLOCK) < TBICE) AT(I,J,K,1,IBLOCK)=TBICE
      ENDDO
      ENDDO
      ENDDO
!
   END DO
        if (mytid == 3 .and. isc < 8) then
           i=83
           j=24
           write(168,*) "ISC=", ISC
           write(168,*) at(i,j,1,1,1),at(i,j,2,1,1)
           write(168,*) at(i,j,1,2,1),at(i,j,2,2,1), licomqice(i,j,1)
           if(isc ==7) close(168)
        end if

      RETURN
      END SUBROUTINE ICESNOW
 
 
