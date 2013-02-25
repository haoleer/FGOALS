!  CVS: $Id: accumm.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ACCUMM
!     =================
 
#include <def-undef.h>
use param_mod
use dyn_mod
use tracer_mod
use output_mod
use pconst_mod
use forc_mod, only: su,sv,lthf,sshf,lwv,swv,fresh,runoff
use domain
      IMPLICIT NONE

      integer :: iblock

 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   do iblock = 1, nblocks_clinic
      DO J = 1,JMT
         DO I = 1,IMT
            Z0MON (I,J,iblock)= Z0MON (I,J,iblock) + H0 (I,J,iblock)
!
!linpf091126
            sumon (I,J,iblock)= sumon (I,J,iblock) + su(I,J,iblock)        !U windstress
            svmon (I,J,iblock)= svmon (I,J,iblock) + sv(I,J,iblock)        !V windstress
            lthfmon (I,J,iblock)= lthfmon (I,J,iblock) + lthf (I,J,iblock) !latent flux
            sshfmon (I,J,iblock)= sshfmon (I,J,iblock) + sshf (I,J,iblock) !sensible flux
            lwvmon (I,J,iblock)= lwvmon (I,J,iblock) + lwv (I,J,iblock)    !long wave flux
            swvmon (I,J,iblock)= swvmon (I,J,iblock) + swv (I,J,iblock)    !shortwave flux
!linpf091126 
            mldmon (I,J,iblock)= mldmon (I,J,iblock) + amld (I,J,iblock)/100.
!
            akmmon (I,J,K,iblock)= akmmon (I,J,K,iblock) + akmu(I,J,K,iblock)
            aktmon (I,J,K,iblock)= aktmon (I,J,K,iblock) + akt(I,J,K,1,iblock)
            aksmon (I,J,K,iblock)= aksmon (I,J,K,iblock) + akt(I,J,K,2,iblock)
            netmon (I,J,1,iblock)= netmon (I,J,1,iblock) + net (I,J,1,iblock)
            netmon (I,J,2,iblock)= netmon (I,J,2,iblock) + net (I,J,2,iblock)
         END DO
      END DO
   end do
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
   do iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               WSMON (I,J,K,iblock)= WSMON (I,J,K,iblock) + WS (I,J,K,iblock)
               TSMON (I,J,K,iblock)= TSMON (I,J,K,iblock) + AT (I,J,K,1,iblock)
               SSMON (I,J,K,iblock)= SSMON (I,J,K,iblock) + (AT (I,J,K,2,iblock)*1000.+35.)
               USMON (I,J,K,iblock)= USMON (I,J,K,iblock) + U (I,J,K,iblock)
               VSMON (I,J,K,iblock)= VSMON (I,J,K,iblock) - V (I,J,K,iblock)
#if (defined SMAG_OUT)
               AM3MON (I,J,K,iblock)= AM3MON (I,J,K,iblock) + AM3 (I,J,K,iblock)
#endif
               penmon (I,J,K,iblock)= penmon (I,J,K,iblock)+penetrate(i,j,k,iblock)
            END DO
         END DO
      END DO
   end do

!$OMP PARALLEL DO PRIVATE (IBLOCK,N,K,J,I)
   do iblock = 1, nblocks_clinic
      DO N = 1,NTRA
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               trendmon (I,J,K,N,IBLOCK)= trendmon (I,J,K,N,IBLOCK)+trend(i,j,k,n,IBLOCK)
               axmon (I,J,K,N,IBLOCK)= axmon (I,J,K,N,IBLOCK)+ax(i,j,k,n,IBLOCK)
               aymon (I,J,K,N,IBLOCK)= aymon (I,J,K,N,IBLOCK)+ay(i,j,k,n,IBLOCK)
               azmon (I,J,K,N,IBLOCK)= azmon (I,J,K,N,IBLOCK)+az(i,j,k,n,IBLOCK)
               dxmon (I,J,K,N,IBLOCK)= dxmon (I,J,K,N,IBLOCK)+dx(i,j,k,n,IBLOCK)
               dymon (I,J,K,N,IBLOCK)= dymon (I,J,K,N,IBLOCK)+dy(i,j,k,n,IBLOCK)
               dzmon (I,J,K,N,IBLOCK)= dzmon (I,J,K,N,IBLOCK)+dz(i,j,k,n,IBLOCK)
!
               ddymon (I,J,K,N,IBLOCK)= ddymon (I,J,K,N,IBLOCK)+ddy(i,j,k,n,IBLOCK)
#ifdef ISO
!              axmon_iso (I,J,K,N,IBLOCK)= axmon_iso (I,J,K,N,IBLOCK)+ax_iso(i,j,k,n,IBLOCK)
!              aymon_iso (I,J,K,N,IBLOCK)= aymon_iso (I,J,K,N,IBLOCK)+ay_iso(i,j,k,n,IBLOCK)
!              azmon_iso (I,J,K,N,IBLOCK)= azmon_iso (I,J,K,N,IBLOCK)+az_iso(i,j,k,n,IBLOCK)
!              dxmon_iso (I,J,K,N,IBLOCK)= dxmon_iso (I,J,K,N,IBLOCK)+dx_iso(i,j,k,n,IBLOCK)
!              dymon_iso (I,J,K,N,IBLOCK)= dymon_iso (I,J,K,N,IBLOCK)+dy_iso(i,j,k,n,IBLOCK)
!              dzmon_iso (I,J,K,N,IBLOCK)= dzmon_iso (I,J,K,N,IBLOCK)+dz_iso(i,j,k,n,IBLOCK)
!
!              aaymon_iso (I,J,K,N,IBLOCK)= aaymon_iso (I,J,K,N,IBLOCK)+aay_iso(i,j,k,n,IBLOCK)
!              ddymon_iso (I,J,K,N,IBLOCK)= ddymon_iso (I,J,K,N,IBLOCK)+ddy_iso(i,j,k,n,IBLOCK)
#endif
            END DO
         END DO
      END DO
      END DO
   end do
!
!
      RETURN
      END SUBROUTINE ACCUMM
 
 
