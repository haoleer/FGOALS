!     =================
      SUBROUTINE READYC
!     =================
!     PREPARATION OF BAROTROPIC AND BAROCLINIC INTEGRATION
 
!     ADVECTION + DIFFUSION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use tracer_mod
use pmix_mod
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL
use domain
use grid
use blocks
use advection
use operators
use hmix_del2
use hmix_del4
use msg_mod 
      IMPLICIT NONE
!      REAL(r8)  :: WKP (KMP1)
      INTEGER   :: IWK,n2, iblock
      REAL(r8)  :: WK1 (KM) ,WK2 (KM), WK3 (KM)
      REAL(r8)  :: WP1 (KM) ,WP2 (KM), WP3 (KM)
      REAL(r8)  :: WP4 (KM) ,WP5 (KM), WP6 (KM)
      REAL(r8)  :: WP7 (KM) ,WP8(KM*2)
      REAL(r8)  :: WP9 ,WP10, WP11
      REAL(r8),dimension(IMT,JMT,KM,MAX_BLOCKS_CLINIC) :: WP12,WP13
      REAL(r8)  :: riv1,riv2,epsln,RKV,RKV1
      REAL(r8)  :: adv_x1,adv_x2,adv_x,adv_y1,adv_y2,adv_z,diff_u1,diff_u2,diff_v1,diff_v2
      REAL(r8)  :: dlux,dlvx,dluy,dlvy,dluz,dlvz,adv_z1,adv_z2,adv_z3,adv_z4
      REAL(r6)  :: xxx, c0
      real(r8)  :: hdvk(imt,jmt), hduk(imt,jmt), adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)
!
      type (block) :: this_block          ! block information for current block
       
#if (defined CANUTO)
      REAL(r8)  :: AIDIF
#endif

!YU 
      real (r8):: akt_back(2*km),aks_back(2*km),akm_back(2*km)
      allocate(riu(imt,jmt,0:km,max_blocks_clinic),stat=ierr)
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---riu'
!Yu      stop
!Yu   end if

      c0 = 0.0D0 
 
      epsln = 1.0D-25
      
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JET
         DO I = 1,IMT
            H0BL (I,J,IBLOCK)= H0BF (I,J,IBLOCK)
            H0BF (I,J,IBLOCK)= H0 (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
!lhl0711
#if (defined CANUTO)
       AIDIF=0.0  
!      if (ISC/=0)  AIDIF=0.5  
#endif
!lhl0711
 
!---------------------------------------------------------------------
!     Calculating Richardson number riu (at U/V-point);
!---------------------------------------------------------------------
      s2t  = c0
      ridt = c0
      riu  = c0
      wp12 = c0
      wp13 = c0


!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         call ugrid_to_tgrid(wp12(:,:,k,iblock),up(:,:,k,iblock),iblock)
         call ugrid_to_tgrid(wp13(:,:,k,iblock),vp(:,:,k,iblock),iblock)
      END DO
   END DO
!
 
!$OMP PARALLEL DO PRIVATE (K,J,I,riv1,riv2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 2, JMT
            DO I = 1,IMT-1
               riv1 = wp12 (I,J,K,iblock) - wp12 (I,J,K +1,iblock)
               riv2 = wp13 (I,J,K,iblock) - wp13 (I,J,K +1,iblock)
               s2t (i,j,k,iblock) =vit(i,j,k+1,iblock)*(riv1*riv1+riv2*riv2)*ODZT(K+1)*ODZT(K+1)
               ridt(i,j,k,iblock) =vit(i,j,k+1,iblock)*ricdt(i,j,k,iblock)/(s2t(i,j,k,iblock)+epsln)
#ifdef CANUTO  
               rit (i,j,k,iblock)= VIT (I,J,K +1,iblock)*rict(i,j,k,iblock)/(s2t(i,j,k,iblock)+epsln)
#else          
               rit (i,j,k,iblock)= rit (i,j,k,iblock,iblock) +VIT (I,J,K +1,iblock,iblock)* &
                                   rict(i,j,k,iblock,iblock)/(s2t(i,j,k,iblock,iblock)+epsln)
#endif
            END DO
         END DO
      END DO
   END DO
!
!$OMP PARALLEL DO PRIVATE (K,J,I,riv1,riv2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 1, JMT
            DO I = 1, IMT
               riv1 = UP (I,J,K,iblock) - UP (I,J,K +1,iblock)
               riv2 = VP (I,J,K,iblock) - VP (I,J,K +1,iblock)
               s2u (i,j,k,iblock) =viv(i,j,k+1,iblock)*(riv1*riv1+riv2*riv2)*ODZT(K+1)*ODZT(K+1)
! calculate the shear square and the T component minus S component of Richardson Number
!lhl1204
               riu (i,j,k,iblock) = VIV (I,J,K +1,iblock)*ric (i,j,k,iblock)/(s2u(i,j,k,iblock)+epsln)
!lhl1204
            END DO
         END DO
      END DO
   END DO
!
!
#ifdef CANUTO
!
     amld = c0
     akmt = c0
     akmu = c0
!

!$OMP PARALLEL DO PRIVATE (iblock,wp1,wp2,wp3,wp4,wp5,wp6,wp7,wp8,wp9,wp10,wp11,akt_back,akm_back,aks_back,iwk,wk1,wk2,wk3,xxx)
   do iblock = 1, nblocks_clinic
      DO J = 2, JMT
         DO I = 1, IMT-1

        if (VIT(I,J,1,iblock).gt.0.5) then
!
         wp1=0.D0
         wp2=0.D0
         wp3=0.D0
         wp4=0.D0
         wp5=0.D0
         wp6=0.D0
         wp7=0.D0
         wp8=0.D0
         wp9=0.D0
         wp10=0.D0
         wp11=0.D0
!
         do k=1,km
            AKM_BACK(K)=0.0D0
            AKT_BACK(K)=0.0D0
            AKS_BACK(K)=0.0D0
         end do
!
         DO K = 1,ITNU(i,j,iblock)-1
               wp1(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,1,iblock)-(AT (I,J,K,1,iblock)-AT (I,J,K+1,1,iblock))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))
               wp2(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,2,iblock)-(AT (I,J,K,2,iblock)-AT (I,J,K+1,2,iblock))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))
               wp3(K)= VIT (I,J,K+1,iblock)* (pdensity (I,J,K,iblock)-(pdensity (I,J,K,iblock)- &
                       pdensity(I,J,K+1,iblock))* DZP (K)/(DZP(K)+DZP(K+1)))*1.d-3
               wp8(k)=-vit(i,j,k+1,iblock)*ZKP(k+1)*1.d+2
         END DO
         DO K = 1,ITNU(i,j,iblock)-1
         wp4(k)=vit(i,j,k+1,iblock)*RIT(i,j,k,iblock)
         wp5(k)=vit(i,j,k+1,iblock)*RIDT(i,j,k,iblock)
         wp6(k)=vit(i,j,k+1,iblock)*S2T(i,j,k,iblock)
         wp7(k)=vit(i,j,k+1,iblock)*RICT(i,j,k,iblock)
         END DO
         wp9=vit(i,j,1,iblock)*USTAR(I,J,iblock)*1.0d+2
         wp10=vit(i,j,1,iblock)*BUOYTUR(I,J,iblock)*1.0d+4
         wp11=vit(i,j,1,iblock)*BUOYSOL(I,J,iblock)*1.0d+4
!
         IWK=ITNU(I,J,iblock)-1
!
!input ZKT in cm, AT(1) in C, AT(2) in (s-35)/1000., PDENSITY in g/cm^3


         CALL  TURB_2(wp8,wp1,wp2,wp3,&
!input RIT and RIDT no unit, S2T in 1/s^2, DFRICMX and DWNDMIX in cm^2/s
                wp4,wp5,wp6, DFRICMX*1.0d+4,DWNDMIX*1.0d+4,&
!output in cm^2/s, so 1d-4 should be multipled
               AKM_BACK,AKT_BACK,AKS_BACK,&
!input  RICT in 1/s^2 USTAR in cm/s, BUOYTUR,BUOYSOL in cm^2/s^3,FF in 1/s
               wp7,wp9,wp10,wp11,FF(J),& !OK
!output amld in cm, akmt,akh, and aks in cm^2/s
!               AMLD(I,J,iblock),AKMT(I,J,1,iblock),AKT(I,J,1,1,iblock),AKT(I,J,1,2,iblock),&
               AMLD(I,J,iblock),WK1,WK2,WK3,&
!input int
               IWK,NA(I,J,iblock),KM,1,0,0,i,j) !OK!

         DO K = 1,KM
         xxx = WK1(K)
         WK1(K) = xxx
         xxx = WK2(K)
         WK2(K) = xxx
         xxx = WK3(K)
         WK3(K) = xxx
         AKMT(I,J,K,iblock)=+(WK1(K)+dmin1(AKM_BACK(K),1d-3))*1.d-4
         AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+(WK2(K)+dmin1(AKT_BACK(K),1d-3))/NCC*1.d-4
         AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+(WK3(K)+dmin1(AKS_BACK(K),1d-3))/NCC*1.d-4
         END DO
!
        endif
         END DO
      END DO
   END DO
!!
! calculate the vertical mixing on U-grid
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   do iblock = 1, nblocks_clinic
      DO K = 1,KMM1
         call tgrid_to_ugrid(akmu(:,:,k,iblock), akmt(:,:,k,iblock),iblock)
      END DO
   end do
!lhl241204
#endif
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM: ZONAL COMPONENT
!---------------------------------------------------------------------
      CALL UPWELL (U,V,H0)
      
 
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
     dlu = c0
     dlv = c0
     wka = c0

!$OMP PARALLEL DO PRIVATE (IBLOCK,K)
   do iblock = 1, nblocks_clinic
       DO K = 1,KM
          call tgrid_to_ugrid(wka(:,:,k,iblock),ws(:,:,k,iblock), iblock)
       END DO
   end do

!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERMS
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,adv_x1,adv_x2,adv_y1,adv_y2,adv_z1,adv_z2,adv_z3,adv_z4, &
!$OMP                      dlux,dlvx,dluy,dlvy, dluz,dlvz)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         call advection_momentum(u(:,:,:,iblock),v(:,:,:,iblock),wka(:,:,:,iblock),adv_uu,adv_vv,iblock)
         dlu(:,:,:,IBLOCK)=adv_uu
         dlv(:,:,:,IBLOCK)=adv_vv
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1, IMT
               WKA (I,J,K,IBLOCK)= C0F * SQRT (UP (I,J,K,IBLOCK)* UP (I,J,K,IBLOCK) + VP (I, &
                            J,K,IBLOCK)* VP (I,J,K,IBLOCK))
            END DO
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,k,J,I,diff_u1,diff_v1,diff_u2,diff_v2,rkv,riv1,rkv1)
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
      DO J = 2, JMT-1
         DO I = 2,IMT-1
!lhl1204
#ifdef CANUTO
!       print*,akmU(i,j,k)
            if (k==1) then
               diff_v1 = SV (I,J,IBLOCK)* OD0*(1-AIDIF)
               diff_u1 = SU (I,J,IBLOCK)* OD0*(1-AIDIF)
            else
               diff_v1= AKMU(I,J,K-1,IBLOCK)*(1-AIDIF)*(VP(I,J,K-1,IBLOCK)- VP(I,J,K,IBLOCK))*ODZT(K)*VIV(I,J,K,IBLOCK)+ &
                        (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK)*(1-AIDIF) &
                      * (-SNLAT(J)*UP(I,J,K-1,IBLOCK)*SAG+VP(I,J,K-1,IBLOCK)*CAG)
               diff_u1= AKMU(I,J,K-1,IBLOCK)*(1-AIDIF)* (UP(I,J,K-1,IBLOCK)-UP(I,J,K,IBLOCK))* ODZT (K)*VIV (I,J,K,IBLOCK) + &
                       (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK)*(1-AIDIF) &
                      *(UP(I,J,K-1,IBLOCK)*CAG+SNLAT(J)*VP(I,J,K-1,IBLOCK)* SAG)
            end if
            if (k==km) then
               diff_v2= WKA (I,J,KM,IBLOCK)* ( - SNLAT (J)* UP (I,J,KM,IBLOCK)        &
                        * SAG + VP (I,J,KM,IBLOCK)* CAG)*(1-AIDIF)
               diff_u2= WKA (I,J,KM,IBLOCK)* ( UP (I,J,KM,IBLOCK)* CAG + SNLAT (J)    &
                        * VP (I,J,KM,IBLOCK)* SAG)*(1-AIDIF)
            else
               diff_v2= AKMU(I,J,K,IBLOCK)*(1-AIDIF)*(VP(I,J,K,IBLOCK)- VP(I,J,K+1,IBLOCK))*ODZT(K+1)*VIV(I,J,K+1,IBLOCK)+ &
                        (1.0D0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK)*(1-AIDIF) &
                      * (-SNLAT(J)*UP(I,J,K,IBLOCK)*SAG+VP(I,J,K,IBLOCK)*CAG)
               diff_u2= AKMU(I,J,K,IBLOCK)*(1-AIDIF)*(UP(I,J,K,IBLOCK)-UP(I,J,K+1,IBLOCK))* ODZT(K+1)*VIV (I,J,K+1,IBLOCK) + &
                       (1.0D0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK)*(1-AIDIF) &
                      *(UP(I,J,K,IBLOCK)*CAG+SNLAT(J)*VP(I,J,K,IBLOCK)* SAG)
            end if
#else
            RKV = AMV
            RKV1= AMV
            IF (J_global(j) >= RUST.AND.J_global(j) <= RUEND.AND.K>=2)THEN
!
!        depended on Richardson number
            IF (riu (i,j,k -1,IBLOCK) < 0.0)THEN
                RKV = visc_cbu_limit
            ELSE
                riv1 = 1.0D0/ (1.0D0+5.0D0* riu (i,j,k -1,IBLOCK))
                RKV = fricmx * riv1* riv1+ visc_cbu_back
            END IF
                IF (k == 2.AND.RKV < wndmix) RKV = wndmix
            END IF
!
            IF (J_global(j) >= RUST.AND.J_global(j) <= RUEND)THEN
!
            IF (riu (i,j,k) < 0.0)THEN
                RKV1= visc_cbu_limit
            ELSE
                riv1 = 1.0D0/ (1.0D0+5.0D0* riu (i,j,k,IBLOCK))
                RKV1= fricmx * riv1* riv1+ visc_cbu_back
            END IF
                IF (k == 1.AND.RKV1< wndmix) RKV1= wndmix
            END IF
!
            AKMU(I,J,K,IBLOCK)=RKV1
!
            if (k==1) then
               diff_v1 = SV (I,J,IBLOCK)* OD0
               diff_u1 = SU (I,J,IBLOCK)* OD0
            else
               diff_v1= RKV*(VP(I,J,K-1,IBLOCK)- VP(I,J,K,IBLOCK))*ODZT(K)*VIV(I,J,K,IBLOCK)+ &
                        (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK) &
                      * (-SNLAT(J)*UP(I,J,K-1,IBLOCK)*SAG+VP(I,J,K-1,IBLOCK)*CAG)
               diff_u1= RKV * (UP(I,J,K-1,IBLOCK)-UP(I,J,K,IBLOCK))* ODZT (K)*VIV (I,J,K,IBLOCK) + &
                       (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK) &
                      *(UP(I,J,K-1,IBLOCK)*CAG+SNLAT(J)*VP(I,J,K-1,IBLOCK)* SAG)   
            end if
            if (k==km) then
               diff_v2= WKA (I,J,KM,IBLOCK)* ( - SNLAT (J)* UP (I,J,KM,IBLOCK)        &
                        * SAG + VP (I,J,KM,IBLOCK)* CAG)
               diff_u2= WKA (I,J,KM,IBLOCK)* ( UP (I,J,KM,IBLOCK)* CAG + SNLAT (J)    &
                        * VP (I,J,KM,IBLOCK)* SAG)
            else
               diff_v2= RKV1*(VP(I,J,K,IBLOCK)- VP(I,J,K+1,IBLOCK))*ODZT(K+1)*VIV(I,J,K+1,IBLOCK)+ &
                        (1.0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK) &
                      * (-SNLAT(J)*UP(I,J,K,IBLOCK)*SAG+VP(I,J,K,IBLOCK)*CAG)
               diff_u2= RKV1* (UP(I,J,K,IBLOCK)-UP(I,J,K+1,IBLOCK))* ODZT(K+1)*VIV (I,J,K+1,IBLOCK) + &
                       (1.0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK) &
                      *(UP(I,J,K,IBLOCK)*CAG+SNLAT(J)*VP(I,J,K,IBLOCK)* SAG)   
            end if
!
!
#endif
!lhl1204
            DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + ODZP (K)* (diff_v1-diff_v2)
            DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) + ODZP (K)* (diff_u1-diff_u2)
        END DO
      END DO
      END DO
   END DO

      deallocate(riu) 
 
 
!---------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL VISCOSITY
!---------------------------------------------------------------------
#if ( defined SMAG)
!
         CALL SMAG2 (K)
!
#if (defined SMAG_FZ )
 
#else
!
#endif
!
#else
      
#if (defined BIHAR)

!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
          call hdiffu_del4(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
          do j = 3, jmt-2
          do j = 3, jmt-2
             dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
             dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
          end do
          end do
      END DO
   END DO
#else
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call hdiffu_del2(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
         DO J = 3, JMT-2
            DO I = 3,IMT-2
               dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
               dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
            END DO
         END DO
      END DO
   END DO
 
#endif
#endif
 
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
      allocate(dlub(imt,jmt,max_blocks_clinic),dlvb(imt,jmt,max_blocks_clinic),stat=ierr)
      CALL VINTEG (DLU,DLUB)
      CALL VINTEG (DLV,DLVB)
!
!---------------------------------------------------------------------
!     VERTICAL INTEGRATION
!---------------------------------------------------------------------
 
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---tmp1,tmp2'
!Yu      stop
!Yu   end if
  call mpi_barrier(mpi_comm_ocn,ierr)

      RETURN
      END SUBROUTINE READYC
 
 
