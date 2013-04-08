!  CVS: $Id: tracer.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE TRACER
!     =================
#include <def-undef.h>
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use work_mod
use dyn_mod
use isopyc_mod
use forc_mod
use pmix_mod
use msg_mod
use smuvh
use advection
      IMPLICIT NONE
 
      integer     :: n2, iblock
      REAL(r8)    :: AIDIF,C2DTTS,AA,FAW,FIW,ALF,RNCC,ABC,fil_lat1,fil_lat2
      REAL(r8)    :: HDTK(imt,jmt), adv_tt(imt,jmt,km)

!Xiao Chan (Hereinafter XC for short)
#if (defined TSPAS)
      real(r8)    :: LAMDA,wt1,wt2,adv_z
      real(r8),dimension(:,:,:,:), allocatable :: adv_xy1,adv_xy2,adv_xy3,adv_xy4
      real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: uaa,vaa, &
                adv_x0,adv_y0,adv_c1,adv_c2,atmax,atmin,at0,adv_xx,adv_yy
      real(r8),dimension(:,:,:,:) , allocatable :: adv_zz,atz, adv_za,adv_zb1,adv_zb2,adv_zc,atmaxz,atminz
#else
      real(r8)    :: LAMDA,wt1,wt2,adv_y,adv_x,adv_z,adv_x1,adv_x2
#endif
!XC
 
!---------------------------------------------------------------------
!     SET LOCAL CONSTANT
!---------------------------------------------------------------------
      deallocate(dlu,dlv,gg,ric,rict)
 
      allocate(stf(imt,jmt,max_blocks_clinic),tf(imt,jmt,km,max_blocks_clinic))
      allocate(wkb(imt,jmt,km,max_blocks_clinic),wkc(imt,jmt,km,max_blocks_clinic),wkd(imt,jmt,km,max_blocks_clinic))
!lhl      AIDIF = 0.5*FLOAT(ISOP)
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=63.0D0
       fil_lat2=63.0D0

#if (defined ISO)
      AIDIF = 0.5D0
#else
      AIDIF = 0.0D0
#endif
 
      LAMDA=1.0D0/(15.D0*86400.D0)
      RNCC = 1.0D0/ FLOAT (NCC)
      IF (IST >= 1)THEN
         C2DTTS = DTS *2.0D0
         AA = 0.5D0
      ELSE
         C2DTTS = DTS
         AA = 0.0D0
      END IF
 

!---------------------------------------------------------------------
!     PREPARATION FOR VERTICAL ADVECTIVE TERM
!---------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (J,I)
   DO iblock = 1, nblocks_clinic
      DO J = JST,JET
         DO I = 1,IMT
            H0F (I,J,iblock)= H0F (I,J,iblock)* ONBC
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO iblock = 1, nblocks_clinic
      DO J = JST,JET
         DO I = 1,IMT
            STF (I,J,IBLOCK)= AA * H0F (I,J,IBLOCK) + (1.0D0- AA)* H0L (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK)* ONCC
               VTF (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK)* ONCC
            END DO
         END DO
      END DO
  END DO

 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
               WKD (I,J,K,IBLOCK)= AA * UTF (I,J,K,IBLOCK) + (1.0D0- AA)* UTL (I,J,K,IBLOCK)
               WKB (I,J,K,IBLOCK)= AA * VTF (I,J,K,IBLOCK) + (1.0D0- AA)* VTL (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO

!XC
#if (defined TSPAS)
!     
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = JSM,JEM
            DO I = 2,IMM
               uaa(i,j,k,iblock)=0.5D0*WKD (i,j,k,iblock) + 0.5D0*WKD (i,j-1,k,iblock)
               vaa(i,j,k,iblock)=0.5D0*WKC (i,j,k,iblock) + 0.5D0*WKC (i+1,j,k,iblock)
            END DO
         END DO
      END DO
   END DO
!XC
#endif
 
      CALL UPWELL (WKD,WKB,STF)
 
#if (defined NODIAG)
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
 
!-----------------------------------------------------------------------
!     PREPARATION FOR ISOPYCNAL DIFFUSION & ADVECTION
!-----------------------------------------------------------------------
#if (defined ISO)
      CALL ISOPYC
#endif
 
!@@@  COMPUTING DIFFUSION COEFFICIENT
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
#if (defined ISO)
               WKC (I,J,K,IBLOCK) = AHV + AHISOP * K3 (I,K,J,3,IBLOCK)
#else
               WKC (I,J,K,IBLOCK) = AHV
#endif
!
#if (!defined CANUTO)
            AKT(I,J,K,1,IBLOCK)=WKC(I,J,K,IBLOCK)
            AKT(I,J,K,2,IBLOCK)=WKC(I,J,K,IBLOCK)
#endif
!
            END DO
         END DO
      END DO
   END DO
 
!
#if (!defined CANUTO)
#ifdef SPMD
      call exch_boundary(akt(1,1,1,1),km)
      call exch_boundary(akt(1,1,1,2),km)
#endif
#endif
!
!-----------------------------------------------------------------------
!     SOLVE FOR ONE TRACER AT A TIME
!-----------------------------------------------------------------------
!     NTRA = 1 => TEMPERATURE
!     NTRA = 2 => SALINITY
 
      DO N = 1,NTRA
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM 
!---------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
      adv_tt = 0.0_r8
      call advection_tracer(wkd(:,:,:,iblock),wkb(:,:,:,iblock),ws(:,:,1:km,iblock),at(:,:,:,n,iblock),adv_tt,iblock)
      do k=1, km
      do j =3, jmt-2
      do i =3, imt-2
         tf(i,j,k,iblock) = adv_tt(i,j,k)
      end do
      end do
      end do
   end do
!
#ifdef CANUTO      
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET          
            DO I = 2,IMM
#if (defined ISO)       
               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK) + AHISOP * K3 (I,K,J,3,IBLOCK) 
#else            
               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK) 
#endif         
            END DO
         END DO
      END DO
   END DO
!
#endif

!-----------------------------------------------------------------------
!     COMPUTE THE ISOPYCNAL/DIPYCNAL MIXING
!-----------------------------------------------------------------------
!     XZ AND YZ ISOPYCNAL DIFFUSIVE FLUX ARE SOLVED EXPLICITLY;
!     WHILE ZZ COMPONENT WILL BE SOLVED IMPLICITLY.
 
 
#if (defined ISO)
         CALL ISOFLUX (N)
#else 
 
#if ( defined SMAG)
         CALL SMAG3
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = JSM,JEM
               WKI (1)= 0.0D0
               DO I = 2,IMT
                  WKI (I) = 0.5D0* (AH3 (I,J,K,iblock) + AH3 (I -1,J,K,iblock))* (ATB ( &
                           I,J,K,N,iblock) - ATB (I -1,J,K,N,iblock))* VIT (I,J,K,iblock)* VIT (I -1,J,K,iblock)
               END DO
               DO I = 2,IMM
                  TF (I,J,K,iblock) = TF (I,J,K) + SOTX (J)*(WKI(I+1)-WKI(I))
!
!                dx(i,j,k,n,iblock)= SOTX (J)*(WKI(I+1)-WKI(I))
!
               END DO
            END DO
         END DO
     END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,wt1,wt2)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
         DO J = JSM,JEM
         DO I = 2,IMM
            wt1= 0.5D0*(AH3(I,J,K,IBLOCK)+AH3(I,J-1,K,IBLOCK))*(ATB(I,J,K,N,IBLOCK)- &
                        ATB(I,J-1,K,N,IBLOCK))*VIT(I,J,K,IBLOCK)*VIT(I,J-1,K,IBLOCK)
            wt2= 0.5D0*(AH3(I,J+1,K,IBLOCK)+AH3(I,J,K,IBLOCK))*(ATB(I,J+1,K,N,IBLOCK)- &
                        ATB(I,J,K,N,IBLOCK))*VIT(I,J+1,K,IBLOCK)*VIT(I,J,K,IBLOCK)
            TF (I,J,K,IBLOCK) = TF (I,J,K) + (R2D (J)* wt2 - R2C(J)* wt1)
!
!           dy(i,j,k,n,iblock)=(R2D(J)*wt2-R2C(J)*wt1)
!          ddy(i,j,k,n,iblock)=(R2A(J)*wt2-R2B(J)*wt1)
!
         END DO
         END DO
         END DO
      END DO
 
#else
 
!-----------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL DIFFUSION TERMS:  ZONAL and Meridional COMPONENTS
!-----------------------------------------------------------------------
 
#if (defined BIHAR)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,K)
   do iblock = 1, nblocks_clinic
   do k =1, km
      call hdifft_del4(k,HDTK,ATB,this_block,ntracer)
      do j= 3, jmt-2
      do i= 3, imt-2
           TF (I,J,K,IBLOCK) = TF (I,J,K,IBLOCK) + HDTK(I,J)
      end do
      end do
   end do
   end do
!
#else
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,K)
   do iblock = 1, nblocks_clinic
   do k =1, km
      call hdifft_del2(k,HDTK,ATB,this_block,ntracer)
      do j= 3, jmt-2
      do i= 3, imt-2
           TF (I,J,K,IBLOCK) = TF (I,J,K,IBLOCK) + HDTK(I,J)
      end do
      end do
   end do
   end do
!
#endif
#endif
#endif

!-----------------------------------------------------------------------
!     VERTICAL COMPONENT
!-----------------------------------------------------------------------
 
         IF (N == 1)THEN
#if (defined SOLAR)
!     SOLAR SHORTWAVE PENETRATION
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
    DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K=2,KM-1
        DO J=JSM,JEM
        DO I=2,IMM
            wt1= SWV (I,J,IBLOCK)*pen(k-1)*VIT(I,J,K,IBLOCK)
            wt2= SWV (I,J,IBLOCK)*pen(k)*VIT(I,J,K+1,IBLOCK)
            TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+(wt1-wt2)*ODZP(K)
!
!lhl0105           penetrate(i,j,k)=(wt1-wt2)*ODZP(k)
            penetrate(i,j,k,IBLOCK)= wt2*ODZP(k)
!
        END DO
        END DO
        END DO
     END DO
!$OMP PARALLEL DO PRIVATE (J,I,wt1,wt2)
     DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J=JSM,JEM
        DO I=2,IMM
           wt1= SWV (I,J,IBLOCK)*pen(1)*VIT(I,J,2,IBLOCK)
           wt2= SWV (I,J,IBLOCK)*pen(km-1)*VIT(I,J,km,IBLOCK)
           TF(I,J,1,IBLOCK)=TF(I,J,1,IBLOCK)-ODZP(1)*wt1
           TF(I,J,km,IBLOCK)=TF(I,J,km,IBLOCK)+ODZP(km)*wt2
!
            penetrate(i,j, 1,IBLOCK)= wt1*ODZP( 1)
            penetrate(i,j,km,IBLOCK)= wt2*ODZP(km)
!
        END DO
        END DO
     END DO
#endif


!==================================
!From here ie. the code about SOLARCHLORO added by linpf
!==================================

#if (defined SOLARCHLORO)
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K=2,KM-1
        DO J=JSM,JEM
        DO I=2,IMM
            wt1= SWV(I,J,IBLOCK)*pen_chl(I,J,K-1,IBLOCK)*VIT(I,J,K,IBLOCK)
            wt2= SWV(I,J)*pen_chl(I,J,K,IBLOCK)*VIT(I,J,K+1,IBLOCK)
            TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+(wt1-wt2)*ODZP(K)
!
            penetrate(I,J,K,IBLOCK)=(wt1-wt2)*ODZP(K)
!
        END DO
        END DO
      END DO
  END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J=JSM,JEM
        DO I=2,IMM
           wt1= SWV(I,J,iblock)*pen_chl(I,J,1,iblock)*VIT(I,J,2,iblock)
           wt2= SWV(I,J,iblock)*pen_chl(I,J,km-1,iblock)*VIT(I,J,km,iblock)
           TF(I,J,1,iblock)=TF(I,J,1,iblock)-ODZP(1)*wt1
           TF(I,J,km,iblock)=TF(I,J,km,iblock)+ODZP(km)*wt2
!
            penetrate(I,J, 1,iblock)=-wt1*ODZP(1)
            penetrate(I,J,km,iblock)= wt2*ODZP(km)
!
        END DO
        END DO
  END DO
#endif

!==================================
!above code about SOLARCHLORO added by linpf
!==================================
         END IF
 
!     EDDY-DIFFUSION
 
        wt1=0
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K=2,KM-1
           DO J=JSM,JEM
           DO I=2,IMM
              wt1= WKC(I,J,K-1,IBLOCK)*(ATB(I,J,K-1,N,IBLOCK)-ATB(I,J,K,N,IBLOCK))*ODZT(K)*VIT(I,J,K,IBLOCK)
              wt2= WKC(I,J,K,IBLOCK)*(ATB(I,J,K,N,IBLOCK)-ATB(I,J,K+1,N,IBLOCK))*ODZT(K+1)*VIT(I,J,K+1,IBLOCK)
              TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+ODZP(K)*(wt1-wt2)*(1.0D0-AIDIF)
!
              dz(i,j,k,N,IBLOCK)=ODZP(K)*(wt1-wt2)*(1.0D0-AIDIF)
!
           END DO
           END DO
        END DO
  END DO
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
       DO J = JSM,JEM
       DO I = 2,IMM
           wt1= WKC(I,J,1,IBLOCK)*(ATB(I,J,1,N,IBLOCK)-ATB(I,J,2,N,IBLOCK))*ODZT(2)*VIT(I,J,2,IBLOCK)
           wt2= WKC(I,J,km-1,IBLOCK)*(ATB(I,J,km-1,N,IBLOCK)-ATB(I,J,km,N,IBLOCK))*ODZT(km)*VIT(I,J,km,IBLOCK)
           TF(I,J,1,IBLOCK)=TF(I,J,1,IBLOCK)-ODZP(1)*wt1*(1.0D0-AIDIF)
           TF(I,J,km,IBLOCK)=TF(I,J,km,IBLOCK)+ODZP(km)*wt2*(1.0D0-AIDIF)
!
              dz(i,j, 1,N,IBLOCK)=-ODZP( 1)*wt1*(1.0-AIDIF)
              dz(i,j,km,N,IBLOCK)= ODZP(km)*wt2*(1.0-AIDIF)
!
       END DO
       END DO
  END DO
 
!-----------------------------------------------------------------------
!     SET NEWTONIAN SURFACE BOUNDARY CONDITION
!-----------------------------------------------------------------------
 
         IF (N == 2)THEN
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)    
        DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = JSM,JEM
               DO I = 2,IMM
                  IF (ITNU (I,J,IBLOCK) > 0)THEN
#ifdef COUP
!                     STF (I,J,IBLOCK) = SSF(I,J,IBLOCK)
                     STF (I,J,IBLOCK) = SSF(I,J,IBLOCK)/ODZP(1)
#else
#ifdef FRC_CORE
                     STF (I,J,IBLOCK) = (fresh(i,j)*35.0/1000.0/1000.0&
                   +GAMMA*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))*seaice(i,j,IBLOCK)/ODZP(1)&
                   +GAMMA*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))/ODZP(1)/12.*(1.0-seaice(i,j,IBLOCK)))
#else
!                     STF (I,J,IBLOCK) = GAMMA * (SSS (I,J,IBLOCK) - ATB (I,J,1,2,IBLOCK))
                      STF (I,J,IBLOCK) = GAMMA * (SSS (I,J,IBLOCK) - ATB (I,J,1,2,IBLOCK))/ODZP(1)
#endif
#endif
!                     TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + STF (I,J,IBLOCK)* (1.0- AIDIF)
                     TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + STF (I,J,IBLOCK)* (1.0- AIDIF)*ODZP(1)
!
!                     NET (I,J,2,IBLOCK) = STF (I,J,IBLOCK)*(1.0-AIDIF)
                     NET (I,J,2,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)
!
                  END IF
               END DO
            END DO
      END DO
!
  ELSE
!$OMP PARALLEL DO PRIVATE (J,I)
        DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = JSM,JEM
               DO I = 2,IMM
                IF (ITNU (I,J,IBLOCK) > 0)THEN
#ifdef COUP
                 STF (I,J,IBLOCK) = TSF(I,J,IBLOCK)
#else
#ifdef FRC_CORE
                 STF (I,J,IBLOCK) = ((SWV(I,J,IBLOCK)+NSWV(I,J,IBLOCK))*OD0CP+ & 
                                      SEAICE(I,J,IBLOCK)*GAMMA*(SST(I,J,IBLOCK)-ATB(I,J,1,1,IBLOCK))/ODZP(1))
#else
                 STF (I,J,IBLOCK) = (SWV(I,J,IBLOCK)+NSWV(I,J,IBLOCK)-DQDT(I,J,IBLOCK)* & 
                                    (SST (I,J,IBLOCK) - ATB (I,J,1,1,IBLOCK)))*OD0CP
#endif
#endif
                 TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + STF (I,J,IBLOCK)* ODZP (1)* (1.0D0- AIDIF)
!
!                 NET (I,J,1,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)*(1.0-AIDIF)
                 NET (I,J,1,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)
!
                END IF
               END DO
            END DO
        END DO
    END IF
 
#if (defined BOUNDARY)
!-----------------------------------------------------------------------
!    boundary condition
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 2,KM
         do j=1,jmt
            if (j_global(j)>=(jst_global+1).and.j_global(j)<=(jst_global+50)) then
               DO I = 2,IMM
                  TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (RESTORE (I,J,K,&
                             N,iblock) - ATB (I,J,K,N,iblock))*LAMDA
               END DO
             endif
            END DO
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 2,KM
         do j=1,jmt
            if (j_global(j)<=(jmt_global-1).and.j_global(j)>=(jmt_global-50)) then
               DO I = 2,IMM
                  TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (RESTORE (I,J,K,&
                             N,iblock) - ATB (I,J,K,N,iblock))*LAMDA
               END DO
             endif
            END DO
         END DO
    END DO
 
#endif
!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 2,IMM
!XC
#if (defined TSPAS)
                  VTL (I,J,K,iblock) = AT (I,J,K,N,iblock) + DTS * TF (I,J,K,iblock)
#else
                  VTL (I,J,K,iblock) = ATB (I,J,K,N,iblock) + C2DTTS * TF (I,J,K,iblock)
#endif
!XC
               END DO
            END DO
         END DO
     END DO
 
!-----------------------------------------------------------------------
!     ADD DT/DT COMPONENT DUE TO IMPLICIT VERTICAL DIFFUSION
!-----------------------------------------------------------------------
 
#if (defined ISO)
         CALL INVTRI (VTL,STF,WKC,AIDIF,C2DTTS)
#endif
 
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
               
         if (mod(ist,180) == 1) then
            CALL SMTS (VTL,VIT,KM,fil_lat2)
         else
             DO IBLOCK = 1, NBLOCKS_CLINIC
             DO J=JSM,JEM
             DO I=1,IMT  
#if (defined TSPAS)
                  VTL (I,J,1,IBLOCK) = VTL(I,J,1,IBLOCK)-AT (I,J,1,N,IBLOCK) - NET(I,J,N,IBLOCK)*DTS
#else
                  VTL (I,J,1,IBLOCK) = VTL(I,J,1,IBLOCK)- ATB (I,J,1,N,IBLOCK) - NET(I,J,N,IBLOCK)*C2DTTS
#endif
             END DO
             END DO
             END DO
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
        DO  IBLOCK = 1, NBLOCKS_CLINIC
            DO K=2,KM
            DO J=JSM,JEM
            DO I=1,IMT  
#if (defined TSPAS)
                  VTL (I,J,K,IBLOCK) = VTL(I,J,K,IBLOCK)-AT (I,J,K,N,IBLOCK)
#else
                  VTL (I,J,K,IBLOCK) = VTL(I,J,K,IBLOCK)- ATB (I,J,K,N,IBLOCK)
#endif
            END DO
            END DO
            END DO
       END DO
                  
           CALL SMTS (VTL,VIT,KM,fil_lat1)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
         DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J=JSM,JEM
            DO I=1,IMT  
#if (defined TSPAS)
                  VTL (I,J,1,IBLOCK) = AT (I,J,1,N,IBLOCK) + VTL(I,J,1,IBLOCK) + NET(I,J,N,IBLOCK)*DTS
#else
                  VTL (I,J,1,IBLOCK) = ATB (I,J,1,N,IBLOCK) + VTL(I,J,1,IBLOCK) + NET(I,J,N,IBLOCK)*C2DTTS
#endif
            END DO
            END DO
        END DO

!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
         DO IBLOCK = 1, NBLOCKS_CLINIC
            DO K=2,KM
            DO J=JSM,JEM
            DO I=1,IMT  
#if (defined TSPAS)
                  VTL (I,J,K,IBLOCK) = AT (I,J,K,N,IBLOCK) + VTL(I,J,K,IBLOCK)
#else
                  VTL (I,J,K,IBLOCK) = ATB (I,J,K,N,IBLOCK) + VTL(I,J,K,IBLOCK)
#endif
            END DO
            END DO
            END DO
        END DO
        end if  
!
!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = JSM,JEM
               DO I = 1,IMT
                 trend(i,j,k,N,iblock)=(VTL(I,J,K,iblock)-ATB(I,J,K,N,iblock))/C2DTTS*VIT(I,J,K,iblock)
               END DO
            END DO
         END DO
     END DO
 
#if (defined TSPAS)
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = JST,JET
               DO I = 1,IMT
                  ATB(I,J,K,N,IBLOCK)=AT (I,J,K,N,IBLOCK)
                  AT (I,J,K,N,IBLOCK) = VTL (I,J,K,IBLOCK)
               END DO
            END DO
         END DO
      END DO
#else
         IF (IST >= 1)THEN
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
            DO K = 1,KM
               DO J = JST,JET
                  DO I = 1,IMT
                     ATB (I,J,K,N,IBLOCK) = AFT2* AT (I,J,K,N,IBLOCK) + AFT1* (ATB (I,&
                                    J,K,N,IBLOCK) + VTL (I,J, K,IBLOCK))
                  END DO
               END DO
            END DO
     END DO
         END IF
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = JST,JET
               DO I = 1,IMT
                  AT (I,J,K,N,BLOCK) = VTL (I,J,K,BLOCK)
               END DO
            END DO
         END DO
      END DO
#endif
   END DO
!XC

#else
         atb(:,:,:,1,:)= 12.0D0
         atb(:,:,:,2,:)= c0
         at (:,:,:,:,:)= atb(:,:,1:km,:,:)
#endif
 
      IST = IST +1
 
#ifdef ISO
      deallocate(K1,K2,K3,adv_vetiso,adv_vbtiso,adv_vntiso)
#endif
      deallocate(stf,tf)
      deallocate(wkb,wkc,wkd)
      deallocate(rit)
      RETURN
      END SUBROUTINE TRACER
 
