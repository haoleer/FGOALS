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
use blocks
use domain
use LICOM_Error_mod
use gather_scatter
use distribution
      IMPLICIT NONE
 
      integer     :: n2, iblock
      REAL(r8)    :: AIDIF,C2DTTS,AA,FAW,FIW,ALF,RNCC,ABC,fil_lat1,fil_lat2
      REAL(r8)    :: HDTK(imt,jmt), adv_tt(imt,jmt,km),ek0, tttt(imt_global,jmt_global)

!Xiao Chan (Hereinafter XC for short)
      real(r8)    :: LAMDA(imt,jmt,km,max_blocks_clinic),wt1,wt2,adv_y,adv_x,adv_z,adv_x1,adv_x2
      type (block) :: this_block          ! block information for current block

   
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
      DO J = 1, JMT
         DO I = 1,IMT
            H0F (I,J,iblock)= H0F (I,J,iblock)* ONBC
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO iblock = 1, nblocks_clinic
      DO J = 1, JMT
         DO I = 1,IMT
            STF (I,J,IBLOCK)= AA * H0F (I,J,IBLOCK) + (1.0D0- AA)* H0L (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 1, JMT
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
         DO J = 1, JMT
            DO I = 1,IMT
               WKD (I,J,K,IBLOCK)= AA * UTF (I,J,K,IBLOCK) + (1.0D0- AA)* UTL (I,J,K,IBLOCK)
               WKB (I,J,K,IBLOCK)= AA * VTF (I,J,K,IBLOCK) + (1.0D0- AA)* VTL (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO

     write(170+mytid,*) "OK-------1"
     close(170+mytid)
      CALL UPWELL (WKD,WKB,STF)
     write(170+mytid,*) "OK-------2"
     close(170+mytid)
 
#if (defined NODIAG)
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
 
!-----------------------------------------------------------------------
!     PREPARATION FOR ISOPYCNAL DIFFUSION & ADVECTION
!-----------------------------------------------------------------------
     write(170+mytid,*) "OK-------3"
     close(170+mytid)
#if (defined ISO)
      CALL ISOPYC
#endif
     write(170+mytid,*) "OK-------4"
     close(170+mytid)
 
!@@@  COMPUTING DIFFUSION COEFFICIENT
    
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 2, JMT-1
            DO I = 2, IMT-1
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
!-----------------------------------------------------------------------
!     SOLVE FOR ONE TRACER AT A TIME
!-----------------------------------------------------------------------
!     NTRA = 1 => TEMPERATURE
!     NTRA = 2 => SALINITY
 
      DO N = 1,NTRA
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM 
!---------------------------------------------------------------------
     write(170+mytid,*) "OK-------5"
     close(170+mytid)

   do iblock = 1, nblocks_clinic
      adv_tt = 0.0_r8
      call advection_tracer(wkd(:,:,:,iblock),wkb(:,:,:,iblock),ws(:,:,:,iblock),at(:,:,:,n,iblock),adv_tt,iblock)
     write(170+mytid,*) "OK-------6"
     close(170+mytid)
      do k=1, km
      do j =3, jmt-2
      do i =3, imt-2
         tf(i,j,k,iblock) = adv_tt(i,j,k)*vit(i,j,k,iblock)
      end do
      end do
      end do
   end do
!
      do k=1, km
         call gather_global(tttt, tf(:,:,k,:), master_task,distrb_clinic)
         if (mytid ==0) write(161,*) ((tttt(i,j), i=1,imt_global), j=1,jmt_global)
      end do
         if (mytid ==0) close(161)
     write(170+mytid,*) "OK-------6"
     close(170+mytid)
     call chk_var3d(wkd,ek0,0,km)
     if ( mytid ==0) write(160,*) "wkd", ek0
     call chk_var3d(wkb,ek0,0,km)
     if ( mytid ==0) write(160,*) "wkb", ek0
     call chk_var3d(ws,ek0,1,km)
     if ( mytid ==0) write(160,*) "ws", ek0
     call chk_var3d(at(:,:,:,1,:),ek0,1,km)
     if ( mytid ==0) write(160,*) "AT", ek0
     call chk_var3d(tf,ek0,1,km)
     if ( mytid ==0) write(160,*) ek0
    
     if (mytid ==0 ) then
        write(123,*)((tf(i,j,1,1), i=3,imt-2),j=6,8)
        close(123)
     end if
#ifdef CANUTO      
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 2, JMT-1
            DO I = 2, IMT-1
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
#else
 
!-----------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL DIFFUSION TERMS:  ZONAL and Meridional COMPONENTS
!-----------------------------------------------------------------------
 
#if (defined BIHAR)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,K)
   do iblock = 1, nblocks_clinic
   this_block = get_block(blocks_clinic(iblock),iblock)
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
   this_block = get_block(blocks_clinic(iblock),iblock)
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
     call chk_var3d(tf,ek0,1,km)
     if ( mytid ==0) write(160,*) ek0
 

         IF (N == 1)THEN
#if (defined SOLAR)
!     SOLAR SHORTWAVE PENETRATION
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
    DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K=2,KM-1
        DO J=3, JMT-2
        DO I=3,IMT-2
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
        DO J=3, JMT-2
        DO I=3,IMT-2
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
        DO J=3,JMT-2
        DO I=3,IMT-2
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
        DO J=3,JMT-2
        DO I=3,IMT-2
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
     call chk_var3d(tf,ek0,1,km)
     if ( mytid ==0) write(160,*) ek0
 
        wt1=0
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K=2,KM-1
           DO J=3,JMT-2
           DO I=3,IMT-2
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
               if (mytid == 0 ) then
                     write(145,*) ((tf(i,j,1,1), i=3,imt-2),j=6,8)
                     close(145)
               end if
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
       DO J = 3,JMT-2
       DO I = 3,IMT-2
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
    
               if (mytid == 0 ) then
                     write(146,*) ((wkc(i,j,1,1), i=3,imt-2),j=6,8)
                     close(146)
               end if
!-----------------------------------------------------------------------
!     SET NEWTONIAN SURFACE BOUNDARY CONDITION
!-----------------------------------------------------------------------
     call chk_var3d(tf,ek0,1,km)
     if ( mytid ==0) write(160,*) ek0
 
         IF (N == 2)THEN
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)    
        DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 3,JMT-2
               DO I = 3,IMT-2
                  IF (KMT(I,J,IBLOCK) > 0)THEN
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
               if (mytid == 0 ) then
                     write(144,*) odzp(1), AIDIF
                     write(144,*) ((stf(i,j,1), i=3,imt-2),j=6,8)
                     close(144)
               end if
!
  ELSE
!$OMP PARALLEL DO PRIVATE (J,I)
        DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 3,JMT-2
               DO I = 3,IMT-2
                IF (KMT(I,J,IBLOCK) > 0)THEN
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
         do j=3,jmt-2
               DO I = 3,jmt-2
                  TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (RESTORE (I,J,K,&
                             N,iblock) - ATB (I,J,K,N,iblock))*LAMDA(i,j,k,iblock)
               END DO
            END DO
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 2,KM
         do j=3,jmt-2
               DO I = 3, imt-2
                  TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (RESTORE (I,J,K,&
                             N,iblock) - ATB (I,J,K,N,iblock))*LAMDA(i,j,k,iblock)
               END DO
            END DO
         END DO
    END DO
 
#endif
!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
 
       if (mytid == 0) then
          write(141,*) ((tf(i,j,1,1),i=3,imt-2),j=6,8)
          close(141)
       end if
       if (mytid == 0) then
          write(142,*) ((tf(i,j,3,1),i=3,imt-2),j=6,8)
          close(142)
       end if
     call chk_var3d(tf,ek0,1,km)
     if ( mytid ==0) write(160,*) ek0
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 3, JMT-2
               DO I = 3, IMT-2
                  if (trim(adv_tracer) == 'tspas') then
                     VTL (I,J,K,iblock) = AT (I,J,K,N,iblock) + DTS * TF (I,J,K,iblock)
                  else if (trim(adv_tracer) == 'centered') then
                     VTL (I,J,K,iblock) = ATB (I,J,K,N,iblock) + C2DTTS * TF (I,J,K,iblock)
                  else 
                     call exit_licom(sigAbort,'The false advection option for tracer')
                  end if
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
!
       if (mytid == 0) then
          write(143,*) ((vtl(i,j,3,1),i=3,imt-2),j=6,8)
          close(143)
       end if
     call chk_var3d(vtl,ek0,1,km)
     if ( mytid ==0) write(160,*) ek0
     call POP_HaloUpdate(VTL , POP_haloClinic, POP_gridHorzLocCenter,&
                         POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
               
         if (mod(ist,180) == 1) then
            CALL SMTS (VTL,VIT,fil_lat2)
         else
             DO IBLOCK = 1, NBLOCKS_CLINIC
             DO J=1, JMT
             DO I=1,IMT  
                  if (trim(adv_tracer) == 'tspas') then
                     VTL (I,J,1,IBLOCK) = VTL(I,J,1,IBLOCK)-AT (I,J,1,N,IBLOCK) - NET(I,J,N,IBLOCK)*DTS
                  else if (trim(adv_tracer) == 'centered') then
                     VTL (I,J,1,IBLOCK) = VTL(I,J,1,IBLOCK)- ATB (I,J,1,N,IBLOCK) - NET(I,J,N,IBLOCK)*C2DTTS
                  else 
                     call exit_licom(sigAbort,'The false advection option for tracer')
                  end if
             END DO
             END DO
             END DO
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
        DO  IBLOCK = 1, NBLOCKS_CLINIC
            DO K=2,KM
            DO J=1, JMT
            DO I=1,IMT  
               if (trim(adv_tracer) == 'tspas') then
                  VTL (I,J,K,IBLOCK) = VTL(I,J,K,IBLOCK)-AT (I,J,K,N,IBLOCK)
               else if (trim(adv_tracer) == 'centered') then
                  VTL (I,J,K,IBLOCK) = VTL(I,J,K,IBLOCK)- ATB (I,J,K,N,IBLOCK)
               else 
                     call exit_licom(sigAbort,'The false advection option for tracer')
               end if
            END DO
            END DO
            END DO
       END DO
                  
           CALL SMTS (VTL,VIT,fil_lat1)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
         DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J=1, JMT
            DO I=1,IMT  
               if (trim(adv_tracer) == 'tspas') then
                  VTL (I,J,1,IBLOCK) = AT (I,J,1,N,IBLOCK) + VTL(I,J,1,IBLOCK) + NET(I,J,N,IBLOCK)*DTS
               else if (trim(adv_tracer) == 'centered') then
                  VTL (I,J,1,IBLOCK) = ATB (I,J,1,N,IBLOCK) + VTL(I,J,1,IBLOCK) + NET(I,J,N,IBLOCK)*C2DTTS
               else 
                     call exit_licom(sigAbort,'The false advection option for tracer')
               end if

            END DO
            END DO
        END DO

!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
         DO IBLOCK = 1, NBLOCKS_CLINIC
            DO K=2,KM
            DO J=1, JMT
            DO I=1,IMT  
              if ( trim(adv_tracer) == 'tspas') then
                  VTL (I,J,K,IBLOCK) = AT (I,J,K,N,IBLOCK) + VTL(I,J,K,IBLOCK)
              else
                  VTL (I,J,K,IBLOCK) = ATB (I,J,K,N,IBLOCK) + VTL(I,J,K,IBLOCK)
              end if
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
            DO J = 1, JMT
               DO I = 1,IMT
                 trend(i,j,k,N,iblock)=(VTL(I,J,K,iblock)-ATB(I,J,K,N,iblock))/C2DTTS*VIT(I,J,K,iblock)
               END DO
            END DO
         END DO
     END DO
 
                     VTL (I,J,K,iblock) = AT (I,J,K,N,iblock) + DTS * TF (I,J,K,iblock)
                     VTL (I,J,K,iblock) = ATB (I,J,K,N,iblock) + C2DTTS * TF (I,J,K,iblock)

  if (trim(adv_tracer) == 'tspas') then
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1, JMT
               DO I = 1,IMT
                  ATB(I,J,K,N,IBLOCK)=AT (I,J,K,N,IBLOCK)
                  AT (I,J,K,N,IBLOCK) = VTL (I,J,K,IBLOCK)
               END DO
            END DO
         END DO
      END DO
  else if (trim(adv_tracer) == 'centered') then
         IF (IST >= 1)THEN
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
            DO K = 1,KM
               DO J = 1, JMT
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
            DO J = 1, JMT
               DO I = 1,IMT
                  AT (I,J,K,N,IBLOCK) = VTL (I,J,K,IBLOCK)
               END DO
            END DO
         END DO
      END DO
  else 
      call exit_licom(sigAbort,'The false advection option for tracer')
  end if
 call chk_var3d(vtl,ek0,1,km)
     if ( mytid ==0) write(160,*) ek0
   stop
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
  call mpi_barrier(mpi_comm_ocn,ierr)
      RETURN
      END SUBROUTINE TRACER
 
