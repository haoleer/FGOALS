!-----------------------------------------------------------------------------
!   Processing some variables from flux coupler.
!-----------------------------------------------------------------------------
!
!        By Yongqiang Yu , 16 Apr., 1999
!
!
!
      SUBROUTINE post_cpl

#include <def-undef.h>


use param_mod
use pconst_mod
use tracer_mod
use forc_mod
use buf_mod
use control_mod
use shr_sys_mod
use output_mod,only:spval
use domain
use grid
use blocks
use POP_HaloMod
use POP_GridHorzMod
use distribution
use gather_scatter

!
      implicit none
      real(r8),dimension(:,:,:),allocatable::tmp_su,tmp_sv
      real(r8) :: ek0, ttt(imt_global,jmt_global)
      integer :: iblock, ErrorCode
!
    type (block) :: this_block          ! block information for current block

         allocate(tmp_su(imt,jmt,max_blocks_clinic))
         allocate(tmp_sv(imt,jmt,max_blocks_clinic))
!
        tsf=0.0_r8
        ssf=0.0_r8
        swv=0.0_r8
        tmp_su=0.0_r8
        tmp_sv=0.0_r8

     n=0
!$OMP PARALLEL DO PRIVATE (iblock,i,j)
    do iblock = 1, nblocks_clinic
        this_block = get_block(blocks_clinic(iblock),iblock)
        do j=this_block%je, this_block%jb, -1
        do i=this_block%ib, this_block%ie
           TSF(i,j,iblock) = (lat1(i,j,iblock)+sen(i,j,iblock)  +lwup(i,j,iblock)+lwdn(i,j,iblock)&
                      +netsw(i,j,iblock)+melth(i,j,iblock)) *OD0CP  ! net heat flux
           SWV(i,j,iblock) =   netsw(i,j,iblock)                               ! net solar radiation
           NSWV(i,j,iblock) = lat1(i,j,iblock)+sen(i,j,iblock)+lwup(i,j,iblock)+lwdn(i,j,iblock)&
                      +melth(i,j,iblock)                                ! none solar radiation !for BUOY
           SSF(i,j,iblock) =  -(prec(i,j,iblock)+evap(i,j,iblock)+&
                         meltw(i,j,iblock)+roff(i,j,iblock)) &
                        *34.7*1.0e-3/DZP(1)*OD0                                     ! P+E+melting !linpf 25->DZP(1)
           tmp_su(i,j,iblock)= taux(i,j,iblock)*cos(anglet(i,j,iblock))-tauy(i,j,iblock)*sin(anglet(i,j,iblock))
           tmp_sv(i,j,iblock)=-tauy(i,j,iblock)*cos(anglet(i,j,iblock))-taux(i,j,iblock)*sin(anglet(i,j,iblock))
        end do
        end do
     end do
        
!!
!!      Update boundary here !!!!
!!
       call POP_HaloUpdate(tsf(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(swv(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(nswv(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(ssf(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(tmp_su(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(tmp_sv(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
#ifdef USE_OCN_CARBON
       call POP_HaloUpdate(pco2(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
#endif         
!!
            lthf = lat1 !latent flux
            sshf = sen !sensible flux
            lwv = lwup + lwdn   !long wave flux
            fresh = ssf
            runoff=roff
!

!!linpf 2012Jul26 !2012Jul28
        where(vit(:,:,1,:)<0.5) tsf=spval
        where(vit(:,:,1,:)<0.5) swv=spval
        where(vit(:,:,1,:)<0.5) nswv=spval
        where(vit(:,:,1,:)<0.5) ssf=spval
        where(vit(:,:,1,:)<0.5) lwv=spval
        where(vit(:,:,1,:)<0.5) sshf=spval
        where(vit(:,:,1,:)<0.5) fresh=spval
        where(vit(:,:,1,:)<0.5) lthf=spval
!
!
        do iblock =1 , nblocks_clinic
            call tgrid_to_ugrid(su(:,:,iblock), tmp_su(:,:,iblock),iblock)
            call tgrid_to_ugrid(sv(:,:,iblock), tmp_sv(:,:,iblock),iblock)
        end do
!
!      Update boundary here !!!!
       call POP_HaloUpdate(su(:,:,:) , POP_haloClinic, POP_gridHorzLocSwcorner,&
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(sv(:,:,:) , POP_haloClinic, POP_gridHorzLocSwcorner,&
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
        ! lihuimin, 2012.7.23, ft. lpf
!!linpf 2012Jul29
        where(viv(:,:,1,:)<0.5) su=spval
        where(viv(:,:,1,:)<0.5) sv=spval
!linpf 2012Jul29
!calculate USTAR

    do iblock = 1, nblocks_clinic
       DO J = 1, jmt
         DO I = 1,imt
          USTAR(I,J,iblock)=sqrt(sqrt(tmp_su(i,j,iblock)*tmp_su(i,j,iblock)+tmp_sv(i,j,iblock)*tmp_sv(i,j,iblock))*OD0)*vit(i,j,1,iblock) 
         END DO
      END DO        
   end do   
!
   licomqice = 0.0_r8
   psa       = 0.0_r8

         deallocate (tmp_su,tmp_sv) !LPF 20120818
        return
        end
