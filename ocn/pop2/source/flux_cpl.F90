!
!-----------------------------------------------------------------------------
!   Preparing oceanic variables for communication between OGCM and coupler.
!-----------------------------------------------------------------------------
!
!        By Yongqiang Yu , 16 Apr., 1999
!

module fluxcpl

#include <def-undef.h>   

      use param_mod
      use pconst_mod
      use buf_mod
      use tracer_mod
      use dyn_mod
      use cdf_mod
      use control_mod
      use msg_mod,only : nproc
      use shr_const_mod,only:SHR_CONST_SPVAL !linpf 20120816
      use domain
      use grid
#ifdef USE_OCN_CARBON      
      use coutput_mod, only : uptake
#endif      
!
!
contains

  SUBROUTINE flux_cpl

!
        implicit none
!
    integer :: iblock
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    do iblock = 1, nblocks_clinic
        do j=1,jmt
        do i=1,imt
           if (licomqice(i,j,iblock) .gt. 0.0) then
               q(i,j,iblock)= licomqice(i,j,iblock)*D0*CP*DZP(1)/86400.
           else
               q(i,j,iblock)= (tbice-at(i,j,1,1,iblock))*D0*CP*DZP(1)/86400.  
           endif
        end do
        end do
    end do
!
        T_CPL = AT(:,:,1,1,:) + 273.15
        S_CPL = AT(:,:,1,2,:)*1000. + 35.
        U_CPL = 0.
        V_CPL = 0.
!linpf 2012Jul26
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    do iblock = 1, nblocks_clinic
        do j=1,jmt
         do i=1,imt
           if(vit(i,j,1,iblock)<0.5) t_cpl(i,j,iblock)=SHR_CONST_SPVAL  
           if(vit(i,j,1,iblock)<0.5) s_cpl(i,j,iblock)=SHR_CONST_SPVAL
         enddo
        enddo
     end do
 
!$OMP PARALLEL DO PRIVATE (IBLOCK)
    do iblock = 1, nblocks_clinic
      call ugrid_to_tgrid(u_cpl(:,:,iblock),u(:,:,1,iblock),iblock)
      call ugrid_to_tgrid(v_cpl(:,:,iblock),v(:,:,1,iblock),iblock)
    end do
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    do iblock = 1, nblocks_clinic
        do j= 2, jmt-1
        do i= 2, imt-1
           dhdx (i,j,iblock)  =   (h0(i+1,j,iblock)-h0(i-1,j,iblock)) /(hun(i,j)+hun(i+1,j))
           dhdy (i,j,iblock)  =   (h0(i,j+1,iblock)-h0(i,j-1,iblock)) /(hue(i,j)+hue(i,j-1))
        end do
        end do
     end do

        ! TODO, consider here
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    do iblock = 1,nblocks_clinic
        do j=1,jmt
         do i=1,imt
           if(vit(i,j,1,iblock)<0.5) dhdx(i,j,iblock)=0.0_r8 !SHR_CONST_SPVAL  
           if(vit(i,j,1,iblock)<0.5) dhdy(i,j,iblock)=0.0_r8 !SHR_CONST_SPVAL
         enddo
        enddo
    enddo
!
!
     licomqice = 0.0
        
#ifdef USE_OCN_CARBON
      co2_cpl(:,:) = uptake(:,:)
#endif          
        return

  END subroutine flux_cpl


end module fluxcpl

