!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module advection

!BOP
! !MODULE: advection
!
! !DESCRIPTION:
!  This module contains arrays and variables necessary for performing
!  advection of momentum and tracer quantities.  Currently, the
!  module supports leapfrog centered advection of momentum and
!  both leapfrog centered advection and third-order upwinding of
!  tracers.
!
! !REVISION HISTORY:
!  SVN:$Id: advection.F90 28439 2011-05-18 21:40:58Z njn01 $

! !USES:
   use precision_mod
   use param_mod
   use pconst_mod
   use constant_mod
   use grid
   use LICOM_Error_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public :: advection_momentum, &
             advection_tracer

!EOP
!BOC

    contains

!-----------------------------------------------------------------------
!EOC
      subroutine advection_momentum(uuu,vvv,www,adv_uu,adv_vv,iblock)
!
      real(r8) ,intent(in) :: uuu(imt,jmt,km), vvv(imt,jmt,km),www(imt,jmt,km)
      real(r8), intent(out):: adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)
      real(r8) :: u_wface(imt,jmt,km), v_sface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2,adv_z3,adv_z4,adv_xx, adv_yy
      integer, intent(in)  :: iblock
!
     adv_uu = c0
     adv_vv = c0
!
      do k=1, km
      do j= 2,jmt-1
      do i= 2,imt-1
         if ( trim(adv_momentum) == 'centered' ) then
            u_wface(i,j,k) = (uuu(i-1,j,k) + uuu(i,j,k))*P25*hue(i-1,j,iblock)
            v_sface(i,j,k) = (vvv(i,j,k) + vvv(i,j+1,k))*P25*hun(i,j+1,iblock)
         else if ( trim(adv_momentum) == 'flux' ) then
            u_wface(i,j,k) = (uuu(i-1,j,k)*dyu(i-1,j,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P25
            v_sface(i,j,k) = (vvv(i,j,k)*dxu(i,j,iblock) + vvv(i,j+1,k)*dxu(i,j+1,iblock))*P25
         end if
      end do
      end do
      end do
!
      if ( trim(adv_momentum) == 'centered' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_uu(i,j,k) = (-u_wface(i  ,j,k)*(uuu(i  ,j,k)-uuu(i-1,j,k))                   &
                            -u_wface(i+1,j,k)*(uuu(i+1,j,k)-uuu(i  ,j,k))                   &
                            -v_sface(i,j  ,k)*(uuu(i,j+1,k)-uuu(i,j  ,k))                   &
                            -v_sface(i,j-1,k)*(uuu(i,j  ,k)-uuu(i,j-1,k)))*uarea_r(i,j,iblock)
           adv_vv(i,j,k) = (-u_wface(i  ,j,k)*(vvv(i  ,j,k)-vvv(i-1,j,k))                   &
                            -u_wface(i+1,j,k)*(vvv(i+1,j,k)-vvv(i  ,j,k))                   &
                            -v_sface(i,j  ,k)*(vvv(i,j+1,k)-vvv(i,j  ,k))                   &
                            -v_sface(i,j-1,k)*(vvv(i,j  ,k)-vvv(i,j-1,k)))*uarea_r(i,j,iblock)
!
           if (k==1 )then
               adv_z1=0.0D0
               adv_z3=0.0D0
           else
               adv_z1=www (I,J,K)* (uuu(I,J,K -1) - uuu(I,J,K))
               adv_z3=www (I,J,K)* (vvv(I,J,K -1) - vvv(I,J,K))
           end if
!
           if (k==km )then
                adv_z2=0.0D0
                adv_z4=0.0D0
           else
                adv_z2=www(I,J,K+1)*(uuu(I,J,K ) - uuu(I,J,K+1))
                adv_z4=www(I,J,K+1)*(vvv(I,J,K ) - vvv(I,J,K+1))
           end if
           adv_uu(i,j,k) = adv_uu(i,j,k) - P5*ODZP(K)* (adv_z1+adv_z2)
           adv_vv(i,j,k) = adv_vv(i,j,k) - P5*ODZP(K)* (adv_z3+adv_z4)
        end do
        end do
        end do
      else if ( trim(adv_momentum) == 'flux' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_uu(i,j,k) = (-u_wface(i  ,j,k)*(uuu(i  ,j,k)+uuu(i-1,j,k))    &
                            +u_wface(i+1,j,k)*(uuu(i+1,j,k)+uuu(i  ,j,k))    &
                            -v_sface(i,j-1,k)*(uuu(i,  j,k)+uuu(i,j-1,k))    &
                            +v_sface(i,j  ,k)*(uuu(i,j+1,k)+uuu(i,  j,k)))*uarea_r(i,j,iblock)
!
           adv_vv(i,j,k) = (-u_wface(i  ,j,k)*(vvv(i,j  ,k)+vvv(i-1,j,k))    &
                            +u_wface(i+1,j,k)*(vvv(i+1,j,k)+vvv(i  ,j,k))    &
                            -v_sface(i,j-1,k)*(vvv(i,  j,k)+vvv(i,j-1,k))    &
                            +v_sface(i,j  ,k)*(vvv(i,j+1,k)+vvv(i,  j,k)))*uarea_r(i,j,iblock)
!
           if (k ==1 ) then
               adv_z1=0.0D0
               adv_z3=0.0D0
           else
               adv_z1=www (I,J,K)* (uuu(I,J,K -1) + uuu(I,J,K))*P5
               adv_z3=www (I,J,K)* (vvv(I,J,K -1) + vvv(I,J,K))*P5
           end if
!
           if (k == km ) then
                adv_z2=0.0D0
                adv_z4=0.0D0
           else
                adv_z2=www(I,J,K+1)*(uuu(I,J,K ) + uuu(I,J,K+1))*P5
                adv_z4=www(I,J,K+1)*(vvv(I,J,K ) + vvv(I,J,K+1))*P5
           end if
           adv_uu(i,j,k) = adv_uu(i,j,k) - ODZP(K)* (adv_z2-adv_z1)
           adv_vv(i,j,k) = adv_vv(i,j,k) - ODZP(K)* (adv_z4-adv_z3)
        end do
        end do
        end do
      else
        write(6,*) "adv_momentum =", adv_momentum
        call exit_licom(sigAbort,'The false advection option for momentum')
      end if
!
      end subroutine advection_momentum
!
!
!
      subroutine advection_tracer(uuu,vvv,www,ttt,adv_tt,iblock)
!
      real(r8) ,intent(in) :: uuu(imt,jmt,km), vvv(imt,jmt,km),www(imt,jmt,km),ttt(imt,jmt,km)
      integer, intent(in)  :: iblock
      real(r8), intent(out):: adv_tt(imt,jmt,km)
      real(r8) :: u_wface(imt,jmt,km),v_sface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2
!
     adv_tt=0.0_r8
!
      do k=1, km
      do j= 2,jmt-1
      do i= 2,imt-1
         if ( trim(adv_momentum) == 'centered' ) then
            u_wface(i,j,k) = (uuu(i,j-1,k) + uuu(i,j,k))*htw(i,j,iblock)*P25
            v_sface(i,j,k) = (vvv(i,j,k) + vvv(i+1,j,k))*hts(i,j,iblock)*P25
         else if ( trim(adv_momentum) == 'flux' ) then
            u_wface(i,j,k) = (uuu(i,j-1,k)*dyu(i,j-1,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P25
            v_sface(i,j,k) = (vvv(i,j,k)*dxu(i,j,iblock) + vvv(i+1,j,k)*dxu(i+1,j,iblock))*P25
         end if
      end do
      end do
      end do
!
      if ( trim(adv_tracer) == 'centered' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_tt(i,j,k) = (-u_wface(i  ,j,k)*(ttt(i  ,j,k)-ttt(i-1,j,k))                      &
                            -u_wface(i+1,j,k)*(ttt(i+1,j,k)-ttt(i  ,j,k))                      &
                            -v_sface(i,j  ,k)*(ttt(i,j+1,k)-ttt(i,j  ,k))                      &
                            -v_sface(i,j-1,k)*(ttt(i,j  ,k)-ttt(i,j-1,k)))*tarea_r(i,j,iblock)
!
           if (k==1 )then
               adv_z1=0.0D0
           else
               adv_z1=www (I,J,K)* (ttt(I,J,K -1) - ttt(I,J,K))
           end if
!
           if (k==km )then
                adv_z2=0.0D0
           else
                adv_z2=www(I,J,K+1)*(ttt(I,J,K ) - ttt(I,J,K+1))
           end if
!
           adv_tt(i,j,k) = adv_tt(i,j,k) - P5*ODZP(K)* (adv_z1+adv_z2)
        end do
        end do
        end do
      else if ( trim(adv_tracer) == 'flux' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_tt(i,j,k) = (-u_wface(i  ,j,k)*(ttt(i  ,j,k)+ttt(i-1,j,k))   &
                            +u_wface(i+1,j,k)*(ttt(i+1,j,k)+ttt(i  ,j,k))   &
                            -v_sface(i,j-1,k)*(ttt(i,  j,k)+ttt(i,j-1,k))   &
                            +v_sface(i,j  ,k)*(ttt(i,j+1,k)+ttt(i,  j,k)))*tarea_r(i,j,iblock)
!
           if (k==1 )then
               adv_z1=0.0D0
           else
               adv_z1=www (I,J,K)* (ttt(I,J,K -1) + ttt(I,J,K))*P5
           end if
!
           if (k==km )then
                adv_z2=0.0D0
           else
                adv_z2=www(I,J,K+1)*(ttt(I,J,K ) + ttt(I,J,K+1))*P5
           end if
!
           adv_tt(i,j,k) = adv_tt(i,j,k) - ODZP(K)* (adv_z2-adv_z1)
!
        end do
        end do
        end do
      else
        call exit_licom(sigAbort,'The false advection option for tracer')
      end if
!

      end subroutine advection_tracer

!***********************************************************************

 end module advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
