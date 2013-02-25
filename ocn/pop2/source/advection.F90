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
   use LICOM_ErrorSet

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public :: advection_momentum, &
             advection_tracer

!EOP
!BOC

!-----------------------------------------------------------------------
!EOC
      subroutine advection_momentum(uuu,vvv,www,adv_uu,adv_vv,iblock)
!
      real(r8) ,intent(in) :: uuu(imt,jmt,km), vvv(imt,jmt,km),www(imt,jmt,km)
      real(r8), intent(out):: adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)
      real(r8) :: u_eface(imt,jmt,km), v_nface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2,adv_z3,adv_z4
      integer, intent(in)  :: iblock
!
     adv_uu = c0
     adv_vv = c0
!
      do k=1, km
      do j= 2,jmt
      do i= 2,imt
         if ( trim(adv_momentum) == 'centered' ) then
            u_eface(i,j,k) = (uuu(i-1,j,k) + uuu(i,j,k))*P5
            v_nface(i,j,k) = (vvv(i,j-1,k) + vvv(i,j,k))*P5
         else if ( trim(adv_momentum) == 'flux' ) then
            u_eface(i,j,k) = (uuu(i-1,j,k)*dyu(i-1,j,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P5
            v_nface(i,j,k) = (vvv(i,j-1,k)*dxu(i,j-1,iblock) + vvv(i,j,k)*dxu(i,j,iblock))*P5
         end if
      end do
      end do
      end do
!
      if ( trim(adv_momentum) == 'centered' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_uu(i,j,k) =  -u_eface(i  ,j,k)*(uuu(i  ,j,k)-uuu(i-1,j,k))*r1a_u(i,j,k)  &
                            -u_eface(i+1,j,k)*(uuu(i+1,j,k)-uuu(i  ,j,k))*r1b_u(i,j,k)  &
                            -v_nface(i,j  ,k)*(uuu(i,  j,k)-uuu(i,j-1,k))*r2a_u(i,j,k)  &
                            -v_nface(i,j+1,k)*(uuu(i,j+1,k)-uuu(i,  j,k))*r2b_u(i,j,k)  
           adv_vv(i,j,k) =  -u_eface(i  ,j,k)*(vvv(i  ,j,k)-vvv(i-1,j,k))*r1a_u(i,j,k)  &
                            -u_eface(i+1,j,k)*(vvv(i+1,j,k)-vvv(i  ,j,k))*r1b_u(i,j,k)  &
                            -v_nface(i,j  ,k)*(vvv(i,  j,k)-vvv(i,j-1,k))*r2a_u(i,j,k)  &
                            -v_nface(i,j+1,k)*(vvv(i,j+1,k)-vvv(i,  j,k))*r2b_u(i,j,k)  
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
           adv_uu(i,j,k) = (-u_eface(i  ,j,k)*(uuu(i  ,j,k)+uuu(i-1,j,k))    &
                            +u_eface(i+1,j,k)*(uuu(i+1,j,k)+uuu(i  ,j,k))    &
                            -v_nface(i,j  ,k)*(uuu(i,  j,k)+uuu(i,j-1,k))    &
                            +v_nface(i,j+1,k)*(uuu(i,j+1,k)+uuu(i,  j,k)))*uarea_r(i,j,iblock)*P5 &
           adv_vv(i,j,k) = (-u_eface(i  ,j,k)*(vvv(i  ,j,k)+vvv(i-1,j,k))    &
                            +u_eface(i+1,j,k)*(vvv(i+1,j,k)+vvv(i  ,j,k))    &
                            -v_nface(i,j  ,k)*(vvv(i,  j,k)+vvv(i,j-1,k))    &
                            +v_nface(i,j+1,k)*(vvv(i,j+1,k)+vvv(i,  j,k)))*uarea_r(i,j,iblock)*P5 
!
           if (k==1 )then
               adv_z1=0.0D0
               adv_z3=0.0D0
           else
               adv_z1=www (I,J,K)* (uuu(I,J,K -1) + uuu(I,J,K))*P5
               adv_z3=www (I,J,K)* (vvv(I,J,K -1) + vvv(I,J,K))*P5
           end if
!
           if (k==km )then
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
         call exit_LICOM
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
      real(r8) :: u_eface(imt,jmt,km),v_nface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2
!
     adv_tt=0.0_r8
!
      do k=1, km
      do j= 2,jmt
      do i= 1,imt-1
         if ( trim(adv_momentum) == 'centered' ) then
            u_eface(i,j,k) = (uuu(i,j-1,k) + uuu(i,j,k))*P5
            v_nface(i,j,k) = (vvv(i,j-1,k) + vvv(i+1,j-1,k))*P5
         else if ( trim(adv_momentum) == 'flux' ) then
            u_eface(i,j,k) = (uuu(i,j-1,k)*dyu(i,j-1,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P5
            v_nface(i,j,k) = (vvv(i,j-1,k)*dxu(i,j-1,iblock) + vvv(i+1,j-1,k)*dxu(i+1,j-1,iblock))*P5
         end if
      end do
      end do
      end do
!
      if ( trim(adv_tracer) == 'centered' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_tt(i,j,k) =  -u_eface(i  ,j,k)*(ttt(i  ,j,k)-ttt(i-1,j,k))*r1a_t(i,j,k)  &
                            -u_eface(i+1,j,k)*(ttt(i+1,j,k)-ttt(i  ,j,k))*r1b_t(i,j,k)  &
                            -v_nface(i,j  ,k)*(ttt(i,  j,k)-ttt(i,j-1,k))*r2a_t(i,j,k)  &
                            -v_nface(i,j+1,k)*(ttt(i,j+1,k)-ttt(i,  j,k))*r2b_t(i,j,k)
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
           adv_tt(i,j,k) = (-u_eface(i  ,j,k)*(ttt(i  ,j,k)+ttt(i-1,j,k))   &
                            +u_eface(i+1,j,k)*(ttt(i+1,j,k)+ttt(i  ,j,k))   &
                            -v_nface(i,j  ,k)*(ttt(i,  j,k)+ttt(i,j-1,k))   &
                            +v_nface(i,j+1,k)*(ttt(i,j+1,k)+ttt(i,  j,k)))*P5*tarea_r(i,j,iblock)
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
        end do
        end do
        end do
      else
         call exit_LICOM
      end if

      end subroutine advection_tracer

!***********************************************************************

 end module advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
