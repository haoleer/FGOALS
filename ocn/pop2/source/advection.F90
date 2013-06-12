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

      real(r8)    :: LAMDA,wt1,wt2,adv_z
      real(r8),dimension(:,:,:),allocatable :: adv_xy1,adv_xy2,adv_xy3,adv_xy4, at00
      real(r8),dimension(:,:,:),allocatable :: adv_zz,atz, adv_za,adv_zb1,adv_zb2,adv_zc,atmaxz,atminz
      real(r8),dimension(:,:,:),allocatable :: adv_x0,adv_y0,adv_c1,adv_c2,atmax,atmin,adv_xx,adv_yy
!
   public :: advection_momentum, advection_tracer

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
      real(r8) :: adv_z1, adv_z2,adv_z3,adv_z4
      integer, intent(in)  :: iblock
!
     adv_uu = c0
     adv_vv = c0
!
      do k=1, km
      do j= 1,jmt-1
      do i= 2,imt
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
      real(r8) ,intent(in) :: uuu(imt,jmt,km), vvv(imt,jmt,km),www(imt,jmt,km+1),ttt(imt,jmt,km)
      integer, intent(in)  :: iblock
      real(r8), intent(out):: adv_tt(imt,jmt,km)
      real(r8) :: u_wface(imt,jmt,km),v_sface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2
!
     adv_tt=0.0_r8
!
      do k=1, km
      do j=1, jmt
      do i=2, imt-1
         if ( trim(adv_tracer) == 'centered' .or. trim(adv_tracer) == 'tspas') then
            v_sface(i,j,k) = (vvv(i,j,k) + vvv(i+1,j,k))*hts(i,j,iblock)*P25
         else if ( trim(adv_tracer) == 'flux' ) then
            v_sface(i,j,k) = (vvv(i,j,k)*dxu(i,j,iblock) + vvv(i+1,j,k)*dxu(i+1,j,iblock))*P25
         end if
      end do
      end do
      do j= 2,jmt-1
      do i= 1,imt
         if ( trim(adv_tracer) == 'centered' .or. trim(adv_tracer) == 'tspas') then
            u_wface(i,j,k) = (uuu(i,j-1,k) + uuu(i,j,k))*htw(i,j,iblock)*P25
         else if ( trim(adv_tracer) == 'flux' ) then
            u_wface(i,j,k) = (uuu(i,j-1,k)*dyu(i,j-1,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P25
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
      else if (trim(adv_tracer) == 'tspas' ) then
!
      allocate ( adv_x0(imt,jmt,km),adv_y0(imt,jmt,km),adv_c1(imt,jmt,km) )
      allocate ( adv_c2(imt,jmt,km), adv_xx(imt,jmt,km), adv_yy(imt,jmt,km))
      allocate ( adv_xy1(imt,jmt,km), adv_xy2(imt,jmt,km), adv_xy3(imt,jmt,km), adv_xy4(imt,jmt,km))
      allocate ( at00(imt,jmt,km), atmax(imt,jmt,km), atmin(imt,jmt,km) )
!
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = 2, jmt-1
          DO I = 2, imt-1
            adv_x0(i,j,k)=((ttt(i+1,j,k)+ttt(i,j,k))*u_wface(i+1,j,k) &
                       -(ttt(i,j,k)+ttt(i-1,j,k))*u_wface(i,j,k))*tarea_r(i,j,iblock)
            adv_y0(i,j,k)=((ttt(i,j+1,k)+ttt(i,j,k))*v_sface(i,j,k) &
                       -(ttt(i,j,k)+ttt(i,j-1,k))*v_sface(i,j-1,k))*tarea_r(i,j,iblock)
            adv_xy1(i,j,k)=-dts*(ttt(i+1,j,k)-ttt(i,j,k))*2.0_r8* &
                                  u_wface(i+1,j,k)*u_wface(i+1,j,k)*tarea_r(i,j,iblock)*tarea_r(i,j,iblock)
            adv_xy2(i,j,k)= dts*(ttt(i,j,k)-ttt(i-1,j,k))*2.0_r8* &
                                  u_wface(i,j,k)*u_wface(i,j,k)*tarea_r(i,j,iblock)*tarea_r(i,j,iblock)
            adv_xy3(i,j,k)=-dts*(ttt(i,j+1,k)-ttt(i,j,k))*2.0_r8*dytr(i,j,iblock)*dyur(i,j,iblock)* &
                                  v_sface(i,j,k)*v_sface(i,j,k)*dxtr(i,j,iblock)*dxtr(i,j,iblock)
            adv_xy4(i,j,k)= dts*(ttt(i,j,k)-ttt(i,j-1,k))*2.0_r8*dytr(i,j,iblock)*dyur(i,j-1,iblock)* &
                                  v_sface(i,j-1,k)*v_sface(i,j-1,k)*dxtr(i,j,iblock)*dxtr(i,j,iblock)
            adv_c1(i,j,k)=-ttt(i,j,k)*(u_wface(i+1,j,k)-u_wface(i,j,k))*tarea_r(i,j,iblock)*2.0_r8
            adv_c2(i,j,k)=-ttt(i,j,k)*(v_sface(i,j,k)-v_sface(i,j-1,k))*tarea_r(i,j,iblock)*2.0_r8

            adv_xx(i,j,k)=-(adv_x0(i,j,k)+adv_xy1(i,j,k)+adv_xy2(i,j,k)+adv_c1(i,j,k))
            adv_yy(i,j,k)=-(adv_y0(i,j,k)+adv_xy3(i,j,k)+adv_xy4(i,j,k)+adv_c2(i,j,k))

            at00(i,j,k)=ttt(i,j,k)+(adv_xx(i,j,k)+adv_yy(i,j,k))*dts
            if (mytid == 5 .and. k==1 .and. i == 35  .and. j ==25 ) then
               write(128,*) adv_x0(i,j,k), adv_y0(i,j,k), adv_xy1(i,j,k), adv_xy2(i,j,k)
               write(128,*) adv_xy3(i,j,k), adv_xy4(i,j,k), adv_c1(i,j,k), adv_c2(i,j,k)
               write(128,*) adv_xx(i,j,k), adv_yy(i,j,k), at00(i,j,k)
            end if
            if (mytid == 0 .and. k==2 .and. i > 2 .and. i < imt-1 .and. j > 5 .and. j< 9) then
               write(124,*) i,j,k
               write(124,*) adv_x0(i,j,k), adv_y0(i,j,k), adv_xy1(i,j,k), adv_xy2(i,j,k)
               write(124,*) adv_xy3(i,j,k), adv_xy4(i,j,k), adv_c1(i,j,k), adv_c2(i,j,k)
               write(124,*) adv_xx(i,j,k), adv_yy(i,j,k), at00(i,j,k)
            end if
          ENDDO
        ENDDO
      ENDDO
      if(mytid == 0) close(124)
      if(mytid == 0) close(128)

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = 2, jmt-1
          DO I = 2,imt-1
            atmax(i,j,k)=max(ttt(i,j,k),ttt(i,j-1,k),ttt(i,j+1,k),ttt(i-1,j,k),ttt(i+1,j,k))
            atmin(i,j,k)=min(ttt(i,j,k),ttt(i,j-1,k),ttt(i,j+1,k),ttt(i-1,j,k),ttt(i+1,j,k))
          ENDDO
        ENDDO
      ENDDO

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = 3, jmt-2
          DO I = 3, imt-2
            if (at00(i,j,k)>atmax(i,j,k).or.at00(i,j,k)<atmin(i,j,k)) then
              adv_xy1(i,j,k)=-(ttt(i+1,j,k)-ttt(i,j,k))*        &
                                 abs(u_wface(i+1,j,k))*tarea_r(i,j,iblock)
              adv_xy2(i,j,k)= (ttt(i,j,k)-ttt(i-1,j,k))*        &
                                 abs(u_wface(i,j,k))*tarea_r(i,j,iblock)
              adv_xy3(i,j,k)=-(ttt(i,j+1,k)-ttt(i,j,k))*        &
                                 abs(v_sface(i,j,k))*tarea_r(i,j,iblock)
              adv_xy4(i,j,k)= (ttt(i,j,k)-ttt(i,j-1,k))*        &
                                 abs(v_sface(i,j-1,k))*tarea_r(i,j,iblock)
            else
              if (at00(i+1,j,k)>atmax(i+1,j,k).or.at00(i+1,j,k)<atmin(i+1,j,k)) then
                 adv_xy1(i,j,k)=-(ttt(i+1,j,k)-ttt(i,j,k))*        &
                                 abs(u_wface(i+1,j,k))*tarea_r(i,j,iblock)
              endif
              if (at00(i-1,j,k)>atmax(i-1,j,k).or.at00(i-1,j,k)<atmin(i-1,j,k)) then
                 adv_xy2(i,j,k)= (ttt(i,j,k)-ttt(i-1,j,k))*        &
                                 abs(u_wface(i,j,k))*tarea_r(i,j,iblock)
              endif
              if (at00(i,j+1,k)>atmax(i,j+1,k).or.at00(i,j+1,k)<atmin(i,j+1,k)) then
                 adv_xy3(i,j,k)=-(ttt(i,j+1,k)-ttt(i,j,k))*        &
                                 abs(v_sface(i,j,k))*tarea_r(i,j,iblock)
              endif
              if (at00(i,j-1,k)>atmax(i,j-1,k).or.at00(i,j-1,k)<atmin(i,j-1,k)) then
                 adv_xy4(i,j,k)= (ttt(i,j,k)-ttt(i,j-1,k))*        &
                                 abs(v_sface(i,j-1,k))*tarea_r(i,j,iblock)
              endif
            endif

            adv_xx(i,j,k)=-(adv_x0(i,j,k)+adv_xy1(i,j,k)+adv_xy2(i,j,k)+adv_c1(i,j,k))
            adv_yy(i,j,k)=-(adv_y0(i,j,k)+adv_xy3(i,j,k)+adv_xy4(i,j,k)+adv_c2(i,j,k))
            adv_tt (I,J,K)= adv_xx(i,j,k)+adv_yy(i,j,k)
            if (mytid == 0 .and. k==2 .and. i > 2 .and. i < imt-1 .and. j > 5 .and. j< 9) then
               write(125,*) i,j,k
               write(125,*) adv_x0(i,j,k), adv_y0(i,j,k), adv_xy1(i,j,k), adv_xy2(i,j,k)
               write(125,*) adv_xy3(i,j,k), adv_xy4(i,j,k), adv_c1(i,j,k), adv_c2(i,j,k)
               write(125,*) adv_xx(i,j,k), adv_yy(i,j,k), at00(i,j,k), adv_tt(i,j,k)
               write(125,*) atmax(i,j,k), atmin(i,j,k)
            end if
            if (mytid == 5 .and. k==1 .and. i == 35  .and. j ==25 ) then
               write(129,*) adv_x0(i,j,k), adv_y0(i,j,k), adv_xy1(i,j,k), adv_xy2(i,j,k)
               write(129,*) adv_xy3(i,j,k), adv_xy4(i,j,k), adv_c1(i,j,k), adv_c2(i,j,k)
               write(129,*) adv_xx(i,j,k), adv_yy(i,j,k), at00(i,j,k), adv_tt(i,j,k)
               write(129,*) atmax(i,j,k), atmin(i,j,k)
            end if
          END DO
        END DO
      END DO
      if (mytid ==0) close(125)
      if (mytid ==0) close(129)

      allocate (adv_zz(imt,jmt,km), adv_za(imt,jmt,km),adv_zb1(imt,jmt,km) )
      allocate (adv_zb2(imt,jmt,km), adv_zc(imt,jmt,km), atmaxz(imt,jmt,km), atminz(imt,jmt,km), atz(imt,jmt,km))
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = 3, jmt-2
          DO I = 3, imt-2
            if (k==1) then
              adv_za (i,j,k)=-0.5*ODZP(1)*WWW(I,J,2)*(ttt(I,J,2)+ttt(I,J,1))
              adv_zb1(i,j,k)=0
              adv_zb2(i,j,k)= 0.5*ODZP(1)*WWW(I,J,2)*WWW(I,J,2)*ODZT(2) &
                                *(ttt(I,J,1)-ttt(I,J,2))
              adv_zc (i,j,k)=     ODZP(1)*ttt(I,J,1)*WWW(I,J,2)
              atmaxz(i,j,k)=max(ttt(i,j,1),ttt(i,j,2))
              atminz(i,j,k)=min(ttt(i,j,1),ttt(i,j,2))
              atz(i,j,k)=ttt(i,j,k)-(adv_za(i,j,k)+adv_zb1(i,j,k) &
                                     +adv_zb2(i,j,k)+adv_zc(i,j,k))*dts
            else if (k==km) then
              adv_za (i,j,k)= 0.5*ODZP(km)*WWW(I,J,km)*(ttt(I,J,km)+ttt(I,J,km-1))
              adv_zb1(i,j,k)=-0.5*ODZP(km)*WWW(I,J,km)*WWW(I,J,km)*ODZT(km  ) &
                                *(ttt(I,J,km-1)-ttt(I,J,km))
              adv_zb2(i,j,k)=0
              adv_zc (i,j,k)=    -ODZP(km)*ttt(I,J,km)*WWW(I,J,km)
              atmaxz(i,j,k)=max(ttt(i,j,km-1),ttt(i,j,km))
              atminz(i,j,k)=min(ttt(i,j,km-1),ttt(i,j,km))
              atz(i,j,k)=ttt(i,j,k)-(adv_za(i,j,k)+adv_zb1(i,j,k) &
                                  +adv_zb2(i,j,k)+adv_zc(i,j,k))*dts
            else
              adv_za (i,j,k)= 0.5*ODZP(k)*WWW(I,J,k  )*(ttt(I,J,k)+ttt(I,J,k-1)) &
                         -0.5*ODZP(k)*WWW(I,J,k+1)*(ttt(I,J,k)+ttt(I,J,k+1))
              adv_zb1(i,j,k)=-0.5*ODZP(k)*WWW(I,J,k  )*WWW(I,J,k  )*ODZT(k  ) &
                                *(ttt(I,J,k-1)-ttt(I,J,k))
              adv_zb2(i,j,k)= 0.5*ODZP(k)*WWW(I,J,k+1)*WWW(I,J,k+1)*ODZT(k+1) &
                                *(ttt(I,J,k)-ttt(I,J,k+1))
              adv_zc (i,j,k)=    -ODZP(k)*ttt(I,J,k)*(WWW(I,J,k)-WWW(I,J,k+1))
              atmaxz (i,j,k)=max(ttt(i,j,k-1),ttt(i,j,k),ttt(i,j,k+1))
              atminz (i,j,k)=min(ttt(i,j,k-1),ttt(i,j,k),ttt(i,j,k+1))
              atz(i,j,k)=ttt(i,j,k)-(adv_za(i,j,k)+adv_zb1(i,j,k) &
                                   +adv_zb2(i,j,k)+adv_zc(i,j,k))*dts
            endif
            if (mytid == 5 .and. k==1 .and. i == 35  .and. j ==25 ) then
!              write(130,*) www(i,j,k),www(i,j,k-1),www(i,j,k+1)
!              write(130,*) ttt(i,j,k),ttt(i,j,k-1),ttt(i,j,k+1)
               write(130,*) www(i,j,k),www(i,j,k+1)
               write(130,*) ttt(i,j,k),ttt(i,j,k+1)
               write(130,*) odzp(1), odzt(2), www(i,j,2), ttt(i,j,1)-ttt(i,j,2)
               write(130,*) adv_za(i,j,k), adv_zb1(i,j,k), adv_zb2(i,j,k)
               write(130,*) adv_zc(i,j,k), atmaxz(i,j,k), atminz(i,j,k), atz(i,j,k)
            end if
            if (mytid == 0 .and. k==2 .and. i > 2 .and. i < imt-1 .and. j > 5 .and. j< 9) then
               write(126,*) i,j,k
               write(126,*) odzp(k), www(i,j,k),www(i,j,k-1),www(i,j,k+1)
               write(126,*) ttt(i,j,k),ttt(i,j,k-1),ttt(i,j,k+1)
               write(126,*) adv_za(i,j,k), adv_zb1(i,j,k), adv_zb2(i,j,k)
               write(126,*) adv_zc(i,j,k), atmaxz(i,j,k), atminz(i,j,k), atz(i,j,k)
               write(126,*) adv_xx(i,j,k), adv_yy(i,j,k), at00(i,j,k)
            end if
          END DO
        END DO
      END DO
      if (mytid ==0) close(126)
      if (mytid ==0) close(130)

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,km
        DO J = 3, jmt-2
          DO I = 3, imt-2
            if (k==1) then
              if(atz(i,j,k+1)>atmaxz(i,j,k+1).or.atz(i,j,k+1)<atminz(i,j,k+1).or. &
                 atz(i,j,k)>atmaxz(i,j,k).or.atz(i,j,k)<atminz(i,j,k)) then
                adv_zb2(i,j,k)= 0.5*abs(WWW(I,J,k+1))*ODZT(k+1)*(ttt(I,J,k)-ttt(I,J,k+1))
              endif
             else if (k==km) then
               if(atz(i,j,k-1)>atmaxz(i,j,k-1).or.atz(i,j,k-1)<atminz(i,j,k-1).or. &
                  atz(i,j,k  )>atmaxz(i,j,k  ).or.atz(i,j,k  )<atminz(i,j,k  )) then
                 adv_zb1(i,j,k)=-0.5*abs(WWW(I,J,k))*ODZT(k)*(ttt(I,J,k-1)-ttt(I,J,k))
               endif
             else
               if(atz(i,j,k)>atmaxz(i,j,k).or.atz(i,j,k)<atminz(i,j,k)) then
                 adv_zb1(i,j,k)=-0.5*abs(WWW(I,J,k  ))*ODZT(k)*(ttt(I,J,k-1)-ttt(I,J,k))
                 adv_zb2(i,j,k)= 0.5*abs(WWW(I,J,k+1))*ODZT(k+1)*(ttt(I,J,k)-ttt(I,J,k+1))
               else
                 if(atz(i,j,k+1)>atmaxz(i,j,k+1).or.atz(i,j,k+1)<atminz(i,j,k+1)) then
                   adv_zb2(i,j,k)= 0.5*abs(WWW(I,J,k+1))*ODZT(k+1)*(ttt(I,J,k)-ttt(I,J,k+1))
                 endif
                 if(atz(i,j,k-1)>atmaxz(i,j,k-1).or.atz(i,j,k-1)<atminz(i,j,k-1)) then
                   adv_zb1(i,j,k)=-0.5*abs(WWW(I,J,k))*ODZT(k)*(ttt(I,J,k-1)-ttt(I,J,k))
                 endif
               endif
             endif
!
             adv_zz(i,j,k)=-(adv_za(i,j,k)+adv_zb1(i,j,k)+adv_zb2(i,j,k)+adv_zc(i,j,k))
             atz(i,j,k)=ttt(i,j,k)+adv_zz(i,j,k)*dts
             adv_tt(I,J,K)= adv_tt(I,J,K)+adv_zz(i,j,k)
!
            if (mytid == 5 .and. k==1 .and. i == 35  .and. j ==25 ) then
               write(131,*) adv_za(i,j,k), adv_zb1(i,j,k), adv_zb2(i,j,k)
               write(131,*) adv_zc(i,j,k), atmaxz(i,j,k), atminz(i,j,k), atz(i,j,k)
               write(131,*) adv_zz(i,j,k), adv_tt(i,j,k)
               write(131,*) atz(i,j,k-1),atmaxz(i,j,k-1), atminz(i,j,k-1)
               write(131,*) atz(i,j,k+1),atmaxz(i,j,k+1), atminz(i,j,k+1)
            end if
            if (mytid == 0 .and. k==2 .and. i > 2 .and. i < imt-1 .and. j > 5 .and. j< 9) then
               write(127,*) i,j,k
               write(127,*) adv_za(i,j,k), adv_zb1(i,j,k), adv_zb2(i,j,k)
               write(127,*) adv_zc(i,j,k), atmaxz(i,j,k), atminz(i,j,k), atz(i,j,k)
               write(127,*) adv_zz(i,j,k), adv_tt(i,j,k)
               write(127,*) atz(i,j,k-1),atmaxz(i,j,k-1), atminz(i,j,k-1)
               write(127,*) atz(i,j,k+1),atmaxz(i,j,k+1), atminz(i,j,k+1)
            end if

          END DO
        END DO
      END DO
!
     if (mytid ==0) close(127)
     if (mytid ==0) close(131)
!
      deallocate ( adv_x0, adv_y0, adv_c1)
      deallocate ( adv_c2, adv_xx, adv_yy)
      deallocate ( adv_xy1, adv_xy2, adv_xy3, adv_xy4)
      deallocate ( at00, atmax, atmin )
      deallocate (adv_zz, adv_za, adv_zb1 )
      deallocate (adv_zb2, adv_zc, atmaxz, atminz, atz)
!
      else
        call exit_licom(sigAbort,'The false advection option for tracer')
      end if
!

      end subroutine advection_tracer

!***********************************************************************

 end module advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
