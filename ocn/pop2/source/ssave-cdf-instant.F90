!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEINS(EClock)
#include <def-undef.h>
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
!LPF 20120815
use buf_mod, only:t_cpl,s_cpl,u_cpl,v_cpl,dhdx,dhdy,q      
!LPF 20120815
use domain
use gather_scatter
use distribution
use mct_mod
use esmf_mod
use seq_flds_mod
use seq_cdata_mod
use seq_infodata_mod
use seq_timemgr_mod

      logical       :: write_restart     ! restart now
      character (len=18) :: fname
      integer :: klevel
      type(ESMF_Clock)            ,intent(in)    :: EClock
      integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD
!
!        call seq_timemgr_EClockGetData( EClock, curr_ymd=curr_ymd)
         write_restart = seq_timemgr_RestartAlarmIsOn(EClock)
         number_day = iday !+1 !LPF 20120816
         number_month = month
         if ( number_day  > imd ) then
             number_day = 1
             number_month = month + 1
         end if
         nwmf= iyfm

!
     if(mytid==0) write(*,*)'in instant,iday=,imd=,rest_freq',iday,imd,rest_freq
     if(mytid==0) write(*,*)'in instant,write_restart', write_restart
!
!    if ( mod(iday,rest_freq) == 0 .or. iday == imd .or. iday == 10 .or. iday ==20) then
     if ( mod(iday,rest_freq) == 0 .or. iday == imd ) then
!
       if ( mytid == 0 ) then
         fname(1:8)='fort.22.'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon0
         write(fname(9:12),'(i4.4)')nwmf
         write(fname(17:18),'(i2.2)') number_day
 
         if ( number_day ==1 .and. mon0 < 12) then
             write(fname(14:15),'(i2.2)')mon0+1
         end if
 
         if ( number_day ==1 .and. mon0 == 12) then
             write(fname(14:15),'(i2.2)')mon0-11
             write(fname(9:12),'(i4.4)')nwmf+1
         end if

         open (17, file="rpointer.ocn", form='formatted')
         write(17,'(a18)') fname
         close(17)
         write(*,*)'fname=',fname
         open(22,file=trim(out_dir)//fname,form='unformatted')
       end if
!
       allocate(buffer(imt_global,jmt_global))
       call gather_global(buffer, h0, master_task,distrb_clinic)
!
       if (mytid==0) then
          WRITE (22)buffer
       end if
!
      do klevel=1,km
          call gather_global(buffer, u(:,:,klevel,:), master_task,distrb_clinic)
          if (mytid==0) then
             WRITE (22)buffer
          end if
       end do
!         write(*,*)'finish U'
!        
      do klevel=1,km
          call gather_global(buffer, v(:,:,klevel,:), master_task,distrb_clinic)
          if (mytid==0) then
             WRITE (22)buffer
          end if
       end do
!         write(*,*)'finish V'
!
      do klevel=1,km
         call gather_global(buffer, at(:,:,klevel,1,:), master_task,distrb_clinic)
         if (mytid==0) then
             WRITE (22)buffer
         end if
      end do
!         write(*,*)'finish at1'
!
      do klevel=1,km
         call gather_global(buffer, at(:,:,klevel,2,:), master_task,distrb_clinic)
         if (mytid==0) then
           WRITE (22)buffer
         end if
      end do
!         write(*,*)'finish at2'
!lhl20110728
      do klevel=1,km
         call gather_global(buffer, ws(:,:,klevel,:), master_task,distrb_clinic)
         if (mytid==0) then
             WRITE (22)buffer
         end if
      end do
!         write(*,*)'finish at2'
!        
      call gather_global(buffer, su, master_task,distrb_clinic)
      if (mytid==0) then
         WRITE (22)buffer
      end if
!         write(*,*)'finish su'

      call gather_global(buffer, sv, master_task,distrb_clinic)
      if (mytid==0) then
         WRITE (22)buffer
      end if

!         write(*,*)'finish sv'
     call gather_global(buffer, swv, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish lwv'
     call gather_global(buffer, sshf, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish sshf'
     call gather_global(buffer, lthf, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish lthf'
     call gather_global(buffer, fresh, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish fresh'
     if (mytid==0) then
         write(22) number_month, number_day
     end if
!
#ifdef COUP
     call gather_global(buffer, t_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish t_cpl'
     call gather_global(buffer, s_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish s_cpl'
     call gather_global(buffer, u_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish u_cpl'
     call gather_global(buffer, v_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish v_cpl'
     call gather_global(buffer, dhdx, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish dhdx'
     call gather_global(buffer, dhdy, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish dhdy'
     call gather_global(buffer, q, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish q'
     deallocate( buffer)
#endif

      if(mytid==0) then
        close(22)
      end if
!
  end if

!$OMP PARALLEL DO PRIVATE (iblock,k,j,i)
   do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
            DO i = 1,imt
               up (i,j,k,iblock) = u (i,j,k,iblock)
               vp (i,j,k,iblock) = v (i,j,k,iblock)
               utf (i,j,k,iblock) = u (i,j,k,iblock)
               vtf (i,j,k,iblock) = v (i,j,k,iblock)
               atb (i,j,k,1,iblock) = at (i,j,k,1,iblock)
               atb (i,j,k,2,iblock) = at (i,j,k,2,iblock)
            END DO
         END DO
      END DO
   END DO

         CALL VINTEG (U,UB)
         CALL VINTEG (V,VB)

!$OMP PARALLEL DO PRIVATE (iblock,j,i)
    do iblock = 1, nblocks_clinic
      DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
         DO i = 1,imt
            h0p (i,j,iblock)= h0 (i,j,iblock)
            ubp (i,j,iblock)= ub (i,j,iblock)
            vbp (i,j,iblock)= vb (i,j,iblock)
            h0f (i,j,iblock)= h0 (i,j,iblock)
            h0bf (i,j,iblock)= h0 (i,j,iblock)
         END DO
      END DO
   end do
 
     ISB = 0
     ISC = 0
     IST = 0
 
!     write(*,*)'ok instant'
      return
      end

#if (defined NETCDF) || (defined ALL)
      SUBROUTINE check_err (iret)
#include <netcdf.inc>
      INTEGER :: iret
      IF (iret /= NF_NOERR) THEN
         PRINT *, nf_strerror (iret)
         STOP
      END IF
      END SUBROUTINE check_err
#endif

