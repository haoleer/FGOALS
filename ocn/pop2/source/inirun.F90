!  CVS: $Id: inirun.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE INIRUN
!     =================
!     INITIALIZING FOR ALL PHYSICAL FIELDS
 
#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_msg_mod
use shr_mpi_mod
use shr_sys_mod
use buf_mod
use control_mod
use domain
use constant_mod
use gather_scatter


      IMPLICIT NONE

  !----- local  ------
  integer            :: fid    ! nc domain file ID
  integer            :: dimid  ! nc dimension id
  integer            :: vid    ! nc variable ID
  integer            :: rcode  ! nc return code
  integer            :: ntim   ! temporary
  integer            :: iblock   ! temporary
 
#include <netcdf.inc>
!
!     Define Variables.
      integer*4   :: ncid, iret
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
      character (len=18) :: fname

      INTEGER :: NMFF

      allocate(h0(imt,jmt,max_blocks_clinic),u(imt,jmt,km,max_blocks_clinic),v(imt,jmt,km,max_blocks_clinic), &
               at(imt,jmt,km,ntra,max_blocks_clinic))
      allocate(buffer(imt_global,jmt_global))


      if (mytid==0)then
          write(6,*)"BEGINNING-----INIRUN !"
          open (17,file='rpointer.ocn',form='formatted')
          read(17,'(a18)') fname
          close(17)
      endif 
 
#ifdef COUP
  ! lihuimin, 2012.7.17, nx,ny --> imt,jmt
  allocate (t_cpl(imt,jmt,max_blocks_clinic))
  allocate (s_cpl(imt,jmt,max_blocks_clinic))
  allocate (u_cpl(imt,jmt,max_blocks_clinic))
  allocate (v_cpl(imt,jmt,max_blocks_clinic))
  allocate (dhdx(imt,jmt,max_blocks_clinic))
  allocate (dhdy(imt,jmt,max_blocks_clinic))
  allocate (Q   (imt,jmt,max_blocks_clinic))

  allocate (taux (imt,jmt,max_blocks_clinic))
  allocate (tauy (imt,jmt,max_blocks_clinic))
  allocate (netsw(imt,jmt,max_blocks_clinic))
  allocate (lat1 (imt,jmt,max_blocks_clinic))
  allocate (sen  (imt,jmt,max_blocks_clinic))
  allocate (lwup (imt,jmt,max_blocks_clinic))
  allocate (lwdn (imt,jmt,max_blocks_clinic))
  allocate (melth(imt,jmt,max_blocks_clinic))
  allocate (salt (imt,jmt,max_blocks_clinic))
  allocate (prec (imt,jmt,max_blocks_clinic))
  allocate (evap (imt,jmt,max_blocks_clinic))
  allocate (meltw(imt,jmt,max_blocks_clinic))
  allocate (roff (imt,jmt,max_blocks_clinic))
  allocate (ifrac(imt,jmt,max_blocks_clinic)) 
  allocate (patm (imt,jmt,max_blocks_clinic)) 
  allocate (duu10n(imt,jmt,max_blocks_clinic))

#endif

     ub = 0.0D0
     vb = 0.0D0
     h0 = 0.0D0
     ubp = 0.0D0
     vbp = 0.0D0
     h0p = 0.0D0
 
     U    = 0.0D0
     V    = 0.0D0
     UP   = 0.0D0
     VP   = 0.0D0
     WS   = 0.0D0
     H0L  = 0.0D0
     H0F  = 0.0D0
     H0BL = 0.0D0
     H0BF = 0.0D0
     UTL  = 0.0D0
     UTF  = 0.0D0
     VTL  = 0.0D0
     VTF  = 0.0D0
     NET  = 0.0D0
     PXB  = 0.0D0
     PYB  = 0.0D0
     PAX  = 0.0D0
     PAY  = 0.0D0
     WHX  = 0.0D0
     WHY  = 0.0D0
     WGP  = 0.0D0
 
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
      CALL YY00
!
 
      MONTH = 1
 
      IF (NSTART == 1) THEN
      number_day = 1
 
!     ------------------------------------------------------------------
!     READ LEVITUS ANNUAL MEAN TEMPERATURE AND SALINITY
!     ------------------------------------------------------------------
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      if (mytid==0) then
       iret=nf_open('TSinitial',nf_nowrite,ncid)
       call check_err (iret)
      end if

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
     
      do k=1,km !km
!
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_double(ncid,   5,start,count, buffer)
        call check_err (iret)
       end if
!
      call scatter_global(at(:,:,k,1,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
       if (mytid == 0 ) then
        iret=nf_get_vara_double(ncid,   6,start,count, buffer)
        call check_err (iret)
       end if
!
       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
       end do !km
 
       call POP_HaloUpdate(at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0_r8)
       call POP_HaloUpdate(at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0_r8)
!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK= 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                  AT (I,J,K,1,IBLOCK) = AT (I,J,K,1,IBLOCK)*VIT(I,J,K,IBLOCK)
                  AT (I,J,K,2,IBLOCK) = (AT (I,J,K,2,IBLOCK)- 35.0D0)*0.001D0*VIT(I,J,K,IBLOCK)
               END DO
            END DO
         END DO
      END DO
!
!$OMP PARALLEL DO PRIVATE (K)
            DO K = 1,KM
                     ATB (:,:,K,N,:) = AT (:,:,K,N,:)
#if (defined BOUNDARY)
                     RESTORE (:,:,K,N,:) = AT (:,:,K,N,:)
#endif
            END DO
 
      ATB (:,:,0,N,:) = 0.0D0
      ELSE
 
!     ------------------------------------------------------------------
!     READ INTERMEDIATE RESULTS (fort.22/fort.21)
!     ------------------------------------------------------------------
 
#if (defined BOUNDARY)
         if (mytid==0) then
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('TSinitial',nf_nowrite,ncid)
      call check_err (iret)

          end if
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
     
      do k=1,km
!
      if (mytid == 0) then
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=k ; count(3)=1
      start(4)=1 ; count(4)=1

      iret=nf_get_vara_double(ncid,   5,start,count, buffer)
      call check_err (iret)
      end if
!
      call scatter_global(at(:,:,k,1.:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
      if (mytid == 0 ) then
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=k ; count(3)=1
      start(4)=1 ; count(4)=1
      iret=nf_get_vara_double(ncid,   6,start,count, buffer)
      call check_err (iret)
      end if
!
      call scatter_global(at(:,:,k,2.:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
       end do
!
      if (mytid == 0 ) then
      iret = nf_close (ncid)
      call check_err (iret)
      end if
!!!!!!!!!!!!!!!!!!!!!!!
       call POP_HaloUpdate(at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0_r8)
       call POP_HaloUpdate(at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0_r8)
 
!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
     DO IBLOCK= 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                  RESTORE (I,J,K,1,IBLOCK) = at (I,J,K,1,IBLOCK)*VIT(I,J,K,IBLOCK) 
                  RESTORE (I,J,K,2,IBLOCK) = (at (I,J,K,2,IBLOCK) - 35.0)*0.001*VIT(I,J,K,IBLOCK) 
               END DO
            END DO
         END DO
      END DO
!
#endif
!
         if (mytid==0) then
         open(22,file=trim(out_dir)//fname,form='unformatted')
         end if

         if (mytid==0) then
         READ (22)buffer
         end if
         call scatter_global(h0,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         call POP_HaloUpdate(h0 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
         call scatter_global(u(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_swcorner, field_type_vector)
         end do
         call POP_HaloUpdate(u , POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
         call scatter_global(v(:,:,k,:),buffer, master_task, distrb_clinic, &
                          field_loc_swcorner, field_type_vector)
         end do
         call POP_HaloUpdate(v , POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
         call scatter_global(at(:,:,k,1,:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         end do
         call POP_HaloUpdate(at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter , &
                       POP_fieldKindScalar, errorCode, fillValue = 0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
         call scatter_global(at(:,:,k,2,:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         end do
         call POP_HaloUpdate(at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter , &
                       POP_fieldKindScalar, errorCode, fillValue = 0_r8)
!lhl20110728 for ws
         do k=1,km
         if (mytid==0) READ (22)buffer
         end do
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
         if (mytid==0) READ (22)buffer
!lhl20110728
         if (mytid == 0) then
          read(22)number_month,number_day
           month= number_month
         endif
            write(*,*) 'number_month =',number_month,'mon0=',mon0,&
                       'number_day=',number_day,'iday=',iday
#ifdef COUP
         if (nstart==2) then
            month=(cdate/10000-1)*12+mod(cdate,10000)/100
!M
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)            
            do iblock = 1, nblocks_clinic
            do j=1,jmt
            do i=1,imt
               ! lihuimin, TODO, to be considered
               ! lihuimin, 2012.7.23, coordinate with flux_cpl, ft. yu
               t_cpl (i,j,iblock)  = 273.15+at(i,j,1,1,iblock)
               s_cpl (i,j,iblock)  = at(i,j,1,2,iblock)*1000.+35.
               ! modi end
               q     (i,j,iblock)  = 0.0
               u_cpl (i,j,iblock)  = 0.0
               v_cpl (i,j,iblock)  = 0.0
               dhdx  (i,j,iblock)  = 0.0
               dhdy  (i,j,iblock)  = 0.0
            end do
            end do
            end do
         else
!LPF 20120815
!for t_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
          call scatter_global(t_cpl,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for s_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
          call scatter_global(s_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for u_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
!for v_cpl
          call scatter_global(u_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_vector)
          if (mytid==0) then
           READ (22)buffer
          end if
          call scatter_global(v_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_vector)

!for dhdx 
          if (mytid==0) then
           READ (22)buffer
          end if
          call scatter_global(dhdx,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for dhdy 
          if (mytid==0) then
           READ (22)buffer
          end if
          call scatter_global(dhdy, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for q
          if (mytid==0) then
           READ (22)buffer
          end if
          call scatter_global(q ,buffer,  master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         end if
!        end if
!LPF 20120815
!Yu
      call mpi_bcast(month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_day,1,mpi_integer,0,mpi_comm_ocn,ierr)

#endif
         CLOSE(22)
 
 
         NMFF = MOD (MONTH -1,12)
 
         CALL VINTEG (U,UB)
         CALL VINTEG (V,VB)
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1,jmt
               DO I = 1,IMT
                  UP (I,J,K,IBLOCK) = U (I,J,K,IBLOCK)
                  VP (I,J,K,IBLOCK) = V (I,J,K,IBLOCK)
                  UTF (I,J,K,IBLOCK) = U (I,J,K,IBLOCK)
                  VTF (I,J,K,IBLOCK) = V (I,J,K,IBLOCK)
                  ATB (I,J,K,1,IBLOCK) = AT (I,J,K,1,IBLOCK)
                  ATB (I,J,K,2,IBLOCK) = AT (I,J,K,2,IBLOCK)
               END DO
            END DO
         END DO
       END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1,jmt
            DO I = 1,IMT
               H0P (I,J,IBLOCK)= H0 (I,J,IBLOCK)
               UBP (I,J,IBLOCK)= UB (I,J,IBLOCK)
               VBP (I,J,IBLOCK)= VB (I,J,IBLOCK)
               H0F (I,J,IBLOCK)= H0 (I,J,IBLOCK)
               H0BF (I,J,IBLOCK)= H0 (I,J,IBLOCK)
            END DO
         END DO
      END DO
      END IF
!
#ifdef BIHAR
      call init_del4t
      call init_del4u
#else
      call init_del2t
      call init_del2u
#endif
!
      if (mytid==0)then
          write(6,*)"END-----------INIRUN !"
#ifdef COUP
          call shr_sys_flush(6)
#endif
      endif 

      deallocate(buffer)

      RETURN

      END SUBROUTINE INIRUN
 

