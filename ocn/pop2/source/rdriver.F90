!  CVS: $Id: rdriver.F90,v 1.7 2003/08/25 07:47:52 lhl Exp $
!     ======================
      SUBROUTINE RDRIVER
!     ======================
 
#include <def-undef.h>
use param_mod
use pconst_mod
use forc_mod
use work_mod
use dyn_mod, only : buffer
use msg_mod
use domain
use constant_mod
use shr_sys_mod

      IMPLICIT NONE
#include <netcdf.inc>
!
!     Define Variables.
      integer*4   :: ncid1, iret, iblock
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
 
!      REAL    :: WCOE (JMT),ABC
 
      allocate(su3(imt,jmt,12,max_blocks_clinic) , sv3(imt,jmt,12,max_blocks_clinic),psa3(imt,jmt,12,max_blocks_clinic), &
               tsa3(imt,jmt,12,max_blocks_clinic),qar3(imt,jmt,12,max_blocks_clinic),uva3(imt,jmt,12,max_blocks_clinic), &
               swv3(imt,jmt,12,max_blocks_clinic),cld3(imt,jmt,12,max_blocks_clinic),sss3(imt,jmt,12,max_blocks_clinic), &
               sst3(imt,jmt,12,max_blocks_clinic),nswv3(imt,jmt,12,max_blocks_clinic),dqdt3(imt,jmt,12,max_blocks_clinic), &
               chloro3(imt,jmt,12,max_blocks_clinic), seaice3(imt,jmt,12,max_blocks_clinic),  &
               runoff3(imt,jmt,12,max_blocks_clinic), wspd3(imt,jmt,12,max_blocks_clinic)  ,  &
               wspdu3(imt,jmt,12,max_blocks_clinic),wspdv3(imt,jmt,12,max_blocks_clinic),lwv3(imt,jmt,12,max_blocks_clinic), &
               rain3(imt,jmt,12,max_blocks_clinic),snow3(imt,jmt,12,max_blocks_clinic))
!
      allocate(buffer(imt_global,jmt_global))
!
      if (mytid==0)then
      write(6,*)"Beginning------RDRIVER! "
#ifdef COUP
      call shr_sys_flush(6)
#endif
      endif 

!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= 0.0
            END DO
         END DO
      END DO
    END DO
 
 
!-----------------------------------------------------------------------
!     SU3  : Sea surface zonal wind stress         (N/M**2)
!     SV3  : Sea surface meridional wind stress    (N/M**2)
!     PSA3 : Sea surface air pressure              (Pa)
!     SWV3 : Total net downward solar radiation    (W/M**2)
!    NSWV3 : None Solar flux                       (Wm-2)
!    DQDT3 : Dq/Dt                                 (WK-1m-2)
!     SST3 : Sea surface temperature               (Celsius)
!     SSS3 : Sea surface salinity                  (psu)
!   chloro3:chlorophll concentration               (mg m-3)
!-----------------------------------------------------------------------
 
!     READ FORCING FIELD
!
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      if(mytid==0)then
      iret=nf_open('MODEL.FRC',nf_nowrite,ncid1)
      call check_err (iret)
      end if
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      do k=1, 12
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=k ; count(4)=1

      if (mytid ==0 ) then
         iret=nf_get_vara_double(ncid1,   5,start,count,buffer)
         call check_err (iret)
      end if
      call scatter_global(swv3(:,:,k,:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
      if (mytid == 0) then
         iret=nf_get_vara_double(ncid1,   6,start,count,buffer)
         call check_err (iret)
      end if
      call scatter_global(nswv3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
      if (mytid == 0) then
          iret=nf_get_vara_double(ncid1,   7,start,count,buffer)
          call check_err (iret)
      end if
      call scatter_global(dqdt3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
      if (mytid == 0 ) then
          iret=nf_get_vara_double(ncid1,   8,start,count,buffer)
          call check_err (iret)
      end if
      call scatter_global(su3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
      if (mytid == 0 ) then
          iret=nf_get_vara_double(ncid1,   9,start,count,buffer)
          call check_err (iret)
      end if
      call scatter_global(sv3(:,:,k,:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
      if (mytid == 0) then
          iret=nf_get_vara_double(ncid1,  10,start,count,buffer)
          call check_err (iret)
      end if
      call scatter_global(sst3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
      if (mytid == 0 ) then
         iret=nf_get_vara_double(ncid1,  11,start,count,buffer)
         call check_err (iret)
      end if
      call scatter_global(sss3(:,:,k,:),buffer,  master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
     end do 
!
      if (mytid == 0 ) then
         iret = nf_close (ncid1)
         call check_err (iret)
      end if

!===============================================
!input seaice
#ifdef FRC_CORE
      if (mytid == 0 ) then
      iret=nf_open('seaice.db.NSIDC.1979-2006.clim.monthmean.modelgrid.01x01.nc',nf_nowrite,ncid1)
      call check_err (iret)
      end if

      do k=1, 12
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=k ; count(4)=1
      if (mytid == 0 ) then
      iret=nf_get_vara_double(ncid1, 5,start,count,buffer)
      call check_err (iret)
      end if
      call scatter_global(seaice3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      end do 

      if (mytid == 0 ) then
      iret = nf_close (ncid1)
      call check_err (iret)
      end if

!===============================================
!input runoff
      if (mytid == 0 ) then
      iret=nf_open('runoff.db.clim.modelgrid.01x01.nc',nf_nowrite,ncid1)
      call check_err (iret)
      end if

      do k=1, 1
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=k ; count(4)=1
      if (mytid == 0 ) then
      iret=nf_get_vara_double(ncid1, 4,start,count,buffer)
      call check_err (iret)
      end if
      call scatter_global(runoff3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      end do 

      if (mytid == 0 ) then
      iret = nf_close (ncid1)
      call check_err (iret)
      end if
#endif

!==============================================
!
      if (mytid==0) then
#ifdef COUP
      call shr_sys_flush(6)
#endif
      end if
 
!-----------------------------------------------------------------------
!     land grids of the forcing fields assigned to 0 
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
                SWV3(I,J,K,IBLOCK)    = SWV3(I,J,K,IBLOCK)*VIT(I,J,1,IBLOCK)
                NSWV3(I,J,K,IBLOCK)   = NSWV3(I,J,K,IBLOCK)*VIT(I,J,1,IBLOCK)
                DQDT3(I,J,K,IBLOCK)   = DQDT3(I,J,K,IBLOCK)*VIT(I,J,1,IBLOCK)
                SU3(I,J,K,IBLOCK)     =  SU3(I,J,K,IBLOCK)*VIV(I,J,1,IBLOCK)
                SV3(I,J,K,IBLOCK)     =  SV3(I,J,K,IBLOCK)*VIV(I,J,1,IBLOCK)
                SST3(I,J,K,IBLOCK)    = SST3(I,J,K,IBLOCK)*VIT(I,J,1,IBLOCK)
                SSS3(I,J,K,IBLOCK)    = SSS3(I,J,K,IBLOCK)*VIT(I,J,1,IBLOCK)
                seaice3(I,J,K,IBLOCK) = seaice3(I,J,K,IBLOCK)*VIT(I,J,1,IBLOCK)
                runoff3(I,J,K,IBLOCK) = runoff3(I,J,K,IBLOCK)*VIT(I,J,1,IBLOCK)
#if (defined SOLARCHLORO)
        chloro3(I,J,K,IBLOCK)= chloro3(I,J,K,IBLOCK)*VIT(I,J,1,IBLOCK) 
#endif
                END DO
           END DO
        END DO
      END DO

!-----------------------------------------------------------------------
!     salinity = (psu-35)*0.001
!     reverse VS (southward is positive)
!     notice: the former program VS is reversed during preparing forcing field
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SSS3 (I,J,K,IBLOCK) = (SSS3 (I,J,K,IBLOCK) -35.0D0)*0.001D0
               SV3  (I,J,K,IBLOCK) = -SV3 (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
    END DO
 
!-----------------------------------------------------------------------
!     CALCULATING THE ANNUAL MEAN FORCING FIELD 
!-----------------------------------------------------------------------
 
#if (defined FRC_ANN)
 
!$OMP PARALLEL DO PRIVATE (IBLOCK )
      DO IBLOCK = 1, NBLOCKS_CLINIC
      DO M = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,1,IBLOCK)= WKA (I,J,1,IBLOCK) + SU3 (I,J,M,IBLOCK)
               WKA (I,J,2,IBLOCK)= WKA (I,J,2,IBLOCK) + SV3 (I,J,M,IBLOCK)
               WKA (I,J,3,IBLOCK)= WKA (I,J,3,IBLOCK) + SSS3 (I,J,M,IBLOCK)
               WKA (I,J,4,IBLOCK)= WKA (I,J,4,IBLOCK) + SWV3 (I,J,M,IBLOCK)
               WKA (I,J,5,IBLOCK)= WKA (I,J,5,IBLOCK) + SST3 (I,J,M,IBLOCK)
               WKA (I,J,6,IBLOCK)= WKA (I,J,6,IBLOCK) + NSWV3 (I,J,M,IBLOCK)
               WKA (I,J,7,IBLOCK)= WKA (I,J,7,IBLOCK) + DQDT3 (I,J,M,IBLOCK)
#if (defined SOLARCHLORO)
               WKA (I,J,8,IBLOCK)= WKA (I,J,8,IBLOCK) + chloro3 (I,J,M,IBLOCK)
#endif
            END DO
         END DO
      END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (M,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
      DO M = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SU3 (I,J,M,IBLOCK) = WKA (I,J,1,IBLOCK)/12.0D0
               SV3 (I,J,M,IBLOCK) = WKA (I,J,2,IBLOCK)/12.0D0
              SSS3 (I,J,M,IBLOCK) = WKA (I,J,3,IBLOCK)/12.0D0
              SWV3 (I,J,M,IBLOCK) = WKA (I,J,4,IBLOCK)/12.0D0
              SST3 (I,J,M,IBLOCK) = WKA (I,J,5,IBLOCK)/12.0D0
             NSWV3 (I,J,M,IBLOCK) = WKA (I,J,6,IBLOCK)/12.0D0
             DQDT3 (I,J,M,IBLOCK) = WKA (I,J,7,IBLOCK)/12.0D0
#if (defined SOLARCHLORO)
             chloro3 (I,J,M,IBLOCK) = WKA (I,J,8,IBLOCK)/12.0D0
#endif
            END DO
         END DO
      END DO
      END DO
 
#endif
!
      if (mytid==0)then
      write(6,*)"END-----------RDRIVER!"
#ifdef COUP
      call shr_sys_flush(6)
#endif
      endif 
 
    deallocate(buffer)

      RETURN
      END SUBROUTINE RDRIVER
 
