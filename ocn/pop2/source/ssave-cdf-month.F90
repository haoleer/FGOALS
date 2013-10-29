!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEMON
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
!     =================
!     output history (netcdf) and restart (binary) files
!     remove pre-compilation NETCDF/ALL/NORMAL, remove yearly mean part 
!     keep the diagnostic part, SPMD
!     written by liu hailong & Lin Pengfei 2012 July
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use output_mod
use dyn_mod
use tracer_mod
use cdf_mod
use diag_mod
use msg_mod
use domain
use distribution
use gather_scatter

      IMPLICIT none
#include <netcdf.inc>
 
      logical :: hist_output,rest_output
      CHARACTER ( LEN =   4 ) :: ftail
      CHARACTER ( LEN =  24 ) :: fname
      CHARACTER ( LEN =  15 ) :: fname1
      CHARACTER ( LEN =   8 ) :: dd
      CHARACTER ( LEN =   10 ) :: tt
      CHARACTER ( LEN =   5 ) :: zz
      INTEGER(r4)             :: vv(8)
      INTEGER :: nwmf, iblock
 
!---------------------------------------------------------------------
!     output monthly results
!---------------------------------------------------------------------
!    file name
 
      nwmf = iyfm
      write (ftail,'(i4.4)') nwmf
      fname1(1:5)='MMEAN'
      fname1(6:9)=ftail
      fname1(10:10)='-'
      write(fname1(11:12),'(i2.2)')mon0
      fname1(13:15)='.nc'
 
!      if (iday==imd) then
      if (mod (iyfm,io_hist)==0 ) then
         hist_output=.true.
      else
         hist_output=.false.
      endif
!
      if (mod ((month-1),io_rest)==0 ) then
         rest_output=.true.
      else
         rest_output=.false.
      endif

       !write(*,*) 'hist_output,io_hist',hist_output,io_hist
       if(mytid==0) write(*,*) 'ok inssavemon'


!--------------------------------------------------------------
!     cdf output
!--------------------------------------------------------------
!     file defination
! enter define mode
       if (mytid==0) then
         iret = nf_create (fname1, NF_CLOBBER, ncid)
         CALL check_err (iret)
          write(*,*)'ok generate,ncid=',ncid, time_len
! define dimensions
         iret = nf_def_dim (ncid, 'lat', lat_len, lat_dim)
         CALL check_err (iret)
         iret = nf_def_dim (ncid, 'lon', lon_len, lon_dim)
         CALL check_err (iret)
 
         IF (mon0 == 12) THEN
            iret = nf_def_dim (ncid, 'lev', lev_len, lev_dim)
            CALL check_err (iret)
         ELSE
            iret = nf_def_dim (ncid, 'lev', klv, lev_dim)
            CALL check_err (iret)
         END IF
 
         iret = nf_def_dim (ncid, 'lev1', lev1_len, lev1_dim)
         CALL check_err (iret)
 
         iret = nf_def_dim (ncid, 'time', NF_UNLIMITED, time_dim)
!      iret = nf_def_dim(ncid, 'time', time_len, time_dim)
         CALL check_err (iret)
! define variables
         lat_dims (1) = lat_dim
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, lat_dims, lat_id)
         CALL check_err (iret)
!
         lon_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, lon_dims, lon_id)
         CALL check_err (iret)
!
         lev_dims (1) = lev_dim
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, lev_dims, lev_id)
         CALL check_err (iret)
!
         lev1_dims (1) = lev1_dim
         iret = nf_def_var (ncid, 'lev1', NF_REAL, lev1_rank, lev1_dims, lev1_id)
         CALL check_err (iret)
!
         time_dims (1) = time_dim
         iret = nf_def_var (ncid, 'time', NF_DOUBLE, time_rank, time_dims, time_id)
         CALL check_err (iret)
 
         z0_dims (3) = time_dim
         z0_dims (2) = lat_dim
         z0_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'z0', NF_REAL, z0_rank, z0_dims, z0_id)
         CALL check_err (iret)
 
         ic1_dims (3) = time_dim
         ic1_dims (2) = lat_dim
         ic1_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ic1', NF_REAL, ic1_rank, ic1_dims, ic1_id)
         CALL check_err (iret)
 
         ic2_dims (3) = time_dim
         ic2_dims (2) = lat_dim
         ic2_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ic2', NF_REAL, ic2_rank, ic2_dims, ic2_id)
         CALL check_err (iret)

         net1_dims (3) = time_dim
         net1_dims (2) = lat_dim
         net1_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'net1', NF_REAL, net1_rank, net1_dims, net1_id)
         CALL check_err (iret)
 
         net2_dims (3) = time_dim
         net2_dims (2) = lat_dim
         net2_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'net2', NF_REAL, net2_rank, net2_dims, net2_id)
         CALL check_err (iret)
!
         mld_dims (3) = time_dim
         mld_dims (2) = lat_dim
         mld_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'mld', NF_REAL, mld_rank, mld_dims, mld_id)
         CALL check_err (iret)
         akm_dims (4) = time_dim
         akm_dims (3) = lev_dim
         akm_dims (2) = lat_dim
         akm_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'akm', NF_REAL, akm_rank, akm_dims, akm_id)
         CALL check_err (iret)
         akt_dims (4) = time_dim
         akt_dims (3) = lev_dim
         akt_dims (2) = lat_dim
         akt_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'akt', NF_REAL, akt_rank, akt_dims, akt_id)
         CALL check_err (iret)
         aks_dims (4) = time_dim
         aks_dims (3) = lev_dim
         aks_dims (2) = lat_dim
         aks_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'aks', NF_REAL, aks_rank, aks_dims, aks_id)
         CALL check_err (iret)
         ts_dims (4) = time_dim
         ts_dims (3) = lev_dim
         ts_dims (2) = lat_dim
         ts_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ts', NF_REAL, ts_rank, ts_dims, ts_id)
         CALL check_err (iret)
         ss_dims (4) = time_dim
         ss_dims (3) = lev_dim
         ss_dims (2) = lat_dim
         ss_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ss', NF_REAL, ss_rank, ss_dims, ss_id)
         CALL check_err (iret)
         us_dims (4) = time_dim
         us_dims (3) = lev_dim
         us_dims (2) = lat_dim
         us_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'us', NF_REAL, us_rank, us_dims, us_id)
         CALL check_err (iret)
         vs_dims (4) = time_dim
         vs_dims (3) = lev_dim
         vs_dims (2) = lat_dim
         vs_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'vs', NF_REAL, vs_rank, vs_dims, vs_id)
         CALL check_err (iret)
         ws_dims (4) = time_dim
         ws_dims (3) = lev_dim
         ws_dims (2) = lat_dim
         ws_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'ws', NF_REAL, ws_rank, ws_dims, ws_id)
         CALL check_err (iret)
         su_dims (3) = time_dim
         su_dims (2) = lat_dim
         su_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'su', NF_REAL, su_rank, su_dims, su_id)
         CALL check_err (iret) !Uwindstress

         sv_dims (3) = time_dim
         sv_dims (2) = lat_dim
         sv_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'sv', NF_REAL, sv_rank, sv_dims, sv_id)
         CALL check_err (iret) !Vwindstress
         lthf_dims (3) = time_dim
         lthf_dims (2) = lat_dim
         lthf_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'lthf', NF_REAL, lthf_rank, lthf_dims, lthf_id)
         CALL check_err (iret) !latent heat flux 
         sshf_dims (3) = time_dim
         sshf_dims (2) = lat_dim
         sshf_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'sshf', NF_REAL, sshf_rank, sshf_dims, sshf_id)
         CALL check_err (iret) !sensible heat flux 
         lwv_dims (3) = time_dim
         lwv_dims (2) = lat_dim
         lwv_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'lwv', NF_REAL, lwv_rank, lwv_dims, lwv_id)
         CALL check_err (iret) !longwave 
         swv_dims (3) = time_dim
         swv_dims (2) = lat_dim
         swv_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'swv', NF_REAL, swv_rank, swv_dims, swv_id)
         CALL check_err (iret) !shortwave 

         if (diag_msf) then
            psi_dims (3) = time_dim
            psi_dims (2) = lev1_dim
            psi_dims (1) = lat_dim
            iret = nf_def_var (ncid, 'psi', NF_REAL, psi_rank, psi_dims, psi_id)
            CALL check_err (iret)
         end if

         if (diag_bsf) then
            bsf_dims (3) = time_dim
            bsf_dims (2) = lat_dim
            bsf_dims (1) = lon_dim
            iret = nf_def_var (ncid, 'bsf', NF_REAL, bsf_rank, bsf_dims, bsf_id)
            CALL check_err (iret)
         end if
 
         if (diag_mth) then
            mth_dims (3) = time_dim
            mth_dims (2) = lat_dim
            mth_dims (1) = lon_dim
            iret = nf_def_var (ncid, 'mth', NF_REAL, mth_rank, mth_dims, mth_id)
            CALL check_err (iret)
         end if

#if (defined SMAG_OUT)
         am_dims (4) = time_dim
         am_dims (3) = lev_dim
         am_dims (2) = lat_dim
         am_dims (1) = lon_dim
         iret = nf_def_var (ncid, 'am', NF_REAL, am_rank, am_dims, am_id)
         CALL check_err (iret)
#endif
 
! assign attributes
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'long_name', 18, 'depth (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev1_id, 'long_name', 18, 'depth (on V grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev1_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, time_id, 'long_name', 4, 'time')
         CALL check_err (iret)
!      iret = nf_put_att_text(ncid, time_id, 'units', 21, 'days since 1001-01-01')
         iret = nf_put_att_text (ncid, time_id, 'units', 23, 'months since 0001-01-01')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, z0_id, 'long_name', 18, 'sea surface height')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, z0_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, z0_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ic1_id, 'long_name', 56, 'total of number of levels involved in convection per day')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ic1_id, 'units', 6, 'levels')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ic1_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ic2_id, 'long_name', 35, 'number of levels ventilated per day')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ic2_id, 'units', 6, 'levels')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ic2_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, net1_id, 'long_name', 21, 'net surface heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net1_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, net1_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, net2_id, 'long_name', 21, 'net surface salt flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net2_id, 'units', 7, 'psu/sec')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, net2_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, mld_id, 'long_name', 17, 'mixed layer depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, mld_id, 'units', 1, 'm')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, mld_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, akm_id, 'long_name', 28, 'turbulent vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, akm_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, akm_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, akt_id, 'long_name', 33, 'turbulent heat vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, akt_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, akt_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, aks_id, 'long_name', 33, 'turbulent salt vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, aks_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, aks_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, ts_id, 'long_name', 11, 'temperature')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ts_id, 'units', 10, 'centigrade')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ts_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ss_id, 'long_name', 8, 'salinity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ss_id, 'units', 3, 'psu')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ss_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, us_id, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, us_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, us_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, vs_id, 'long_name', 18, 'meridional current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, vs_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, vs_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ws_id, 'long_name', 16, 'vertical current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ws_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ws_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, su_id, 'long_name', 11, 'Uwindstress')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, su_id, 'units', 2, 'Pa')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, su_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, sv_id, 'long_name', 11, 'Vwindstress')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, sv_id, 'units', 2, 'Pa')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, sv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, lthf_id, 'long_name', 16, 'latent heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lthf_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, lthf_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)


         iret = nf_put_att_text (ncid, sshf_id, 'long_name', 18, 'sensible heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, sshf_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, sshf_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, lwv_id, 'long_name', 8, 'Longwave')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lwv_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, lwv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, swv_id, 'long_name', 9, 'Shortwave')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, swv_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, swv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
!
         if (diag_msf) then
             iret = nf_put_att_text (ncid, psi_id, 'long_name', 27, 'Meridioanl Stream Function')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, psi_id, 'units', 8, 'Sverdrup')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, psi_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
         end if

         if (diag_bsf) then
             iret = nf_put_att_text (ncid, bsf_id, 'long_name', 26, 'Barotropic Stream Function')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, bsf_id, 'units', 8, 'Sverdrup')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, bsf_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
         end if
!
         if (diag_mth) then
             iret = nf_put_att_text (ncid, mth_id, 'long_name', 27, 'Meridional Tracer Transport')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, mth_id, 'units', 9, 'PW or m/s')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, mth_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
         end if
!
#if (defined SMAG_OUT)
         iret = nf_put_att_text (ncid, am_id, 'long_name', 20, 'horizontal viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, am_id, 'units', 6, 'm**2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, am_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
#endif
!   define global attribute
         CALL date_and_time (dd,tt,zz,vv)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', 4, 'test')
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', 20, tt //'  '//dd)
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
         CALL check_err (iret)
! leave define mode
         iret = nf_enddef (ncid)
         CALL check_err (iret)

!----------------------------------------------------------
!     prepare data for storing
!----------------------------------------------------------
 
         t0_cdf = month -1 !need to test 
         iret = nf_put_var_real (ncid, lon_id, lon)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lat_id, lat)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev1_id, lev1)
         CALL check_err (iret)
         start1 (1)= 1
         count1 (1)= 1
         iret = nf_put_vara_double (ncid, time_id,start1,count1,t0_cdf)
         CALL check_err (iret)
 
         endif !mytid==0

! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= 1
 
         allocate(buffer_r4_global(imt_global,jmt_global), buffer_r4_local(imt,jmt,max_blocks_clinic))
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = z0mon/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
         if(mytid==0) then         
            iret = nf_put_vara_real (ncid,z0_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = icmon(:,:,1,:)/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,ic1_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0

!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = icmon(:,:,2,:)/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,ic2_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = netmon(:,:,1,:)/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,net1_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0

         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = netmon(:,:,2,:)/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,net2_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = mldmon/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,mld_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!3D output
         
         start4 (1)= 1
         start4 (2)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (4)= 1

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( viv(:,:,k,:) > 0.5D0 )
               buffer_r4_local = akmmon(:,:,k,:)/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,akm_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do
!
         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = aktmon(:,:,k,:)/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,akt_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = aksmon(:,:,k,:)/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,aks_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = tsmon(:,:,k,:)/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,ts_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = ssmon(:,:,k,:)/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,ss_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = wsmon(:,:,k,:)/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,ws_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( viv(:,:,k,:) > 0.5D0 )
               buffer_r4_local = usmon(:,:,k,:)/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,us_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( viv(:,:,k,:) > 0.5D0 )
               buffer_r4_local = vsmon(:,:,k,:)/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,vs_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

!Taux 

         where ( viv(:,:,1,:) > 0.5D0 )
            buffer_r4_local = sumon/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,su_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( viv(:,:,1,:) > 0.5D0 )
            buffer_r4_local = svmon/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,sv_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = lthfmon/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,lthf_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = sshfmon/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,sshf_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = lwvmon/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,lwv_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = swvmon/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,swv_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
     if (diag_bsf) then
        call barosf
        if (mytid ==0 ) then
            iret = nf_put_vara_real (ncid,bsf_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
        end if
     end if
!
     if (diag_msf) then
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lat_len
         count3 (2)= lev1_len
         count3 (3)= 1
        call msf
        if (mytid ==0 ) then
            iret = nf_put_vara_real (ncid,psi_id,start3, count3, psi(:,:,1))
            CALL check_err (iret)
        end if
     end if
!
     if (diag_mth) then
        call diag_heat_transport(1)
        call diag_heat_transport(2)
        if (mytid ==0 ) then
            iret = nf_put_vara_real (ncid,mth_id,start3, count3, mth)
            CALL check_err (iret)
        end if
     end if
!
     if(mytid==0) then 
         iret = nf_CLOSE (ncid)
         CALL check_err (iret)
     endif

 
         CALL mm00 (klv)
!lhl20120728      IF (mod ( (month -1),12) == 0)THEN
      IF (mod ( (month -1),io_rest) == 0)THEN
         CALL yy00
      END IF
!

! 
      deallocate (buffer_r4_local,buffer_r4_global)
      RETURN
      END 
