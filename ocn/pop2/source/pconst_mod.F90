!  CVS: $Id: pconst_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module pconst_mod
!
#include <def-undef.h>
use precision_mod
use param_mod
!     -----------------------------------------------------------
!     Index Fields
!     -----------------------------------------------------------
!YU   real,dimension(:),allocatable:: dyr_global
      integer,dimension(jmt):: j_global
      integer,dimension(imt):: i_global
      integer :: ix,iy
      real(r8),dimension(imt,jmt,km,max_blocks_clinic):: vit,viv 
!     real(r8),dimension(imt_global,jmt,km):: vit_1d,viv_1d
      real(r8),dimension(jmt_global):: ahv_back
      integer,dimension(imt,jmt,max_blocks_clinic):: basin
!     integer,dimension(imt,jmt,max_blocks_clinic):: itnu
!lhl1204
      integer,dimension(imt,jmt,max_blocks_clinic):: na
      real(r8) :: dfricmx,dwndmix

      ! lihuimin, 2012.7.15
!     integer :: i_num   ! actual grid number in i direction in this process
!     integer :: j_num   ! actual grid number in j direction in this process
!     integer :: i_f_num ! formal grid number in i dircetion 
!     integer :: j_f_num ! formal grid number in j direction !lhl1204
!
!
!     -----------------------------------------------------------
!     Grids
!     -----------------------------------------------------------
#if (defined NETCDF) || (defined ALL)
      real(r4),dimension(imt_global):: lon
      real(r4),dimension(jmt_global):: lat
      real(r4),dimension(km):: lev
      real(r4),dimension(km+1):: lev1
#endif
!lhl090729
      real(r8),dimension(s_imt):: s_lon
      real(r8),dimension(s_jmt):: s_lat
!lhl090729
      real(r8),dimension(km):: zkt,dzp,odzp,odzt
      real(r8),dimension(kmp1):: zkp
      real(r8),dimension(imt,jmt,max_blocks_clinic):: EBEA,EBEB,EBLA,EBLB,EPEA,EPEB,EPLA,EPLB
      real(r8),dimension(imt,jmt,max_blocks_clinic):: RRD1,RRD2
!lhl060506

#if ( defined SMAG)

#endif

!Yu

      real(r8),dimension(imt,jmt,max_blocks_clinic)::ohbt,ohbu,dzph,hbx,hby
      real(r8),dimension(imt,jmt,max_blocks_clinic):: SNLAT
!     real(r8),dimension(jmt):: COSU,COST
!     real(r8),dimension(:),allocatable::COSU_global,COST_global
      integer,dimension(:),allocatable:: i_start,j_start
!     real(r8),dimension(imt)::CF1,CF2,SF1,SF2
!
!
!     -----------------------------------------------------------
!     Reference T S & coefficients for calculation of d(density)
!     -----------------------------------------------------------

      REAL(r8):: TO(KM),SO(KM),C(KM,9),PO(KM)
!lhl1204

!YU
!
!
!     -----------------------------------------------------------
!     Control Parameter
!     -----------------------------------------------------------
      INTEGER:: ISOP
!
!
!     -----------------------------------------------------------
!     Phycical Parameter
!     -----------------------------------------------------------
      real(r8)::amv,ahv,ahice
!lhl1204
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::akmu,akmt
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::akt
!lhl1204
      real(r8),dimension(jmt)::AM,AH
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::am3,ah3
!lhl
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::amx,amy
!lhl
      real(r8)::gamma
#if ( defined SMAG)
      real(r8):: D0,CP,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
                    VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2,KARMAN,RR
#else
      real(r8)::  D0,CP,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
                    VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2
#endif
      logical :: diag_msf, diag_bsf, diag_mth, diag_budget
!
!
      CHARACTER (LEN=3):: ABMON(12)
      CHARACTER (LEN=3):: ABMON1(12)
      CHARACTER (LEN=80):: out_dir
!
      REAL(r8):: DTB,DTC,DTS,DTB2,DTC2,ONBB,ONBC,ONCC
      INTEGER:: NBB,NCC,NSS,ISB,ISC,IST,MONTH
      INTEGER:: NMONTH(12),NNMONTH(12)
      INTEGER:: number_day, number_month
!
      ! lihuimin 2012.6.18, add REFDATE
      INTEGER :: NUMBER,NSTART,IY0,IYFM,MON0,MEND,IMD,IDAY,II,JJ,IO_HIST,IO_REST,rest_freq,hist_freq,REFDATE
      integer :: klv
!
      character (len=80) :: adv_momentum, adv_tracer
      integer(kind(1))   :: curr_ymd_licom     ! Current date YYYYMMDD
      integer(kind(1))   :: curr_ymd_cpl     ! Current date YYYYMMDD
      integer(kind(1))   :: yy_licom,mm_licom,dd_licom,tod_licom     ! year, month, day
      integer(kind(1))   :: yy_cpl,mm_cpl,dd_cpl,tod_cpl     ! year, month, day
      integer:: ocn_cpl_dt
end module pconst_mod
