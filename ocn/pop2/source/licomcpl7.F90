#define LOGMSG()
!write(mytid+600,'(a,i4)')"LICOM",__LINE__
!#define LOGMSGCarbon()
!write(600+mytid,'(a,i4)')"Licom CARBON",__LINE__


module licom_comp_mct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !MODULE:    LICOM_COMP_MCT
!   !AUTHOR:    Huimin Li
!   !Date:      2012/5/30
!
!   !DESCRIPTION:
!               This is the main driver for the licom coupled in CPL7
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use mct_mod
   use esmf_mod
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use shr_file_mod
   use shr_cal_mod, only : shr_cal_date2ymd
   use shr_sys_mod
   use shr_const_mod,only:SHR_CONST_SPVAL !linpf 2012Jul26
   use perf_mod
   use fluxcpl
   use POP_CplIndices
   use shr_dmodel_mod
   use POP_CommMod
   use domain
   use grid
   use blocks
   use gather_scatter
   use distribution

#include <def-undef.h>
use param_mod
use pconst_mod
use shr_msg_mod
use shr_sys_mod
use control_mod
use constant_mod
use shr_cal_mod,       only: shr_cal_date2ymd
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use tracer_mod
use pmix_mod
use forc_mod
!
#ifdef USE_OCN_CARBON
use carbon_mod
use cforce_mod
#endif


  implicit none
#include <netcdf.inc>
  public :: licom_init_mct
  public :: licom_run_mct
  public :: licom_final_mct
  SAVE
  private



  private :: licom_export_mct
  private :: licom_import_mct
  private :: licom_SetGSMap_mct
  private :: licom_domain_mct


  type(seq_infodata_type), pointer :: &
     infodata

!==========================================================================
  contains
!==========================================================================


!**************************************************************************
!   !ROUTINE:   licom_init_mct
!   !AUTHOR:    Huimin Li
!   !Date:      2012/5/30
!
!   !INTERFACE:
  subroutine licom_init_mct(EClock, cdata_o, x2o_o, o2x_o, NLFilename)
!
!   !DESCRIPTION:
!               This is the initialize routine for licom coupled in CPL7
!   !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!   !LOCAL VARIABLES
    type(mct_gsMap), pointer :: &
       gsMap_o

    type(mct_gGrid), pointer :: &
       dom_o

    integer (kind(1)) :: &
       nThreads

    integer (kind(1)) :: &
       OCNID,     &
       lsize

    real (r8) ::  &
       precadj

    integer (kind(1)) :: iam,ierr
    character(len=32)  :: starttype          ! infodata start type

    integer :: i_temp, j_temp, i_comp, j_comp ! used in specifying the grid number in each process

  !------------------------------------------------------
  ! ROUTINE BEGIN
  !------------------------------------------------------
  mpi_comm_ocn=0

  ! lihuimin 2012.7.16
  ! OCNID, errcode
  call seq_cdata_setptrs(cdata_o, ID=OCNID, mpicom=mpi_comm_ocn, &
       gsMap=gsMap_o, dom=dom_o, infodata=infodata)

  ! five parameters initialized in msg_pass('init') in CPL6 coupled version
  ! get infodata from drv
  cdate = 0
  sec = 0
  ierr = 0
  info_time = 0
  call seq_infodata_GetData( infodata, info_debug=info_dbug)

  ! send initial state to drv
  call seq_infodata_PutData( infodata, ocn_nx=imt_global, ocn_ny=jmt_global)
!  call seq_infodata_PutData( infodata, ocn_prognostic=.true.) !LPF 20120829
  ! lihuimin, 2012.7.25, not use ocnrof_p
  call seq_infodata_PutData( infodata, ocn_prognostic=.true.,ocnrof_prognostic=.true.,rof_present=.true.)
!open ocnrof. !LPF 20120829

  !-----------------------------
  ! namelist & logfile setup
  !-----------------------------
  ! TODO: redirect the log file


  mytid=0

! He - 2010-10-08 | NDay: Number of days since the first day when the simulation begins.
! lihuimin cancle , 2012.6.14
!     NDAY = 0



!
! mpicom_o <-> mpi_comm_ocn
! TODO
      if (mytid==0) write(6,*)"Begin mpi_comm_rank"
      call mpi_comm_rank (mpi_comm_ocn, mytid, ierr)
      if (mytid==0)  write(6,*)"End mpi_comm_rank"
      call mpi_comm_size (mpi_comm_ocn, nproc, ierr)
      write(6,*) "MYTID=",mytid,"Number of Processors is",nproc
!YU

      my_task = mytid
      master_task = 0
      POP_mytask = mytid
      POP_mastertask = 0
!YU


   ! lihuimin 2012.6.14, initial indices
   call POP_CplIndicesSet()


!---------------------------------------------------------------------
!     SET THE CONSTANTS USED IN THE MODEL
!---------------------------------------------------------------------
#ifdef COUP
      call shr_sys_flush(6)
#endif
    LOGMSG()
      CALL CONST
    LOGMSG()
      if (mytid == 0) then
      write(111,*)"OK------3"
      close(111)
      end if
#ifdef COUP
      call shr_sys_flush(6)
#endif
#ifdef SHOW_TIME
      call run_time('CONST')
#endif

!***********************************************************************
!          SET SOME CONSTANTS FOR THE BOGCM
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL CTRLC
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('CTRLC')
#endif
#endif

!---------------------------------------------------------------------
!     SET MODEL'S RESOLUTION,TOPOGRAPHY AND THE CONSTANT
!     PARAMETERS RELATED TO LATITUDES (J)
!---------------------------------------------------------------------
    LOGMSG()
   call init_domain_blocks
   call init_grid1
   call init_domain_distribution(KMT_G)
   call init_grid2
   CALL GRIDS
   call calc_coeff
      if (mytid == 0) then
      write(111,*)"OK------4"
      close(111)
      end if
#ifdef SHOW_TIME
      call run_time('GRIDS')
#endif


!---------------------------------------------------------------------
!     SET SURFACE FORCING FIELDS (1: Annual mean; 0: Seasonal cycle)
!---------------------------------------------------------------------
    LOGMSG()
#ifdef SHOW_TIME
      call run_time('RDRIVER')
#endif
      if (mytid == 0) then
      write(111,*)"OK------5"
      close(111)
      end if

!***********************************************************************
!      SET FORCING DATA USED IN CARBON CYCLYE
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
       CALL CFORCE
      LOGMSGCarbon()
#ifdef SHOW_TIME
       call run_time('CFORCE')
#endif
#endif

!---------------------------------------------------------------------
!     INITIALIZATION
!---------------------------------------------------------------------
    LOGMSG()
      CALL INIRUN
      if (mytid == 0) then
      write(111,*)"OK------5.0"
      close(111)
      end if
#ifdef SHOW_TIME
      call run_time('INRUN')
#endif

!----------------------------------------------------------------------
!     Inialize mct attribute vectors
!     imitate CPL7/POP
!     lihuimin 2012.7.16
!----------------------------------------------------------------------

   call licom_SetGSMap_mct(mpi_comm_ocn, OCNID, GSMap_o)
      if (mytid == 0) then
      write(111,*)"OK------6.0"
      close(111)
      end if

   lsize = mct_gsMap_lsize(gsMap_o, mpi_comm_ocn)
   write(6,*) "mct_gsMap_lsize = ",lsize

   call licom_domain_mct(lsize, gsMap_o, dom_o)

   call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
   call mct_aVect_zero(x2o_o)


   call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize)
   call mct_aVect_zero(o2x_o)
      if (mytid == 0) then
      write(111,*)"OK------7.0"
      close(111)
      end if



!**********************************************************************
!      INITIALIZATION CARBON MODEL
!**********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL INIRUN_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INIRUN_PT')
#endif
#endif
#ifdef CANUTO
      call turb_ini
#endif
!---------------------------------------------------------------------
!     INITIALIZATION OF ISOPYCNAL MIXING
!---------------------------------------------------------------------
    LOGMSG()
#ifdef ISO
      CALL ISOPYI
#ifdef SHOW_TIME
      call run_time('ISOPYI')
#endif
#endif

#ifdef COUP
         call shr_sys_flush(6)
#endif

    call licom_export_mct(o2x_o)

      if (mytid == 0) then
      write(111,*)"OK------8.0"
      close(111)
      end if

  end subroutine licom_init_mct


!**************************************************************************
!   !ROUTINE:   licom_run_mct
!   !AUTHOR:    Huimin Li
!   !Date:      2012/6/2
!
!   !INTERFACE:
  subroutine licom_run_mct( EClock, cdata_o, x2o_o, o2x_o)
!
!   !DESCRIPTION:
!   !run licom for a coupling interval
!
!   !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o

    integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD
    integer(kind(1))   :: yy,mm,dd     ! year, month, day

!----------------------------------------------------------------------
!   get Eclock data information
!   lihuimin, 2012.7.20
!----------------------------------------------------------------------
    call seq_timemgr_EClockGetData( EClock, curr_ymd=curr_ymd)
    call shr_cal_date2ymd(curr_ymd,yy,mm,dd)
!    imd = 30
    iyfm = yy  ! nwmf = iyfm
    mon0 = mm  ! month = (iyfm-1)*12 + mon0
    iday = dd  ! number_day = iday + 1
    
    if(mytid == 0) then
       write(6,*) "From CPL7-Time yy=",yy,"mm=",mm,"dd=",dd
    endif

!linpf 2012Jul27
!=====================time control for licom output===============
!lhl20120728
       month=(IYFM-1)*12+mon0
       IMD = NMONTH (MON0)
!lhl20120728
!      month=mm
!      IY0 = (MONTH -1)/12
!      IYFM = IY0+1
!!     IYFM is the number of the current year
!
!      MON0 = MONTH - IY0*12
!      IMD = NMONTH (MON0)
!!=================================================================
!linpf 2012Jul27

!----------------------------------------------------------------------
!   call msg_pass('recv')
!----------------------------------------------------------------------
    LOGMSG()
      if (mytid == 0) then
      write(111,*)"OK------9.0"
      close(111)
      end if
    call licom_import_mct(x2o_o)
      if (mytid == 0) then
      write(111,*)"OK------10.0"
      close(111)
      end if

    LOGMSG()
    call post_cpl

      if (mytid == 0) then
      write(111,*)"OK------11.0"
      close(111)
      end if
!---------------------------------------------------------------------
!     THERMAL CYCLE
!---------------------------------------------------------------------
    LOGMSG()
         DO II = 1,NSS


!     COMPUTE DENSITY, BAROCLINIC PRESSURE AND THE RELAVANT VARIABLES
    LOGMSG()
      if (mytid == 0) then
      write(111,*)"OK------12.0"
      close(111)
      end if
            CALL READYT
      if (mytid == 0) then
      write(111,*)"OK------13.0"
      close(111)
      end if

!----------------------------------------------------
!     BAROCLINIC & BAROTROPIC CYCLE
!---------------------------------------------------------------------
    LOGMSG()
            DO JJ = 1,NCC

!     COMPUTE MOMENTUM ADVECTION, DIFFUSION & THEIR VERTICAL INTEGRALS
    LOGMSG()
               CALL READYC
      if (mytid == 0) then
      write(111,*)"OK------14.0"
      close(111)
      end if

!     PREDICTION OF BAROTROPIC MODE
    LOGMSG()
      call energy
               CALL BAROTR
      call energy
      if (mytid == 0) then
      write(111,*)"OK------15.0"
      close(111)
      end if


!     PREDICTION OF BAROCLINIC MODE
    LOGMSG()
               CALL BCLINC
      call energy
      if (mytid == 0) then
      write(111,*)"OK------16.0"
      close(111)
      end if



            END DO

!*******************************************************************
!     PREDICTION OF PASSIVE TRACER
!*******************************************************************

    LOGMSG()
      if (mytid == 0) then
      write(111,*)"OK------17.0"
      close(111)
      end if
            CALL TRACER
      call energy
      stop
      if (mytid == 0) then
      write(111,*)"OK------18.0"
      close(111)
      end if

    LOGMSG()
            CALL ICESNOW

      if (mytid == 0) then
      write(111,*)"OK------19.0"
      close(111)
      end if

!***********************************************************************
!     PERFORM CONVECTIVE ADJUSTMENT IF UNSTABLE STRATIFICATION OCCURS
!************************************************************************
!lhl1204
!#if (!defined CANUTO)
    LOGMSG()
            CALL CONVADJ
!#endif
      if (mytid == 0) then
      write(111,*)"OK------20.0"
      close(111)
      end if

!*************************************************************************
!
         END DO
    LOGMSG()
      if (mytid == 0) then
      write(111,*)"OK------21.0"
      close(111)
      end if
         CALL ENERGY

      if (mytid == 0) then
      write(111,*)"OK------22.0"
      close(111)
      end if
!     COMPENSATE THE LOSS OF GROSS MASS

    LOGMSG()
         CALL ADDPS
      if (mytid == 0) then
      write(111,*)"OK------23.0"
      close(111)
      end if


!
!     MONITOR THE MODEL INTEGRATION

!         CALL ACCUMM

!*********************************************************************
!      ACCUMULATE SOME VARIABLES IN CARBON MODEL
!*********************************************************************

#ifdef COUP
         call shr_sys_flush(6)
#endif
!
#ifdef COUP
    LOGMSG()
         call flux_cpl
#endif
      if (mytid == 0) then
      write(111,*)"OK------24.0"
      close(111)
      end if

    LOGMSG()
!   CALL SSAVEINS
      if (mytid == 0) then
      write(111,*)"OK------25.0"
      close(111)
      end if

!     ACCUMULATE SOME VARIABLES FOR MONTHLY OUTPUT

    CALL ACCUMM
      if (mytid == 0) then
      write(111,*)"OK------26.0"
      close(111)
      end if

    if (iday==imd) then
!     CALL SSAVEMON
    endif
      if (mytid == 0) then
      write(111,*)"OK------27.0"
      close(111)
      end if
!----------------------------------------------------------------------
!   call msg_pass('send')
!----------------------------------------------------------------------
    LOGMSG()
    call licom_export_mct(o2x_o)
    if ( iday == 5) stop
    LOGMSG()
      if (mytid == 0) then
      write(111,*)"OK------28.0"
      close(111)
      end if

  end subroutine licom_run_mct


  subroutine licom_final_mct

    call cleanup
    LOGMSG()

  end subroutine licom_final_mct


!*************************************************************************
!ROUTINE: licom_import_mct
!         lihuimin 2012.7.16
!
  subroutine licom_import_mct(x2o_o)
    type(mct_aVect)   , intent(inout) :: x2o_o
    ! local
    integer :: j_begin, j_end, iblock
    type (block) :: this_block          ! block information for current block

     n=0
    do iblock = 1, nblocks_clinic
        this_block = get_block(blocks_clinic(iblock),iblock)
        do j=this_block%jb, this_block%je
        do i=this_block%ib, this_block%ie
           n=n+1
           !--- states ---
           ifrac(i,j,iblock) = x2o_o%rAttr(index_x2o_Si_ifrac,n) ! ice fraction
           patm (i,j,iblock) = x2o_o%rAttr(index_x2o_Sa_pslv,n)  ! sea level pressure index_x2o_Sa_pslv 
           !--- fluxes ---
           taux (i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_taux,n)  ! surface stress, zonal
           tauy (i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_tauy,n)  ! surface stress, merid
           netsw(i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_swnet,n) ! net sw rad
           sen  (i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_sen,n)   ! sensible
           lwup (i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_lwup,n)  ! long-wave up
           lwdn (i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_lwdn,n)  ! long-wave down
           melth(i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_melth,n) ! melt heat
           salt (i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_salt,n)  ! salinity flux
           prec (i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_prec,n)  !index_x2o_Foxx_prec 
           evap (i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_evap,n)  ! evaporation
           meltw(i,j,iblock) = x2o_o%rAttr(index_x2o_Foxx_meltw,n) ! melt water
           roff (i,j,iblock) = x2o_o%rAttr(index_x2o_Forr_roff,n)  ! runoff
           duu10n(i,j,iblock) = x2o_o%rAttr(index_x2o_So_duu10n,n)  ! 10m wind speed squared
        end do
        end do
    end do

        lat1= LATVAP*evap ! latent (derive from evap)
!       write(*,*) 'n= import',n

  end subroutine licom_import_mct



!*************************************************************************
!ROUTINE: licom_export_mct
!         lihuimin 2012.7.16
!
  subroutine licom_export_mct(o2x_o)
    type(mct_aVect)   , intent(inout) :: o2x_o
    ! local
    integer :: j_begin,j_end, iblock
    type (block) :: this_block          ! block information for current block
    real (r8) :: ek0

       n=0
! lihuimin, 2012.8.7, consider jst_global
  do iblock =1 , nblocks_clinic
        this_block = get_block(blocks_clinic(iblock),iblock)
        do j=this_block%jb, this_block%je
        do i=this_block%ib, this_block%ie
          n=n+1
          o2x_o%rAttr(index_o2x_So_t,n)    = T_cpl   (i,j,iblock) ! temperature
          o2x_o%rAttr(index_o2x_So_s,n)    = S_cpl   (i,j,iblock) ! salinity
          o2x_o%rAttr(index_o2x_So_u,n)    = U_cpl   (i,j,iblock) ! velocity, zonal
          o2x_o%rAttr(index_o2x_So_v,n)    = V_cpl   (i,j,iblock) ! velocity, meridional
          o2x_o%rAttr(index_o2x_So_dhdx,n) = dhdx(i,j,iblock) ! surface slope, zonal
          o2x_o%rAttr(index_o2x_So_dhdy,n) = dhdy(i,j,iblock) ! surface slope, meridional
          o2x_o%rAttr(index_o2x_Fioo_q,n)    = q   (i,j,iblock) ! heat of fusion xor melt pot
#ifdef USE_OCN_CARBON
!         buffs(n,cpl_fields_o2c_co2)    = co2_cpl(i,j) ! state: air-sea CO2 of ocean  ~ mol
#endif
       enddo
       enddo
   end do
!       write(*,*) 'n= export',n
  end subroutine licom_export_mct


!*************************************************************************
!ROUTINE: licom_SetGSMap_mct
!         lihuimin 2012.7.16
!
  subroutine licom_SetGSMap_mct(mpicom_ocn, OCNID, gsMap_ocn)

    implicit none
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: OCNID
    type(mct_gsMap), intent(inout) :: gsMap_ocn

    integer,allocatable :: &
      gindex(:)

    integer (kind(1)) ::   &
      i,j, k, n, iblock, &
      lsize, gsize,   &
      ier

    type (block) :: this_block          ! block information for current block

!-----------------------------------------------------------------------
!  Build the LICOM grid numbering for MCT
!  NOTE:  Numbering scheme is: West to East and South to North starting
!  at the south pole.  Should be the same as what's used in SCRIP
!-----------------------------------------------------------------------

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
       enddo
       enddo
    enddo
    lsize = n

! not correct for padding, use "n" above
!    lsize = block_size_x*block_size_y*nblocks_clinic
    gsize = imt_global*jmt_global
    allocate(gindex(lsize),stat=ier)

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          gindex(n) = (jmt_global- this_block%j_glob(j))*(imt_global) + this_block%i_glob(i)
       enddo
       enddo
    enddo


    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )

    deallocate(gindex)

  end subroutine licom_SetGSMap_mct

!*************************************************************************
!ROUTINE: licom_domain_mct
!         lihuimin 2012.7.16
!
  subroutine licom_domain_mct(lsize, gsMap_o, dom_o)

    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o
    ! local
    integer :: j_begin,j_end
    integer, pointer :: &
      idata(:)

    real(r8), pointer :: &
      data(:)

    integer (kind(1)) ::    i,j, k, n, iblock, ier

    type (block) :: this_block          ! block information for current block

!
  integer            :: fid    ! nc domain file ID
  integer            :: dimid  ! nc dimension id
  integer            :: vid    ! nc variable ID
  integer            :: rcode  ! nc return code
  integer            :: ntim   ! temporary
 
!
!     Define Variables.
      integer*4   :: ncid, iret
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
      character (len=18) :: fname
      real(r8), allocatable :: data1(:,:),data2(:,:), temp1(:,:,:),temp2(:,:,:),temp3(:,:,:)
      INTEGER :: NMFF

!
  if (trim(horiz_grid_opt) == 'lat_lon') then
  if (mytid==0 ) then
     call wrap_open ('domain_licom.nc', NF_NOWRITE, fid)

     ! obtain dimensions
     call wrap_inq_dimid  (fid, 'ni' , dimid)
     call wrap_inq_dimlen (fid, dimid, nx   )
     call wrap_inq_dimid  (fid, 'nj' , dimid)
     call wrap_inq_dimlen (fid, dimid, ny   )
     call wrap_inq_dimid  (fid, 'nv' , dimid)
     call wrap_inq_dimlen (fid, dimid, nv   )

     if (nx/=imt_global) then
        write(6,*)"nx=",nx,"imt=",imt_global
        stop
     end if
  end if


  call shr_mpi_bcast(nx,mpi_comm_ocn,"nx")
  call shr_mpi_bcast(ny,mpi_comm_ocn,"ny")

     ! obtain grid variables
     allocate(data1(nx,ny), data2(nx,ny))
     allocate(temp1(imt,jmt,max_blocks_clinic),temp2(imt,jmt,max_blocks_clinic),temp3(imt,jmt,max_blocks_clinic))

  if (mytid==0 ) then
     call wrap_inq_varid(fid, 'xc' ,vid)
     call wrap_get_var_realx(fid, vid, data1)
     do j=1, ny
     do i=1, nx
        data2(i,j) = data1(i,1+ny-j)
     end do
     end do
  end if
      call scatter_global(temp1, data2, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

  if (mytid==0 ) then
     call wrap_inq_varid(fid, 'yc' ,vid)
     call wrap_get_var_realx(fid, vid, data1)
     do j=1, ny
     do i=1, nx
        data2(i,j) = data1(i,1+ny-j)
     end do
     end do
  end if
      call scatter_global(temp2, data2, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
  if (mytid==0 ) then
     call wrap_inq_varid(fid, 'area' ,vid)
     call wrap_get_var_realx(fid, vid, data1)
     do j=1, ny
     do i=1, nx
        data2(i,j) = data1(i,1+ny-j)
     end do
     end do
     call wrap_close(fid)
  end if
      call scatter_global(temp3, data2, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

     deallocate(data1,data2)
  end if
!
    call mct_gGrid_init( GGrid=dom_o, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_o%data)
    allocate(data(lsize))


    call mct_gsMap_orderedPoints(gsMap_o, mytid, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)



    data(:) = shr_const_spval !-9999.0_R8 !linpf 2012Jul27
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize)
    data(:) = shr_const_spval !0.0_R8  !linpf 2012Jul27
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize)
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize)


      write(6,*)"domain_mct  special, id =  ", mytid,"   imt=",imt,"   jmt=",jmt,"   nx=",nx,"   ny=",ny,"lsize=",lsize

!
!-------------------------------------------------------------------
!
! Fill in correct values for domain components
!
!-------------------------------------------------------------------

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          if (trim(horiz_grid_opt) == 'lat_lon') then
            data(n) = temp1(i,j,iblock)
          else
            data(n) = TLON(i,j,iblock)/DEGtoRAD
          end if
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lon",data,lsize)
!
    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          if (trim(horiz_grid_opt) == 'lat_lon') then
            data(n) = temp2(i,j,iblock)
          else
            data(n) = TLAT(i,j,iblock)/DEGtoRAD
          end if
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lat",data,lsize)

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          if (trim(horiz_grid_opt) == 'lat_lon') then
            data(n) = temp3(i,j,iblock)
          else
            data(n) = TAREA(i,j,iblock)/(radius*radius)
          end if
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"area",data,lsize)
!
     n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
!         data(n) = float(mask_r(i,j,iblock))
          data(n) = float(KMT(i,j,iblock))
          if (data(n) > 1.0_r8) data(n) = 1.0_r8
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"mask",data,lsize)
    call mct_gGrid_importRattr(dom_o,"frac",data,lsize)

    deallocate(data)
    deallocate(idata)
    if (trim(horiz_grid_opt) == 'lat_lon') deallocate (temp1, temp2, temp3)

  end subroutine licom_domain_mct

!*****************************************************************************
!ROUTINE:   CLEANUP
!           lihuimin, 2012.7.22
!
  subroutine cleanup

    ! allocated in inirun.F90
    deallocate(t_cpl, s_cpl, u_cpl, v_cpl, dhdx, dhdy, Q)
    deallocate(taux, tauy, netsw, lat1, sen, lwup, lwdn, melth, salt, prec, evap, meltw, roff, ifrac, patm, duu10n)

    ! allocated in inirun.F90
    deallocate(h0, u, v, at)

    ! allocated in rdriver.F90
    deallocate(su3,sv3,psa3,tsa3,qar3,uva3,swv3,cld3,sss3,sst3 ,nswv3,dqdt3,chloro3)
    deallocate(seaice3,runoff3)
    deallocate(wspd3,wspdu3,wspdv3,lwv3,rain3,snow3) 

  end subroutine cleanup



 end module licom_comp_mct
