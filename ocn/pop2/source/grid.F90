!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module grid

!BOP
! !MODULE: grid
!
! !DESCRIPTION:
!  This module contains grid info and routines for setting up the
!  POP grid quantities.
!
! !REVISION HISTORY:
!  SVN:$Id: grid.F90 17212 2009-07-20 23:01:42Z njn01 $

! !USES:

   use POP_KindsMod
   use precision_mod
   use param_mod
   use LICOM_Error_mod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod

   use blocks
   use distribution
   use domain
   use broadcast
   use gather_scatter
   use constant_mod
   use msg_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: init_grid1,     &
              init_grid2,     &
              read_horiz_grid,&
              read_topography
!             tgrid_to_ugrid, &
!             ugrid_to_tgrid, &

! !PUBLIC DATA MEMBERS:

   real (POP_r8), public :: &
      area_u, area_t       ,&! total ocean area of U,T cells
      volume_u, volume_t   ,&! total ocean volume of U,T cells
      volume_t_marg        ,&! volume of marginal seas (T cells)
      area_t_marg          ,&! area of marginal seas (T cells)
      uarea_equator          ! area of equatorial cell

   real (POP_r8), dimension(km), public :: &
      area_t_k             ,&! total ocean area (T cells) at each dpth
      volume_t_k           ,&! total ocean volume (T cells) at each dpth
      volume_t_marg_k        ! tot marginal seas vol (T cells) at each dpth

   integer (POP_i4), public :: &
      sfc_layer_type,       &! choice for type of surface layer
      kmt_kmin, 	    &! minimum allowed non-zero KMT value
      n_topo_smooth 	     ! number of topo smoothing passes

   integer (POP_i4), parameter, public :: &
      sfc_layer_varthick = 1,  &! variable thickness surface layer
      sfc_layer_rigid    = 2,  &! rigid lid surface layer
      sfc_layer_oldfree  = 3    ! old free surface form

   logical (POP_logical), public ::    &
      partial_bottom_cells   ! flag for partial bottom cells

   real (POP_r8), dimension(:,:), allocatable, public :: &
      BATH_G           ! Observed ocean bathymetry mapped to global T grid
                       ! for use in computing KMT internally

   integer (POP_i4), dimension(:,:), allocatable, public :: &
      KMT_G            ! k index of deepest grid cell on global T grid
                       ! for use in performing work distribution

!-----------------------------------------------------------------------
!
!  grid information for all local blocks
!  the local blocks are by default in baroclinic distribution
!
!-----------------------------------------------------------------------

   !*** dimension(1:km)

   real (POP_r8), dimension(km), public :: &
      dz                ,&! thickness of layer k
      c2dz              ,&! 2*dz
      dzr, dz2r         ,&! reciprocals of dz, c2dz
      zt                ,&! vert dist from sfc to midpoint of layer
      zw                  ! vert dist from sfc to bottom of layer

   !*** dimension(0:km)

   real (POP_r8), dimension(0:km), public :: &
      dzw, dzwr          ! midpoint of k to midpoint of k+1
                         !   and its reciprocal

   !*** geometric 2d arrays

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic), public :: &
      DXU, DYU            ,&! {x,y} spacing centered at U points
      DXT, DYT            ,&! {x,y} spacing centered at T points
      DXUR, DYUR          ,&! reciprocals of DXU, DYU
      DXTR, DYTR          ,&! reciprocals of DXT, DYT
      HTN, HTE            ,&! cell widths on {N,E} sides of T cell
      HUS, HUW            ,&! cell widths on {S,W} sides of U cell
      ULAT, ULON          ,&! {latitude,longitude} of U points
      TLAT, TLON          ,&! {latitude,longitude} of T points
      ANGLE, ANGLET       ,&! angle grid makes with latitude line
      FCOR, FCORT         ,&! coriolis parameter at U,T points
      UAREA, TAREA        ,&! area of U,T cells
      UAREA_R, TAREA_R    ,&! reciprocal of area of U,T cells
      HT, HU, HUR           ! ocean depth at T,U points

   !*** 3d depth fields for partial bottom cells

   real (POP_r8), dimension(:,:,:,:), allocatable, public :: &
      DZU, DZT               ! thickness of U,T cell for pbc

   !*** 2d landmasks

   integer (POP_i4), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      KMT            ,&! k index of deepest grid cell on T grid
      KMU            ,&! k index of deepest grid cell on U grid
      KMTOLD           ! KMT field before smoothing

   logical (POP_logical), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      CALCT          ,&! flag=true if point is an ocean point
      CALCU            !   at the surface

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic), public :: &
      RCALCT         ,&! real equiv of CALCT,U to use as more
      RCALCU           !   efficient multiplicative mask

   integer (POP_i4), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      KMTN,KMTS,KMTE,KMTW   ,&! KMT field at neighbor points
      KMUN,KMUS,KMUE,KMUW     ! KMU field at neighbor points

   integer (POP_i4), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      KMTEE,KMTNN      ! KMT field 2 cells away for upwind stencil
                       ! allocated and computed in advection module

   integer (POP_i4), dimension(:,:,:), allocatable, public :: &
      REGION_MASK      ! mask defining regions, marginal seas

!-----------------------------------------------------------------------
!
!     define types used with region masks and marginal seas balancing
!
!-----------------------------------------------------------------------
   integer (POP_i4), parameter, public :: &
      max_regions =   15, &              ! max no. ocean regions
      max_ms      =    7                 ! max no. marginal seas
 
   integer (POP_i4), public :: &
      num_regions, &
      num_ms

   type, public :: ms_bal
      real    (POP_r8)         :: lat         ! transport latitude
      real    (POP_r8)         :: lon         ! transport longitude
      real    (POP_r8)         :: area        ! total distribution area
      real    (POP_r8)         :: transport   ! total excess/deficit (E+P+M+R)
      integer (POP_i4)   :: mask_index  ! index of m-s balancing mask
   end type ms_bal

   type, public :: regions                ! region-mask info
      integer   (POP_i4) :: number 
      character (char_len) :: name          
      logical   (POP_logical) :: marginal_sea
      real      (POP_r8      ) :: area
      real      (POP_r8      ) :: volume
      type      (ms_bal)   :: ms_bal
   end type regions       

   type (regions),dimension(max_regions), public :: region_info

   integer (POP_i4),public ::       &
      nocean_u, nocean_t,      &! num of ocean U,T points
      nsurface_u, nsurface_t    ! num of ocean U,T points at surface

!#################### temporary kludge for overflows ####################
   character (char_len), public ::  &
      topography_filename    ! copy of input topography filename 
!########################################################################

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module private data
!
!-----------------------------------------------------------------------

   !*** geometric scalars

   integer (POP_i4) ::       &
      jeq                       ! j index of equatorial cell

   logical (POP_logical) ::  &
      flat_bottom,        &! flag for flat-bottom topography
      lremove_points       ! flag for removing isolated points

   real (POP_r8), dimension(:,:), allocatable :: &
      ULAT_G, ULON_G        ! {latitude,longitude} of U points
                            ! in global-sized array

!-----------------------------------------------------------------------
!
!     area-weighted averaging coefficients
!     AT{0,S,W,SW} = {central,s,w,sw} coefficients for area-weighted
!       averaging of four U points surrounding a T point
!     AU{0,N,E,NE} = {central,n,e,ne} coefficients for area-weighted
!       averaging of four T points surrounding a U point
!
!-----------------------------------------------------------------------

   real (POP_r8), dimension (nx_block,ny_block,max_blocks_clinic) :: &
      AT0,ATS,ATW,ATSW,AU0,AUN,AUE,AUNE

!-----------------------------------------------------------------------
!
!  variables which are shared between init_grid1,init_grid2
!
!-----------------------------------------------------------------------

   character (char_len) ::  &
      horiz_grid_opt,       &! horizontal grid option
      vert_grid_opt,        &! vertical grid option
      sfc_layer_opt,        &! choice for surface layer type
      topography_opt,       &! topography (KMT) option
      horiz_grid_file,      &! input file for reading horiz grid info
      vert_grid_file,       &! input file for reading horiz grid info
      topography_file,      &! input file for reading horiz grid info
      bathymetry_file,      &! input file for reading horiz grid info
      region_mask_file,     &! input file for region mask
      region_info_file,     &! input file with region identification
      bottom_cell_file,     &! input file for thickness of pbc
      topography_outfile     ! output file for writing horiz grid info

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_grid1
! !INTERFACE:

 subroutine init_grid1

! !DESCRIPTION:
!  Initializes only grid quantities necessary for completing
!  decomposition setup (ULAT, ULON, KMT).
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   namelist /grid_nml/horiz_grid_opt, vert_grid_opt, topography_opt,   &
                      horiz_grid_file, vert_grid_file, topography_file,&
		      topography_outfile,bathymetry_file,             &
                      n_topo_smooth, flat_bottom, lremove_points,      &
                      region_mask_file, region_info_file,sfc_layer_opt,&
                      partial_bottom_cells, bottom_cell_file, kmt_kmin

   integer (POP_i4) :: &
      nml_error           ! namelist i/o error flag

!-----------------------------------------------------------------------
!
!  read input namelist for grid setup options
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  get global ULAT,ULON
!
!-----------------------------------------------------------------------

      call broadcast_scalar(horiz_grid_file, master_task)
      call read_horiz_grid(horiz_grid_file,.true.)

!-----------------------------------------------------------------------
!  set up topography by getting global KMT field (used for
!  creating a load balanced block distribution).
!
!-----------------------------------------------------------------------

      call broadcast_scalar(topography_file, master_task)
      call read_topography(topography_file,.true.)

 end subroutine init_grid1

!***********************************************************************
!BOP
! !IROUTINE: init_grid2
! !INTERFACE:

 subroutine init_grid2

! !DESCRIPTION:
!  Initializes all grid quantities
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (POP_i4) :: errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      reclength,ioerr,nu,i,j,k,n,iblock,    &! dummy loop index variables
      range_count         ! counter for angle out of range

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      WORK                 ! local temp space

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      DZBC                 ! thickness of bottom T cell for pbc

   character (*), parameter ::  &! output formats
      vgrid_fmt1 = "(3x,' k ',3x,'Thickness (cm)',3x,' Depth (cm) ')", &
      vgrid_fmt2 = "(3x,'---',3x,'--------------',3x,'------------')", &
      vgrid_fmt3 = "(3x,i3,4x,1pe12.5,4x,1pe12.5)"                   , &
      topo_fmt1  = "(' # surface (T,U) points',2x,i10,2x,i10)"       , &
      topo_fmt2  = "(' # ocean   (T,U) points',2x,i10,2x,i10)"       , &
      topo_fmt3  = "(' T-area,   U-area   (km^2)',2(2x,1pe23.15))"   , &
      topo_fmt4  = "(' T-volume, U-volume (km^3)',2(2x,1pe23.15))"

   type (block) :: &
      this_block  ! block info for current block

   real (POP_r8) ::         &
      angle_0, angle_w, &! temporaries for computing angle at T points
      angle_s, angle_sw 

!-----------------------------------------------------------------------
!
!  output grid setup options to log file
!
!-----------------------------------------------------------------------

   errorCode = 0

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a13)') ' Grid options'
      write(stdout,blank_fmt)
   endif

!-----------------------------------------------------------------------
!
!  set up horizontal grid
!
!-----------------------------------------------------------------------

      call broadcast_scalar(horiz_grid_file, master_task)
      call read_horiz_grid(horiz_grid_file,.false.)


!-----------------------------------------------------------------------
!
!  if boundaries are closed, extend physical domain values into ghost
!  cells
!  compute other derived quantities like areas and reciprocals
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,i,j,this_block)
   do iblock=1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      if (this_block%i_glob(1) == 0) then ! closed western bndy
         do j=1,ny_block
         do i=1,this_block%ib-1
            DXU(i,j,iblock) = DXU(this_block%ib,j,iblock)
            DYU(i,j,iblock) = DYU(this_block%ib,j,iblock)
            DXT(i,j,iblock) = DXT(this_block%ib,j,iblock)
            DYT(i,j,iblock) = DYT(this_block%ib,j,iblock)
         end do
         end do
      endif

      if (this_block%i_glob(this_block%ie+1) == 0) then ! closed east bndy
         do j=1,ny_block
         do i=this_block%ie+1,nx_block
            DXU(i,j,iblock) = DXU(this_block%ie,j,iblock)
            DYU(i,j,iblock) = DYU(this_block%ie,j,iblock)
            DXT(i,j,iblock) = DXT(this_block%ie,j,iblock)
            DYT(i,j,iblock) = DYT(this_block%ie,j,iblock)
         end do
         end do
      endif

      if (this_block%j_glob(1) == 0) then ! closed southern bndy
         do j=1,this_block%jb-1
         do i=1,nx_block
            DXU(i,j,iblock) = DXU(i,this_block%jb,iblock)
            DYU(i,j,iblock) = DYU(i,this_block%jb,iblock)
            DXT(i,j,iblock) = DXT(i,this_block%jb,iblock)
            DYT(i,j,iblock) = DYT(i,this_block%jb,iblock)
         end do
         end do
      endif

      if (this_block%j_glob(this_block%je+1) == 0) then ! closed north bndy
         do j=this_block%je+1,ny_block
         do i=1,nx_block
            DXU(i,j,iblock) = DXU(i,this_block%je,iblock)
            DYU(i,j,iblock) = DYU(i,this_block%je,iblock)
            DXT(i,j,iblock) = DXT(i,this_block%je,iblock)
            DYT(i,j,iblock) = DYT(i,this_block%je,iblock)
         end do
         end do
      endif

      DXUR(:,:,iblock) = c1/DXU(:,:,iblock)
      DYUR(:,:,iblock) = c1/DYU(:,:,iblock)

      UAREA(:,:,iblock) = DXU(:,:,iblock)*DYU(:,:,iblock)
      UAREA_R(:,:,iblock) = c1/UAREA(:,:,iblock)

      DXTR(:,:,iblock) = c1/DXT(:,:,iblock)
      DYTR(:,:,iblock) = c1/DYT(:,:,iblock)

      TAREA(:,:,iblock) = DXT(:,:,iblock)*DYT(:,:,iblock)
      TAREA_R(:,:,iblock) = c1/TAREA(:,:,iblock)

   end do
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  calculate stencil coefficients for area-averaging
!
!-----------------------------------------------------------------------

!  call cf_area_avg  ! coefficients for area-weighted averages

!-----------------------------------------------------------------------
!
!  calculate lat/lon of T points and calculate ANGLET from ANGLE
!
!-----------------------------------------------------------------------

!  call calc_tpoints(errorCode)

!  if (errorCode /= 0) then
!     call POP_ErrorSet(errorCode, &
!        'init_grid2: error in calc_tpoints')
!     return
!  endif

!  !***
!  !*** first, ensure that -pi <= ANGLE <= pi
!  !***

!  range_count = global_count ((ANGLE < - pi .or. ANGLE > pi), &
!                              distrb_clinic, field_loc_NEcorner)

!  if (range_count > 0) call exit_POP(sigAbort, &
!                       'ERROR: ANGLE is outside its expected range')

!  !***
!  !*** compute ANGLE on T-grid
!  !***

!  ANGLET = c0


!  !$OMP PARALLEL DO PRIVATE (n,i,j,angle_0,angle_w,angle_s,angle_sw, &
!  !$OMP                      this_block)

!  do n=1,nblocks_clinic
!     this_block = get_block(blocks_clinic(n),n)

!     do j=this_block%jb,this_block%je
!        do i=this_block%ib,this_block%ie

!           angle_0  = ANGLE(i,  j  ,n)
!           angle_w  = ANGLE(i-1,j  ,n)
!           angle_s  = ANGLE(i,  j-1,n)
!           angle_sw = ANGLE(i-1,j-1,n)

!           if ( angle_0 < c0 ) then
!              if ( abs(angle_w -angle_0) > pi ) &
!                 angle_w  = angle_w  - pi2 
!              if ( abs(angle_s -angle_0) > pi ) &
!                 angle_s  = angle_s  - pi2 
!              if ( abs(angle_sw-angle_0) > pi ) &
!                 angle_sw = angle_sw - pi2 
!           endif

!           ANGLET(i,j,n) =  angle_0 *AT0 (i,j,n) + &
!                            angle_w *ATW (i,j,n) + &
!                            angle_s *ATS (i,j,n) + &
!                            angle_sw*ATSW(i,j,n)

!        enddo

!        !***
!        !*** set ANGLET to zero for all of (global) j=1 row 
!        !*** (bottom row of ANGLET is not used, but is written to file)
!        !***

!        if (this_block%j_glob(j) == 1) ANGLET(:,j,n) = c0

!     enddo
!  enddo
!  !$OMP END PARALLEL DO

!  call POP_HaloUpdate(ANGLET, POP_haloClinic, POP_gridHorzLocCenter, &
!                              POP_fieldKindAngle, errorCode,         &
!                              fillValue = 0.0_POP_r8)

!  if (errorCode /= 0) then
!     call POP_ErrorSet(errorCode, &
!        'init_grid2: error updating angleT halo')
!     return
!  endif

!-----------------------------------------------------------------------
!
!  set up vertical grid
!
!-----------------------------------------------------------------------

      if (my_task == master_task) then
         write(stdout,*) ' Reading vertical grid from file:', &
                         trim(vert_grid_file)
      endif
      call broadcast_scalar(vert_grid_file, master_task)
      call read_vert_grid(vert_grid_file)

!  !***
!  !*** calculate other vertical grid quantities
!  !***

!  dzw(0)  = p5*dz(1)
!  dzw(km) = p5*dz(km)
!  dzwr(0) = c1/dzw(0)
!  zw(1) = dz(1)
!  zt(1) = dzw(0)

!  do k = 1,km-1
!     dzw(k) = p5*(dz(k) + dz(k+1))
!     zw(k+1) = zw(k) + dz(k+1)
!     zt(k+1) = zt(k) + dzw(k)
!  enddo

!  do k = 1,km
!     c2dz(k) = c2*dz(k)
!     dzr(k)  = c1/dz(k)
!     dz2r(k) = c1/c2dz(k)
!     dzwr(k) = c1/dzw(k)
!  enddo

!  if (my_task == master_task) then
!     write(stdout,blank_fmt)
!     write(stdout,'(a15)') ' Vertical grid:'
!     write(stdout,vgrid_fmt1)
!     write(stdout,vgrid_fmt2)
!     do k=1,km
!        write(stdout,vgrid_fmt3) k,dz(k),zt(k)
!     end do
!     write(stdout,blank_fmt)
!  endif

!-----------------------------------------------------------------------
!
!  set up topography
!
!-----------------------------------------------------------------------

      if (my_task == master_task) write(stdout,'(a30,a)') &
         ' Reading topography from file:', trim(topography_file)
      call broadcast_scalar(topography_file, master_task)
      call read_topography(topography_file,.false.)

!-----------------------------------------------------------------------
!
!  landmasks
!
!-----------------------------------------------------------------------

!  call landmasks

!-----------------------------------------------------------------------
!
!  calculate area, volume, # surface points, # ocean points
!
!-----------------------------------------------------------------------

!  area_t   = global_sum(TAREA, distrb_clinic, field_loc_center, RCALCT)
!  area_u   = global_sum(UAREA, distrb_clinic, field_loc_NEcorner, RCALCU)
!  WORK = TAREA*HT
!  volume_t = global_sum(WORK, distrb_clinic, field_loc_center, RCALCT)
!  WORK = UAREA*HU
!  volume_u = global_sum(WORK, distrb_clinic, field_loc_NEcorner, RCALCU)
!  area_t_k(1) = area_t
!  volume_t_k(1) = global_sum(TAREA*dz(1), distrb_clinic, &
!                             field_loc_center, RCALCT)
!  do k=2,km
!     WORK = merge(TAREA, c0, k <= KMT)
!     area_t_k(k) = global_sum(WORK, distrb_clinic, field_loc_center)
!     WORK = merge(TAREA*dz(k), c0, k <= KMT)
!     volume_t_k(k) = global_sum(WORK, distrb_clinic, field_loc_center)
!  end do

!  nsurface_t = global_count(RCALCT, distrb_clinic, field_loc_center)
!  nsurface_u = global_count(RCALCU, distrb_clinic, field_loc_NEcorner)
!  nocean_t = global_sum(KMT, distrb_clinic, field_loc_center)
!  nocean_u = global_sum(KMU, distrb_clinic, field_loc_NEcorner)

!  if (my_task == master_task) then
!     write(stdout,blank_fmt)
!     write(stdout,topo_fmt1) nsurface_t, nsurface_u
!     write(stdout,topo_fmt2) nocean_t, nocean_u
!     write(stdout,topo_fmt3) area_t*1.0e-10_POP_r8, &
!                             area_u*1.0e-10_POP_r8
!     write(stdout,topo_fmt4) volume_t*1.0e-15_POP_r8, &
!                             volume_u*1.0e-15_POP_r8
!     write(stdout,blank_fmt)
!  endif

!
!  compute coriolis parameter 2*omega*sin(true_latitude)
!
!-----------------------------------------------------------------------

!  FCOR  = c2*omega*sin(ULAT)    ! at u-points
!  FCORT = c2*omega*sin(TLAT)    ! at t-points

!-----------------------------------------------------------------------
!EOC

!call flushm (stdout)

 end subroutine init_grid2

!***********************************************************************
! !IROUTINE: read_horiz_grid
! !INTERFACE:

 subroutine read_horiz_grid(horiz_grid_file, latlon_only)

! !DESCRIPTION:
!  Reads horizontal grid information from input grid file
!
! !REVISION HISTORY:
!  same as module

!  !INPUT PARAMETERS:

   character (*), intent(in) :: &
      horiz_grid_file     ! filename of file containing grid data

   logical (POP_logical), intent(in) :: &
      latlon_only       ! flag requesting only ULAT, ULON

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,iblock        ,&! loop counters
      ip1, im1, jp1, jm1,&! shift indexes
      nu                ,&! i/o unit number
      ioerr             ,&! i/o error flag
      reclength           ! record length

   type (block) :: &
      this_block

!-----------------------------------------------------------------------
!
!  if only lat,lon are requested, read only these
!
!-----------------------------------------------------------------------


      allocate (ULAT_G(imt_global,jmt_global), &
                ULON_G(imt_global,jmt_global))



      if (my_task == master_task) then
         read(nu,rec=1,iostat=ioerr) ULAT_G
         read(nu,rec=2,iostat=ioerr) ULON_G
      endif

      call scatter_global(ULAT, ULAT_G, master_task, distrb_clinic, &
                          field_loc_NEcorner, field_type_scalar)
      call scatter_global(ULON, ULON_G, master_task, distrb_clinic, &
                          field_loc_NEcorner, field_type_scalar)

      if (my_task == master_task) then
         read(nu,rec=3,iostat=ioerr) ULAT_G  ! holds HTN
      endif

      call scatter_global(HTN, ULAT_G, master_task, distrb_clinic, &
                          field_loc_Nface, field_type_scalar)

      do j=1,jmt_global
      do i=1,imt_global
         ip1 = i+1
         if (i == imt_global) ip1 = 1 ! assume cyclic. non-cyclic
                                     ! will be handled during scatter
         !DXU
         ULON_G(i,j) = p5*(ULAT_G(i,j) + ULAT_G(ip1,j))
      end do
      end do
      call scatter_global(DXU, ULON_G, master_task, distrb_clinic, &
                          field_loc_NEcorner, field_type_scalar)

      do j=1,jmt_global
         jm1 = j-1
         if (j == 1) jm1 = jmt_global ! assume cyclic. non-cyclic
                                     ! will be handled during scatter
         do i=1,imt_global

            !DXT = p5(HTN(i,j)+HTN(i,j-1))
            ULON_G(i,j) = p5*(ULAT_G(i,j) + ULAT_G(i,jm1))
         end do
      end do
      call scatter_global(DXT, ULON_G, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

      if (my_task == master_task) then
         read(nu,rec=4,iostat=ioerr) ULAT_G  ! holds HTE
      endif

      call scatter_global(HTE, ULAT_G, master_task, distrb_clinic, &
                          field_loc_Eface, field_type_scalar)

      do j=1,jmt_global
      do i=1,imt_global
         im1 = i-1
         if (i == 1) im1 = imt_global ! assume cyclic. non-cyclic
                                     ! will be handled during scatter
         !DYT
         ULON_G(i,j) = p5*(ULAT_G(i,j) + ULAT_G(im1,j))
      end do
      end do
      call scatter_global(DYT, ULON_G, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

      do j=1,jmt_global
         jp1 = j+1
         if (j == jmt_global) jp1 = 1 ! assume cyclic. non-cyclic
                                     ! will be handled during scatter
         do i=1,imt_global

            !DYU = p5(HTE(i,j)+HTN(i,j+1))
            ULON_G(i,j) = p5*(ULAT_G(i,j) + ULAT_G(i,jp1))
         end do
      end do

!-----------------------------------------------------------------------
 
!   Add tripole-grid correction (KL and GD)
!
!-----------------------------------------------------------------------

      if (ltripole_grid) then
         j=jmt_global
         do i=1,imt_global
            ULON_G(i,j) = ULAT_G(i,j)
         end do
      endif

      call scatter_global(DYU, ULON_G, master_task, distrb_clinic, &
                          field_loc_NEcorner, field_type_scalar)

      if (my_task == master_task) then
         read(nu,rec=5,iostat=ioerr) ULAT_G
         read(nu,rec=6,iostat=ioerr) ULON_G
      endif

      call scatter_global(HUS, ULAT_G, master_task, distrb_clinic, &
                          field_loc_Eface, field_type_scalar)
      call scatter_global(HUW, ULON_G, master_task, distrb_clinic, &
                          field_loc_Nface, field_type_scalar)

      if (my_task == master_task) then
         read(nu,rec=7,iostat=ioerr) ULAT_G
         close(nu)
      endif

      call scatter_global(ANGLE, ULAT_G, master_task, distrb_clinic, &
                          field_loc_NEcorner, field_type_angle)
      deallocate(ULAT_G,ULON_G)

      where (HTN <= c0) HTN = c1
      where (HTE <= c0) HTE = c1
      where (HUS <= c0) HUS = c1
      where (HUW <= c0) HUW = c1
      where (DXU <= c0) DXU = c1
      where (DYU <= c0) DYU = c1
      where (DXT <= c0) DXT = c1
      where (DYT <= c0) DYT = c1

!-----------------------------------------------------------------------
!EOC

 end subroutine read_horiz_grid

!***********************************************************************
!BOP
! !IROUTINE: vert_grid_internal
! !INTERFACE:

!***********************************************************************
!BOP
! !IROUTINE: compute_dz
! !INTERFACE:

 subroutine compute_dz(depth,zlength,dz_sfc,dz_deep)

! !DESCRIPTION:
!  Computes a thickness profile and total depth given the
!  parameters for the thickness function
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), intent(in) :: &
      zlength,          &! gaussian parameter for thickness func
      dz_sfc,           &! thickness of surface layer
      dz_deep            ! thickness of deep ocean layers

! !OUTPUT PARAMETERS:

   real (POP_r8), intent(out) :: &
      depth               ! depth based on integrated thicknesses

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: k

!-----------------------------------------------------------------------

!  depth = c0

!  do k=1,km
!     dz(k) = dz_deep - (dz_deep - dz_sfc)*exp(-(depth/zlength)**2)
!     depth = depth + dz(k)
!  end do

!-----------------------------------------------------------------------
!EOC

 end subroutine compute_dz

!***********************************************************************
!BOP
! !IROUTINE: read_vert_grid
! !INTERFACE:

 subroutine read_vert_grid(vert_grid_file)

! !DESCRIPTION:
!  Reads in layer thicknesses from grid input file
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (char_len), intent(in) :: &
      vert_grid_file

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      k,                 &! vertical level index
      nu,                &! i/o unit number
      ioerr               ! i/o error flag

!-----------------------------------------------------------------------
!
!  read vertical layer thickness from file
!
!-----------------------------------------------------------------------

!  call get_unit(nu)
!  if (my_task == master_task) then
!     open(nu,file=vert_grid_file,status='old',form='formatted', &
!             iostat=ioerr)
!  endif
!  call broadcast_scalar(ioerr, master_task)
!  if (ioerr /= 0) call exit_POP(sigAbort, &
!                                'Error opening vert_grid_file')

!  if (my_task == master_task) then
!     if (ioerr == 0) then ! successful open
!        grid_read: do k = 1,km
!           read(nu,*,iostat=ioerr) dz(k)
!           if (ioerr /= 0) exit grid_read
!        end do grid_read
!        close(nu)
!     endif
!  endif
!  call release_unit(nu)

!  call broadcast_scalar(ioerr, master_task)
!  if (ioerr /= 0) call exit_POP(sigAbort, &
!                                'Error reading vert_grid_file')

!  call broadcast_array(dz, master_task)

!-----------------------------------------------------------------------
!EOC

 end subroutine read_vert_grid

!***********************************************************************

 subroutine read_topography(topography_file,kmt_global)

! !DESCRIPTION:
!  Reads in KMT field from file
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (char_len), intent(in) :: &
      topography_file     ! input file containing KMT field

   logical(POP_logical), intent(in) :: &
      kmt_global       ! flag for generating only global KMT field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      nu                ,&! i/o unit number
      ioerr             ,&! i/o error flag
      reclength           ! record length

!-----------------------------------------------------------------------
!
!  read global KMT field (for use in setting up domain)
!
!-----------------------------------------------------------------------

   if (kmt_global) then

      allocate(KMT_G(imt_global,jmt_global))


      if (my_task == master_task) then
         read(nu, rec=1, iostat=ioerr) KMT_G
         close(nu)
      endif
      call broadcast_array(KMT_G, master_task)

!-----------------------------------------------------------------------
!
!  otherwise read KMT field from file
!
!-----------------------------------------------------------------------

   else
      call scatter_global(KMT, KMT_G, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      deallocate(KMT_G)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine read_topography

!***********************************************************************
!BOP
! !IROUTINE: landmasks
! !INTERFACE:

 subroutine landmasks

! !DESCRIPTION:
!  Calculates additional masks for land points at each depth level.
!  These include real masks for applying multiplicative masks
!  instead of logical masks and also KMT arrays for neighbor points.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  masks for surface ocean T or U points
!
!-----------------------------------------------------------------------

!  where (KMT >= 1)
!     CALCT = .true.
!     RCALCT = c1
!  elsewhere
!     CALCT = .false.
!     RCALCT = c0
!  endwhere

!  where (KMU >= 1)
!     CALCU = .true.
!     RCALCU = c1
!  elsewhere
!     CALCU = .false.
!     RCALCU = c0
!  endwhere

!-----------------------------------------------------------------------
!
!  depth level fields (KMT,KMU) to north,south,east,west
!
!-----------------------------------------------------------------------

!  KMTN = eoshift(KMT,dim=2,shift=+1)
!  KMTS = eoshift(KMT,dim=2,shift=-1)
!  KMTE = eoshift(KMT,dim=1,shift=+1)
!  KMTW = eoshift(KMT,dim=1,shift=-1)

!  KMUN = eoshift(KMU,dim=2,shift=+1)
!  KMUS = eoshift(KMU,dim=2,shift=-1)
!  KMUE = eoshift(KMU,dim=1,shift=+1)
!  KMUW = eoshift(KMU,dim=1,shift=-1)

!  KMTEE = eoshift(KMT,dim=1,shift=2)
!  KMTNN = eoshift(KMT,dim=2,shift=2)

!-----------------------------------------------------------------------
!EOC

 end subroutine landmasks

!***********************************************************************
!BOP
! !IROUTINE: area_masks
! !INTERFACE:

 subroutine area_masks(mask_filename,info_filename)

! !DESCRIPTION:
!   This subroutine reads in file with regional area mask and
!   marginal seas defined
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in) :: &
      mask_filename           ,&! name of file containing region masks
      info_filename             ! name of file containing region names

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      k, n,              &! loop counters
      nu,                &! i/o unit number
      reclength,         &! record length of file
      ioerr,             &! i/o error flag
      region              ! region counter

   real (POP_r8) :: &
      tmp_vol,  &! temporary volume
      sea_area, &! total volume of a particular sea
      sea_vol    ! total volume of a particular sea

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      WORK                ! temporary space

   integer (POP_i4), dimension(:,:), allocatable :: &
      REGION_G            ! global-sized region mask

!-----------------------------------------------------------------------
!
!  read in regional area masks, including marginal seas.  then
!  calculate related variables.
!
!-----------------------------------------------------------------------

!  allocate(REGION_MASK(nx_block,ny_block,max_blocks_clinic))

!  call get_unit(nu)
!  if (my_task == master_task) then
!     allocate(REGION_G(imt_global,jmt_global))
!     inquire (iolength=reclength) REGION_G
!     open(nu, file=mask_filename,status='old',form='unformatted', &
!              access='direct', recl=reclength, iostat=ioerr)
!  endif

!  call broadcast_scalar(ioerr, master_task)
!  if (ioerr /= 0) call exit_POP(sigAbort, &
!                                'Error opening region mask file')

!  if (my_task == master_task) then
!     read(nu, rec=1, iostat=ioerr) REGION_G
!     close(nu)
!  endif
!  call release_unit(nu)

!  call broadcast_scalar(ioerr, master_task)
!  if (ioerr /= 0) call exit_POP(sigAbort, &
!                                'Error reading region mask file')

!  call scatter_global(REGION_MASK, REGION_G, master_task,distrb_clinic, &
!                      field_loc_center, field_type_scalar)
!  if (my_task == master_task) deallocate(REGION_G)

!  num_regions = global_maxval(abs(REGION_MASK), &
!                              distrb_clinic, field_loc_center)
!
!---------------------------------------------------------------------
!
!       open and read file which contains region names and marginal-
!         seas balancing information
!
!---------------------------------------------------------------------


!  if(info_filename == 'unknown_region_info') then
!    call exit_POP (sigAbort,'ERROR: unknown region_info_filename')
!  endif
!
!  region_info(:)%name   = 'unknown_region_name'
!  region_info(:)%area   = c0
!  region_info(:)%volume = c0
!
!  call get_unit(nu)
 
!  if (my_task == master_task) then
 
!    open(nu,file=info_filename,form='formatted', &
!         status='unknown', iostat=ioerr)
!
!    do n = 1,num_regions
!      read (nu,*) region_info(n)%number     , &
!                  region_info(n)%name       , &
!                  region_info(n)%ms_bal%lat , &
!                  region_info(n)%ms_bal%lon , &
!                  region_info(n)%ms_bal%area
!    enddo
!
!    close(nu)
 
!  endif 
!
!  call release_unit(nu)
!
!  call broadcast_scalar(ioerr, master_task)
!  if (ioerr /= 0) call exit_POP(sigAbort, &
!                                'Error reading region name file')

!  do n = 1,num_regions
!    call broadcast_scalar(region_info(n)%name       ,master_task)
!    call broadcast_scalar(region_info(n)%number     ,master_task)
!    call broadcast_scalar(region_info(n)%ms_bal%lat ,master_task)
!    call broadcast_scalar(region_info(n)%ms_bal%lon ,master_task)
!    call broadcast_scalar(region_info(n)%ms_bal%area,master_task)
!  enddo
 
 
!---------------------------------------------------------------------
!
!       determine if region is a marginal sea
!
!---------------------------------------------------------------------
!  num_ms = 0
!  do n = 1,num_regions
!    if (region_info(n)%number < 0 )  then
!        num_ms = num_ms + 1
!        region_info(n)%marginal_sea = .true.
!    else
!        region_info(n)%marginal_sea = .false.
!    endif
!  enddo
!
!  if (num_ms > max_ms) then
!    if (my_task == master_task) then
!      write(stdout,*)'area_masks: maximum number of marginal seas exceeded'   
!    endif
!    call exit_POP(sigAbort, &
!                 'ERROR: must increase max_ms in module grid.F')
!  endif


!-----------------------------------------------------------------------
!
!  a negative value of REGION_MASK designates a marginal sea.  if
!  a region is a marginal sea, calculate the area and volume
!
!-----------------------------------------------------------------------

!  area_t_marg = c0
!  volume_t_marg = c0
!  volume_t_marg_k = c0

!  do region = 1, num_regions

!     WORK = merge(TAREA, c0, 1 <= KMT .and. REGION_MASK == -region)
!     sea_area = global_sum(WORK, distrb_clinic, field_loc_center)
!     area_t_marg = area_t_marg + sea_area

!     sea_vol = c0
!     do k = 1, km
!        WORK = merge(TAREA*dz(k), c0, &
!                     k <= KMT .and. REGION_MASK == -region)
!        tmp_vol = global_sum(WORK, distrb_clinic, field_loc_center)
!        sea_vol = sea_vol + tmp_vol
!        volume_t_marg_k(k) = volume_t_marg_k(k) + tmp_vol
!     enddo
!     volume_t_marg = volume_t_marg + sea_vol

!     if(sea_area /= c0) then
!       if (.not. region_info(region)%marginal_sea) then
!         call exit_POP (sigAbort,'ERROR: marginal-sea mismatch')
!       endif
!       region_info(region)%area   = sea_area
!       region_info(region)%volume = sea_vol
!     endif

!     if (my_task == master_task .and. sea_area /= c0) then
!        write(stdout,"('Region #',i2,' is a marginal sea')") region
!        write(stdout, &
!            "('  area (km^2) = ',e12.5, '  volume (km^3) = ',e12.5)") &
!            sea_area*1.0e-10_POP_r8, sea_vol*1.0e-15_POP_r8
!     endif

!  enddo


!---------------------------------------------------------------------
!
!       document regional information
!
!---------------------------------------------------------------------
 
!  if (my_task == master_task) then
!    write(stdout,blank_fmt)
!    write(stdout,1003)
!
!    do n = 1,num_regions
!     if (region_info(n)%marginal_sea) then
!       write(stdout,1004) region_info(n)%number                   , &
!                          trim(region_info(n)%name)               , &
!                          region_info(n)%area  *1.0e-10_POP_r8        , &
!                          region_info(n)%volume*1.0e-15_POP_r8      
!     else
!       write(stdout,1004) region_info(n)%number, trim(region_info(n)%name)
!     endif
!    enddo
!
!    write(stdout,blank_fmt)
!    write(stdout,1005)
 
!    do n = 1,num_regions
!     if (region_info(n)%marginal_sea) then
!       write(stdout,1006) region_info(n)%number        , &
!                          trim(region_info(n)%name)    , &
!                          region_info(n)%ms_bal%lat    , &
!                          region_info(n)%ms_bal%lon    , &
!                          region_info(n)%ms_bal%area
!     endif
!    enddo
!
!  endif 
!
!  if (my_task == master_task) then
!    write(stdout,blank_fmt)
!  endif

!  do n = 1,num_regions
!    if (region_info(n)%marginal_sea) then
!        region_info(n)%area  =region_info(n)%area  *1.0e-4_POP_r8
!        region_info(n)%volume=region_info(n)%volume*1.0e-6_POP_r8
!    endif
!  enddo


1003  format (  30x, '+', 23('-'),'+'                                     , &
              /,30x, '|',  2x,  'Marginal Seas Only   |'                  , &
              /,2x,'Region', 8x, 'Region ',   7x, '|',  2x,'Area'         , &
              8x,'Volume   |'                                             , &
              /,2x,'Number', 8x,  'Name',10x, '| (km^2)', 7x, '(km^3)'    , &
              3x, '|'                                                     , &
              /, 2x,6('-'), 1x,19('-'),2x,11('-'),2x                      , &
              11('-') )
1004  format (2x, i4, a22, 2(1pe13.5) )
1005  format (/,3x, ' Marginal Sea (E+P+M+R) Balancing Information'       , &
              /,2x,'Region', 8x, 'Region ',  8x, 'Bal.'                   , &
              4x,'Bal.',  4x, 'Search'                                    , &
              /,2x,'Number', 8x,  'Name',11x, 'Lat'  , 5x, 'Lon'          , &
              5x, 'Area'                                                  , &
              /,47x, '(cm^2)'                                             , &
              /, 2x,6('-'), 1x,19('-'),2x, 6('-'),2x                      , &
               6('-'), 2x, 11('-')  )                      
1006  format (1x, i4, a20, 3x,2(0pf8.2),   1pe13.5  )
!-----------------------------------------------------------------------
!EOC

 end subroutine area_masks

!***********************************************************************

 subroutine cf_area_avg

!-----------------------------------------------------------------------
!
!  calculate the coefficients of the 4-point stencils for
!  area-weighted averaging at T and U points
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  calculate {central,s,w,sw} coefficients for area-averaging
!  to T points.
!
!-----------------------------------------------------------------------


!  AT0  = p25
!  ATS  = p25
!  ATW  = p25
!  ATSW = p25

!-----------------------------------------------------------------------
!
!  calculate {central,n,e,ne} coefficients for area-averaging
!   to U points.
!
!-----------------------------------------------------------------------

!  AU0  = TAREA
!  AUN  = eoshift(TAREA,dim=2,shift=+1)
!  AUE  = eoshift(TAREA,dim=1,shift=+1)
!  AUNE = eoshift(AUE  ,dim=2,shift=+1)

!  AU0  = AU0 *p25*UAREA_R
!  AUN  = AUN *p25*UAREA_R
!  AUE  = AUE *p25*UAREA_R
!  AUNE = AUNE*p25*UAREA_R

!-----------------------------------------------------------------------

 end subroutine cf_area_avg

!***********************************************************************
!BOP
! !IROUTINE: calc_tpoints
! !INTERFACE:

!subroutine calc_tpoints(errorCode)

! !DESCRIPTION:
!  Calculates lat/lon coordinates of T points from U points
!  using a simple average of four neighbors in Cartesian 3d space.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

!  integer (POP_i4), intent(out) :: &
!     errorCode                ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!  integer (POP_i4) :: i,j,n

!  real (POP_r8) ::                   &
!     xc,yc,zc,xs,ys,zs,xw,yw,zw, &! Cartesian coordinates for
!     xsw,ysw,zsw,tx,ty,tz,da      !    nbr points

!  type (block) :: &
!     this_block    ! block info for this block

!-----------------------------------------------------------------------
!
!  TLAT, TLON are southwest 4-point averages of ULAT,ULON
!  for general grids, must drop into 3-d Cartesian space to prevent
!  problems near the pole
!
!-----------------------------------------------------------------------

!  !$OMP PARALLEL DO PRIVATE(n,i,j,this_block,xsw,ysw,zsw,xw,yw,zw,xs,ys,zs,xc,yc,zc, &
!  !$OMP                     tx,ty,tz,da)

!  do n=1,nblocks_clinic
!     this_block = get_block(blocks_clinic(n),n)

!     do j=2,ny_block
!     do i=2,nx_block

!        !***
!        !*** set up averaging weights
!        !***

!        !wt0  = AT0 (i,j,n)
!        !wts  = ATS (i,j,n)
!        !wtw  = ATW (i,j,n)
!        !wtsw = ATSW(i,j,n)

!        !***
!        !*** convert neighbor U-cell coordinates to 3-d Cartesian coordinates 
!        !*** to prevent problems with averaging near the pole
!        !***

!        zsw = cos(ULAT(i-1,j-1,n))
!        xsw = cos(ULON(i-1,j-1,n))*zsw
!        ysw = sin(ULON(i-1,j-1,n))*zsw
!        zsw = sin(ULAT(i-1,j-1,n))

!        zs  = cos(ULAT(i  ,j-1,n))
!        xs  = cos(ULON(i  ,j-1,n))*zs
!        ys  = sin(ULON(i  ,j-1,n))*zs
!        zs  = sin(ULAT(i  ,j-1,n))

!        zw  = cos(ULAT(i-1,j  ,n))
!        xw  = cos(ULON(i-1,j  ,n))*zw
!        yw  = sin(ULON(i-1,j  ,n))*zw
!        zw  = sin(ULAT(i-1,j  ,n))

!        zc  = cos(ULAT(i  ,j  ,n))
!        xc  = cos(ULON(i  ,j  ,n))*zc
!        yc  = sin(ULON(i  ,j  ,n))*zc
!        zc  = sin(ULAT(i  ,j  ,n))

!        !***
!        !*** straight 4-point average to T-cell Cartesian coords
!        !***

!        tx = p25*(xc + xs + xw + xsw)
!        ty = p25*(yc + ys + yw + ysw)
!        tz = p25*(zc + zs + zw + zsw)

!        !***
!        !*** convert to lat/lon in radians
!        !***

!        da = sqrt(tx**2 + ty**2 + tz**2)

!        TLAT(i,j,n) = asin(tz/da)

!        if (tx /= c0 .or. ty /= c0) then
!           TLON(i,j,n) = atan2(ty,tx)
!        else
!           TLON(i,j,n) = c0
!        endif
!          
!     end do
!     end do

      !***
      !*** for bottom row of domain where sw 4pt average is not valid,
      !*** extrapolate from interior
      !*** NOTE: THIS ASSUMES A CLOSED SOUTH BOUNDARY - WILL NOT
      !***       WORK CORRECTLY FOR CYCLIC OPTION
      !***

!     if (this_block%j_glob(this_block%jb) == 1) then
!        do i=this_block%ib,this_block%ie
!           TLON(i,this_block%jb,n) = TLON(i,this_block%jb+1,n)
!           TLAT(i,this_block%jb,n) = c2*TLAT(i,this_block%jb+1,n) - &
!                                        TLAT(i,this_block%jb+2,n)
!        end do
!     endif

!     where (TLON(:,:,n) > pi2) TLON(:,:,n) = TLON(:,:,n) - pi2
!     where (TLON(:,:,n) < c0 ) TLON(:,:,n) = TLON(:,:,n) + pi2

!  end do
!  !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  Update boundaries
!
!-----------------------------------------------------------------------

!  call POP_HaloUpdate(TLAT, POP_haloClinic, POP_gridHorzLocCenter,  & 
!                      POP_fieldKindScalar, errorCode,               &
!                      fillValue = 0.0_POP_r8)

!  if (errorCode /= 0) then
!     call POP_ErrorSet(errorCode, &
!        'calc_tpoints: error updating tlat halo')
!     return
!  endif

!  call POP_HaloUpdate(TLON, POP_haloClinic, POP_gridHorzLocCenter,  & 
!                      POP_fieldKindScalar, errorCode,               &
!                      fillValue = 0.0_POP_r8)

!  if (errorCode /= 0) then
!     call POP_ErrorSet(errorCode, &
!        'calc_tpoints: error updating tlon halo')
!     return
!  endif

!-----------------------------------------------------------------------

!end subroutine calc_tpoints

!***********************************************************************
!BOP
! !IROUTINE: ugrid_to_tgrid
! !INTERFACE:

!subroutine ugrid_to_tgrid(ARRAY_TGRID, ARRAY_UGRID, iblock)

! !DESCRIPTION:
!  Interpolates values at U points on a B grid to T points.
!  Note that ghost cells are not updated.
!  Also note that the input array is assumed to be in the baroclinic
!  distribution (where the stencil weights are defined).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

!  integer (POP_i4), intent(in) :: &
!    iblock                  ! index for block in baroclinic distrb

!  real (POP_r8), dimension(nx_block,ny_block), intent(in) :: &
!    ARRAY_UGRID             ! field on U points

! !OUTPUT PARAMETERS:

!  real (POP_r8), dimension(nx_block,ny_block), intent(out) :: &
!    ARRAY_TGRID             ! field on T points

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!  integer (POP_i4) :: &
!    i,j                 ! dummy indices

!-----------------------------------------------------------------------
!
!  southwest 4pt average
!
!-----------------------------------------------------------------------

!  do j=2,ny_block
!  do i=2,nx_block

!    ARRAY_TGRID(i,j) = AT0 (i,j,iblock)*ARRAY_UGRID(i  ,j  ) + &
!                       ATS (i,j,iblock)*ARRAY_UGRID(i  ,j-1) + &
!                       ATW (i,j,iblock)*ARRAY_UGRID(i-1,j  ) + &
!                       ATSW(i,j,iblock)*ARRAY_UGRID(i-1,j-1)

!  end do
!  end do

!  ARRAY_TGRID(:,1) = c0
!  ARRAY_TGRID(1,:) = c0

!-----------------------------------------------------------------------
!EOC

!end subroutine ugrid_to_tgrid

!***********************************************************************
!BOP
! !IROUTINE: tgrid_to_ugrid
! !INTERFACE:

!subroutine tgrid_to_ugrid(ARRAY_UGRID, ARRAY_TGRID, iblock)

! !DESCRIPTION:
!  Interpolates values at T points on a B grid to U points.
!  Note that ghost cells are not updated.
!  Also note that the input array is assumed to be in the baroclinic
!  distribution (where the stencil weights are defined).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

!  integer (POP_i4), intent(in) :: &
!    iblock                  ! index for block in baroclinic distrb

!  real (POP_r8), dimension(nx_block,ny_block), intent(in) :: &
!    ARRAY_TGRID    ! field on T points

! !OUTPUT PARAMETERS:

!  real (POP_r8), dimension(nx_block,ny_block), intent(out) :: &
!    ARRAY_UGRID    ! field on U points

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!  integer (POP_i4) :: &
!    i,j                 ! dummy indices

!-----------------------------------------------------------------------
!
!  northeast 4pt average
!
!-----------------------------------------------------------------------

!  do j=1,ny_block-1
!  do i=1,nx_block-1

!    ARRAY_UGRID(i,j) = AU0 (i,j,iblock)*ARRAY_TGRID(i  ,j  ) + &
!                       AUN (i,j,iblock)*ARRAY_TGRID(i  ,j+1) + &
!                       AUE (i,j,iblock)*ARRAY_TGRID(i+1,j  ) + &
!                       AUNE(i,j,iblock)*ARRAY_TGRID(i+1,j+1)

!  end do
!  end do

!  ARRAY_UGRID(:,ny_block) = c0
!  ARRAY_UGRID(nx_block,:) = c0

!-----------------------------------------------------------------------
!EOC

!end subroutine tgrid_to_ugrid

!***********************************************************************

 end module grid

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

