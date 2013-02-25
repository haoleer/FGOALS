!  CVS: $Id: param_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module param_mod
#include <def-undef.h>
use precision_mod
!-------------------------------------------------------------------------------
!
! Author: Yongqiang YU  ( 12 Nov, 2002)
!
!-------------------------------------------------------------------------------
!YU  01/02/2013
      integer, parameter:: LICOM_blockSizeX = BLCKX ! size of block in first  horizontal dimension
      integer, parameter:: LICOM_blockSizeY = BLCKY ! size of block in second horizontal dimension
      integer, parameter:: max_blocks_clinic = MXBLCKS, max_blocks_tropic = MXBLCKS    !   in each distribution
      integer , parameter, public :: nghost = 2       ! number of ghost cells around each block
      integer , parameter, public :: nx_block = licom_BlockSizeX + 2*nghost,   &!  x,y dir including ghost
                                     ny_block = licom_BlockSizeY + 2*nghost     !  cells
!YU
      integer,parameter:: jmt_global=JMT_GLOBAL  ! Number of the End Grid for Tracer in Latitude.
      integer,parameter:: jmm_global=jmt_global-1
      integer,parameter:: jstart    = 3     ! Satrting grid for Tracer in Latitude.
      integer,parameter:: imt_global=360   ! Number of Grid Points in Longitude
      integer,parameter:: km=30     ! Number of Grid Points in Vertical Direction
!
      integer,parameter:: nx_proc=NX_PROC ! Number of MPI tasks in zonal direction
      integer,parameter:: ny_proc=NY_PROC ! Number of MPI tasks in meridional direction
      integer,parameter:: n_proc=nx_proc*ny_proc ! Total number of Processors for MPI
      integer,parameter:: num_overlap=2 ! Number of overlapping grids for subdomain.
      integer,parameter:: jst_global=jstart ! Number of the Strating Grid for Tracer in Latitude.

!Nummber of grids in the each subdomain
      integer,parameter:: jst=1     ! Number of the Strating Grid for Tracer in Latitude.
      integer,parameter:: jsm=jst+1 ! Number of the Strating Grid for Momentum in Latitude.
!ZHW should Add 1
      integer,parameter:: jet= ny_block
      integer,parameter:: jem=jet-1 ! Number of the End Grid for Momentum in Latitude.
      integer,parameter:: jmt= ny_block
      integer,parameter:: imt= nx_block
      integer :: j_loop             ! Loop index of J cycle for the each subdomain.
!
      integer,parameter:: imm_global=imt_global-1,imm=imt-1,jmm=jmt-1,kmp1=km+1,kmm1=km-1
      integer,parameter:: ntra=2    ! Number of Tracers
      integer:: i,j,k,m,n,ierr, mytid
      integer::jj_start,jj_end ! No Overlaped j ,North to South
      integer:: my_task, master_task

      real(r8) :: dlam   ! Zonal grid distance in degree
      real(r8) :: am_tro ! Horizontal viscosity in the tropics.
      real(r8) :: am_ext ! Horizontal viscosity in the extra-tropics.

!lhl090729
      integer,parameter:: s_imt=192,s_jmt=94
!lhl090729
!YU  01/02/2013



end module param_mod
