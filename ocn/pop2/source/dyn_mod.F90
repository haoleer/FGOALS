!  CVS: $Id: dyn_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module dyn_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt,max_blocks_clinic)::ub,vb,ubp,vbp,h0p
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::up,vp
      real(r8),dimension(imt,jmt,kmp1,max_blocks_clinic)::ws
      real(r8),dimension(imt,jmt,max_blocks_clinic)::h0l,h0f,h0bl,h0bf
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::utl,utf,vtl,vtf
      real(r8),allocatable,dimension(:,:) :: buffer
!
      real(r8),allocatable,dimension(:,:,:)::h0
      real(r8),allocatable,dimension(:,:,:,:)::u,v
!
#ifdef COUP
      real(r4),dimension(imt_global,jmt_global)::t_cpl_io,s_cpl_io,u_cpl_io,v_cpl_io,dhdx_io,dhdy_io,q_io
#endif
!
!     ------------------------------------------------------------------
!     Pressure gradient
!     ------------------------------------------------------------------
      real(r8),dimension(:,:,:,:),allocatable::gg,dlu,dlv
      real(r8),dimension(:,:,:),allocatable::dlub,dlvb
end module dyn_mod

