!  CVS: $Id: tracer_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module tracer_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------

implicit none

      real(r8),dimension(imt,jmt,0:km,NTRA,max_blocks_clinic)::atb
      real(r8),dimension(imt,jmt,NTRA,max_blocks_clinic)::net
!mohr
      real(r8),allocatable,dimension(:,:,:,:,:)::at
!
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::pdensity
      real(r8),dimension(imt,jmt,max_blocks_clinic)::amld
!lhl1204
!
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::trend
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ax,ay,az
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dx,dy,dz
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::penetrate
!
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ddy
!
#ifdef ISO
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::aay_iso,ddy_iso
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ax_iso,ay_iso,az_iso
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dx_iso,dy_iso,dz_iso
#endif
!
!     ------------------------------------------------------------------
!     Sea Ice
!     ------------------------------------------------------------------
!
      real(r8),dimension(imt,jmt,max_blocks_clinic),public:: licomqice
!
end module tracer_mod
