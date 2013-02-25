!  CVS: $Id: local_to_global.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      subroutine local_to_global(hist_output,rest_output)
!     ================
!     To transfer global 1-d data to local processor.
 
#include <def-undef.h>
#ifdef SPMD
use precision_mod
use param_mod
use pconst_mod
use output_mod
use dyn_mod
use tracer_mod
use work_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn

      IMPLICIT NONE
      integer :: itag
      logical :: hist_output,rest_output 

     if (hist_output) then
!        call local_to_global_4d(z0mon,z0mon_io,1,1)
!        call local_to_global_4d(himon,himon_io,1,1)
!        call local_to_global_4d(hdmon,hdmon_io,1,1)
!        call local_to_global_4d(netmon,netmon_io,2,1)
!        call local_to_global_4d(tsmon,tsmon_io,km,1)
!        call local_to_global_4d(ssmon,ssmon_io,km,1)
!        call local_to_global_4d(usmon,usmon_io,km,1)
!        call local_to_global_4d(vsmon,vsmon_io,km,1)
!        call local_to_global_4d(wsmon,wsmon_io,km,1)
!        call local_to_global_4d(icmon,icmon_io,2,1)
#if (defined SMAG_OUT)
!        call local_to_global_4d(am3mon,am3mon_io,km,1)
#endif
!
!        call local_to_global_4d(axmon,axmon_io,km,2)
!        call local_to_global_4d(aymon,aymon_io,km,2)
!        call local_to_global_4d(azmon,azmon_io,km,2)
!        call local_to_global_4d(dxmon,dxmon_io,km,2)
!        call local_to_global_4d(dymon,dymon_io,km,2)
!        call local_to_global_4d(dzmon,dzmon_io,km,2)
!        call local_to_global_4d(ddymon,ddymon_io,km,2)
#ifdef ISO
!        call local_to_global_4d(axmon_iso,axmon_iso_io,km,2)
!        call local_to_global_4d(aymon_iso,aymon_iso_io,km,2)
!        call local_to_global_4d(azmon_iso,azmon_iso_io,km,2)
!        call local_to_global_4d(dxmon_iso,dxmon_iso_io,km,2)
!        call local_to_global_4d(dymon_iso,dymon_iso_io,km,2)
!        call local_to_global_4d(dzmon_iso,dzmon_iso_io,km,2)
!
!        call local_to_global_4d(aaymon_iso,aaymon_iso_io,km,2)
!        call local_to_global_4d(ddymon_iso,ddymon_iso_io,km,2)
#endif
!        call local_to_global_4d(trendmon,trendmon_io,km,2)
!        call local_to_global_4d(penmon,penmon_io,km,1)
!        call local_to_global_4d(mldmon,mldmon_io,1,1)
!        call local_to_global_4d(akmmon,akmmon_io,km,1)
!        call local_to_global_4d(aktmon,aktmon_io,km,1)
!        call local_to_global_4d(aksmon,aksmon_io,km,1)
!
     end if
!
     if (rest_output) then
!        call local_to_global_4d_double(h0,h0_io,1,1)
!        call local_to_global_4d_double(u,u_io,km,1)
!        call local_to_global_4d_double(v,v_io,km,1)
!        call local_to_global_4d_double(at,at_io,km,2)
     end if

!           end do
#endif
      return
      end subroutine local_to_global
! 
