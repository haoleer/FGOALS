!  CVS: $Id: output_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module output_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
      real(r4),dimension(imt,jmt,max_blocks_clinic)::z0mon,himon,hdmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::lthfmon,sshfmon,lwvmon,swvmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::sumon,svmon
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wsmon,tsmon,ssmon,usmon,vsmon
      real(r4),dimension(imt,jmt,2,max_blocks_clinic):: icmon
      real(r4),dimension(imt,jmt,NTRA,max_blocks_clinic)::netmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::mldmon
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::akmmon,aktmon,aksmon
      real(r4),dimension(jmt_global,km+1,2,max_blocks_clinic):: psi
      real(r4),dimension(imt,jmt_global):: bsf
      real(r4),dimension(jmt_global,2,NTRA):: mth,mth_adv,mth_dif,mth_adv_iso
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::trendmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::axmon,aymon,azmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dxmon,dymon,dzmon
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::penmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ddymon

#if (defined SMAG_OUT)
      real(r4),dimension(imt,jmt,km)::a3mon
#endif
      real(r4), parameter :: spval =1.0e+35
end module output_mod
