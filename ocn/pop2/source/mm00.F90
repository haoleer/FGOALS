!  CVS: $Id: mm00.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ===================
      SUBROUTINE MM00 (KLV)
!     ===================
 
#include <def-undef.h>
use param_mod
use output_mod
      IMPLICIT NONE
 
      INTEGER :: KLV
 
            Z0MON = 0.0
            HIMON = 0.0
            HDMON = 0.0
            ICMON = 0.0
            ICMON = 0.0
            sumon = 0.0   !U windstress
            svmon = 0.0   !V windstress
            lthfmon = 0.0 !latent flux
            sshfmon = 0.0 !sensible flux
            lwvmon = 0.0  !long wave flux
            swvmon = 0.0  !shortwave flux
!
               TSMON = 0.0
               SSMON = 0.0
               USMON = 0.0
               VSMON = 0.0
               WSMON = 0.0
#if (defined SMAG_OUT)
               AM3MON = 0.0
#endif
               trendmon = 0.0
               axmon = 0.0
               aymon = 0.0
               azmon = 0.0
               dxmon = 0.0
               dymon = 0.0
               dzmon = 0.0
!
               ddymon = 0.0
#ifdef ISO
!              axmon_iso = 0.0
!              aymon_iso = 0.0
!              azmon_iso = 0.0
!              dxmon_iso = 0.0
!              dymon_iso = 0.0
!              dzmon_iso = 0.0
!              aaymon_iso = 0.0
!              ddymon_iso = 0.0
#endif
               netmon = 0.0
               penmon = 0.0
               mldmon = 0.0
               akmmon = 0.0
               aktmon = 0.0
               aksmon = 0.0
 
      RETURN
      END SUBROUTINE MM00
 
 
