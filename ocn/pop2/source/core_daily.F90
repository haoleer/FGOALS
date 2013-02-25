!     =================
      SUBROUTINE CORE_DAILY(TNUM)
!     =================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use forc_mod
use dyn_mod, only: u,v,buffer
use tracer_mod, only: at
use msg_mod
      IMPLICIT NONE
#include <netcdf.inc>

    integer :: TNUM
    return
end subroutine CORE_DAILY
