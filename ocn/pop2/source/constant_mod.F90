!=======================================================================
! CVS $Id: constant_mod.F90,v 1.2 2001/12/11 22:52:37 kauff Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/ocn/docn5/constant_mod.F90,v $
! CVS $Name: ccsm2_0_beta58 $
!=======================================================================

module constant_mod

  use precision_mod
  use shr_kind_mod
  use shr_const_mod
  use netcdf

  implicit none

  !----- constants -----
  real(SHR_KIND_R8), parameter ::  latvap   = SHR_CONST_LATVAP ! latent heat of evap   ~ J/kg
  real(SHR_KIND_R8), parameter ::  Tzro     = SHR_CONST_TKFRZ  ! 0 degrees C                       ~ kelvin
  real(SHR_KIND_R8), parameter ::  Tfrz     = Tzro   - 1.8     ! temp of saltwater freezing ~ kelvin
  real(SHR_KIND_R8), parameter ::  pi       = SHR_CONST_PI     ! a famous math constant
  real(SHR_KIND_R8), parameter ::  omega    = SHR_CONST_OMEGA  ! earth's rotation  ~ rad/sec
  real(SHR_KIND_R8), parameter ::  g        = SHR_CONST_G      ! gravity ~ m/s^2
  real(SHR_KIND_R8), parameter ::  DEGtoRAD = PI/180.0         ! PI/180


   real (r8), parameter, public :: &
      c0     =    0.0_r8   ,&
      c1     =    1.0_r8   ,&
      c2     =    2.0_r8   ,&
      c3     =    3.0_r8   ,&
      c4     =    4.0_r8   ,&
      c5     =    5.0_r8   ,&
      c8     =    8.0_r8   ,&
      c10    =   10.0_r8   ,&
      c16    =   16.0_r8   ,&
      c1000  = 1000.0_r8   ,&
      c10000 =10000.0_r8   ,&
      c1p5   =    1.5_r8   ,&
      p33    = c1/c3       ,&
      p5     = 0.500_r8    ,&
      p25    = 0.250_r8    ,&
      p125   = 0.125_r8    ,&
      p001   = 0.001_r8    ,&
      eps    = 1.0e-10_r8  ,&
      eps2   = 1.0e-20_r8  ,&
      bignum = 1.0e+30_r8

   real (r4), parameter, public ::         &
      undefined_nf_r4  = NF90_FILL_FLOAT,  &
      undefined        = -12345._r4

   real (r8), parameter, public ::         &
      undefined_nf_r8  = NF90_FILL_DOUBLE

   real (r8), public ::  &
      undefined_nf = NF90_FILL_DOUBLE

   integer (int_kind), parameter, public ::   &
      undefined_nf_int = NF90_FILL_INT

   !*** location of fields for staggered grids

   integer (int_kind), parameter, public ::   &
      field_loc_unknown  =  0, &
      field_loc_noupdate = -1, &
      field_loc_center   =  1, &
      field_loc_NEcorner =  2, &
      field_loc_Nface    =  3, &
      field_loc_Eface    =  4

   !*** field type attribute - necessary for handling
   !*** changes of direction across tripole boundary

   integer (int_kind), parameter, public ::   &
      field_type_unknown  =  0, &
      field_type_noupdate = -1, &
      field_type_scalar   =  1, &
      field_type_vector   =  2, &
      field_type_angle    =  3

end module constant_mod

