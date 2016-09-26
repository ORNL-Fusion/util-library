!-----------------------------------------------------------------------------
!
!   Routines/modules related to rmp COIL calculations
!   
!-----------------------------------------------------------------------------
Module rmpcoil_module
  Use kind_mod, Only: real64, int32
  Use DIIID_routines_mod, Only : build_d3d_ccoils_jl,build_d3d_icoils_jl
  Use NSTX_routines_mod, Only : build_NSTX_rwmcoils_jl
  Use biotsavart_module, Only : bfield_bs_cyl, bfield_bs_jdl
  Use coil_typedef, Only : coil_type
  Implicit None
  Public
End Module rmpcoil_module
