!-----------------------------------------------------------------------------
!
!   Routines/modules related to rmp COIL calculations
!   
!-----------------------------------------------------------------------------
Module rmpcoil_module
Use kind_mod, Only: real64, int32
Use DIIID_routines_mod, Only : build_d3d_ccoils_jl,build_d3d_icoils_jl
Use NSTX_routines_mod, Only : build_NSTX_rwmcoils_jl
Use bfield_module, Only : bfield_bs_cyl, bfield_bs_jdl
Implicit None

Real(real64), Allocatable :: rmp_coil(:,:)  ! (ncoil_pts,3)
Real(real64), Allocatable ::  rmp_coil_current(:)  ! ncoil_pts
Integer(int32) :: rmp_ncoil_pts

End Module rmpcoil_module
