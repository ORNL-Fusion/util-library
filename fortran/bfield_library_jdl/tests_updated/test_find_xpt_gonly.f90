!-----------------------------------------------------------------------------
!
! Tests xpoint finding routine
! -- JDL
!
Program test_find_xpt
  ! 6/4/15
  ! Author(s): J. Lore - current

Use kind_mod, Only: real64
Use g3d_module, Only: readg_g3d, g_type
Use util_routines, Only: find_xpt_jdl
Use bfield, Only : bfield_type
Implicit None
Character(Len=100) :: gfilename
Type(g_type) :: g
Type(bfield_type) :: bfield
Real(real64) :: rx,zx,rx2,zx2
!- End of header -------------------------------------------------------------


!gfilename = './g160884.03014_251'
gfilename = '/home/jjl/DIII-D/164723/g164723.03059_410'


Call readg_g3d(gfilename,g)

bfield%method = 1
bfield%g = g

!Subroutine find_xpt_jdl(second,refine,tol,quiet,rx,zx,rx2,zx2,phi_eval_deg,dx)
Call find_xpt_jdl(bfield,.true.,.true.,1.d-10,.false.,rx,zx,rx2,zx2)


End program test_find_xpt

