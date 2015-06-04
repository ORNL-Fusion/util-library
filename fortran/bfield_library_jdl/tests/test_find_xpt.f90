!-----------------------------------------------------------------------------
!
! Tests xpoint finding routine
! -- JDL
!
Program test_read_gfile
  ! 6/4/15
  ! Author(s): J. Lore - current

Use kind_mod                     ! Import rknd, iknd specifications
Use g3d_module
Implicit None


Character(Len=100) :: gfilename
Real(rknd) :: rx,zx,rx2,zx2
!- End of header -------------------------------------------------------------


gfilename = '/home/jjl/gfiles/DIII-D/g160884.03014_251'

Call readg_g3d(gfilename)
Call find_xpt_jdl(.true.,.true.,1.d-8,.false.,rx,zx,rx2,zx2)


End program test_read_gfile

