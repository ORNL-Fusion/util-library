!-----------------------------------------------------------------------------
!
! Tests reading of a gfile and outputs read info to screen.
! -- JDL
!
Program test_read_gfile
 
! Description: 
!
! History:
! Version   Date        Comment
! -------   ----        -------
! 1.0     05/24/2030   Original code. JL 
! 
! Author(s): J. Lore - current

Use kind_mod                     ! Import rknd, iknd specifications
Use gfile_var_pass
Use g3d_module
Use bspline
!Use g3df_functions_mod, Only : &
!  bfield_geq_bicub
Implicit None

! Local variables (scalar)
Real(real64) :: Btest(2,3)
Real  :: tarray(2),tres
Real(real64) :: tmp,pn_test(100),fp_test(100)
Type(gdata) :: g
! Local variables (array)
Character(Len=100) :: gfilename

Integer(int32),parameter :: npts_line  = 10000
Integer :: i

Real(real64),dimension(npts_line+1) :: rout,zout,phiout

!- End of header -------------------------------------------------------------



gfilename = './g160884.03014_251'

Call readg_g3d(gfilename,g)

call etime(tarray,tres)
print *,'Time: ',tres

!Call display_gfile

Call close_gfile
Write(*,*) 'Done with test_read_gfile'


End program test_read_gfile

