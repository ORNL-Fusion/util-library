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
Use g3d_module
Implicit None

! Local variables (scalar)
Real  :: tarray(2),tres
Type(g_type) :: g
! Local variables (array)
Character(Len=100) :: gfilename

!- End of header -------------------------------------------------------------



gfilename = './g160884.03014_251'

Call readg_g3d(gfilename,g)

call etime(tarray,tres)
print *,'Time: ',tres

!Call display_gfile

Call close_gfile(g)
Write(*,*) 'Done with test_read_gfile'


End program test_read_gfile

