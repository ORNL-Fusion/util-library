!-----------------------------------------------------------------------------
!
! Tests bfield_geq_bicub
! -- JDL
!
Program test_bfield_geq_bicub
 
! Description: 
!
! History:
! Version   Date        Comment
! -------   ----        -------
! 1.0     05/24/2013   Original code. JL 
! 
! Author(s): J. Lore - current

Use kind_mod                     ! Import rknd, iknd specifications
Use g3d_module
Use bspline
Implicit None

! Local variables (scalar)
Real(real64),Allocatable :: Btest(:,:),Rtest(:),Ztest(:)
Real  :: tarray(2),tres
! Local variables (array)
Character(Len=100) :: gfilename

Integer(int32),parameter :: npts_line  = 10000
Integer :: i,ierr_b,Ntest

Real(real64),dimension(npts_line+1) :: rout,zout,phiout

!- End of header -------------------------------------------------------------


gfilename = '../bfield_files/g135183.00433'

Call readg_g3d(gfilename)

Write(*,*) 'Test call to bfield_geq_bicub'
Ntest = 2
Allocate(Rtest(Ntest),Ztest(Ntest),Btest(Ntest,3))
Rtest = (/1.2d0,0.4d0/)
Ztest = (/-0.8d0,-1.45d0/)
Call bfield_geq_bicub(Rtest,Ztest,Ntest,Btest,ierr_b)
Write(*,*) 'And the answer: ',Btest(1,:),ierr_b
Write(*,*) 'And the answer: ',Btest(2,:),ierr_b

Deallocate(g_fpol,g_pres,g_ffprim,g_pprime,g_psirz,g_qpsi)
Deallocate(g_r,g_z,g_bdry,g_lim,g_pn)
Deallocate(g_bicub_coeffs)

Write(*,*) 'Done with test_read_gfile'

End program test_bfield_geq_bicub

