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
Real(rknd) :: Btest(2,3)
Real  :: tarray(2),tres
Real(rknd) :: tmp,pn_test(100),fp_test(100)
! Local variables (array)
Character(Len=100) :: gfilename

Integer(iknd),parameter :: npts_line  = 10000
Integer :: i

Real(rknd),dimension(npts_line+1) :: rout,zout,phiout

!- End of header -------------------------------------------------------------


gfilename = '../bfield_files/g135183.00433'

Call readg_g3d(gfilename)

call etime(tarray,tres)
print *,'Time: ',tres

!do i=1,100
!  pn_test(i) = dble(i-1)/(100.d0-1.d0)
!  fp_test(i) = dbsval(pn_test(i),g_bspl_ord,g_pnknot,g_mw,g_fpol_bscoef)
!enddo
!writE(*,*) pn_test
!write(*,*) '------'
!write(*,*) fp_test

write(*,*) dbsval(0.d0,g_bspl_ord,g_pnknot,g_mw,g_fpol_bscoef)

!Write(*,*) 'Test call to bfield_geq_bicub'
!Btest = bfield_geq_bicub((/1.2d0,1.2d0/),(/-0.5d0,0.5d0/),2)
!Write(*,*) 'And the answer: ',Btest(1,:)
!Write(*,*) 'And the answer: ',Btest(2,:)
!
!Write(*,*) 'test fieldline following'
!Call follow_fieldlines_rzphi((/1.3d0/),(/0.5d0/),(/0.d0/),1,0.01d0,npts_line,rout,zout,phiout)
!call etime(tarray,tres)
!print *,'Time: ',tres

!Write(*,*) 'done fieldline following'
!Write(*,*) rout(npts_line+1)

Deallocate(g_fpol,g_pres,g_ffprim,g_pprime,g_psirz,g_qpsi)
Deallocate(g_r,g_z,g_bdry,g_lim,g_pn)
Deallocate(g_bicub_coeffs)

Write(*,*) 'Done with test_read_gfile'

End program test_read_gfile

