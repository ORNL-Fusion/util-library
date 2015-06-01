!-----------------------------------------------------------------------------
!
! Tests fieldline following for 2d gfile eq.
! -- JDL
!
Program test_fieldline_follow
 
! Description: 
!
! History:
! Version   Date        Comment
! -------   ----        -------
! 1.0     05/24/2013   Original code. JL 
! 
! Author(s): J. Lore - current

Use kind_mod                     ! Import rknd, iknd specifications
Use gfile_var_pass
Use g3d_module
Use rmp_module
Use fieldline_follow_mod
Use bspline
Implicit None

! Local variables (scalar)
Real(rknd),Allocatable :: Btest(:,:),Rtest(:),Ztest(:),Phitest(:)
Real  :: tarray(2),tres
Real  :: tarray0(2),tres0
! Local variables (array)
Character(Len=100) :: gfilename

Integer(iknd),parameter :: nsteps  = 10000
Integer :: i,ntest
Integer,Allocatable :: ierr_b(:),ilg(:)
Real(rknd) :: dphi
Real(rknd),Allocatable,Dimension(:,:) :: rout,zout,phiout

!- End of header -------------------------------------------------------------


gfilename = '../bfield_files/g135183.00433'

Call readg_g3d(gfilename)



Write(*,*) 'test fieldline following'
Ntest = 2
Allocate(Rtest(Ntest),Ztest(Ntest),Phitest(Ntest))
Allocate(rout(Ntest,nsteps+1))
Allocate(zout(Ntest,nsteps+1))
Allocate(phiout(Ntest,nsteps+1))
Allocate(ierr_b(Ntest),ilg(Ntest))
Rtest = (/1.1d0,1.1d0/)
Ztest = (/0.d0,0.d0/)
phitest = (/0.0d0,0.0d0/)
dphi = -0.5d0*3.1415d0/180.d0
if (.true.) Then
call etime(tarray0,tres0)
Call follow_fieldlines_rzphi(Rtest,Ztest,phitest,Ntest,dphi,nsteps,rout,zout,phiout,ierr_b,ilg,1)
call etime(tarray,tres)
print *,'Time: ',tres - tres0
!!$open(99,file='outline.out')
!!$write(99,*) rout
!!$write(99,*) zout
!!$write(99,*) phiout
!!$close(99)
write(*,*) 'ierr_b',ierr_b
Write(*,*) 'i_last_good',ilg,nsteps+1
endif



ierr_b = 0
ilg = 0
Write(*,*)
Write(*,*) 'test fieldline following  -- 2 --'
rmp_method = 0
call etime(tarray0,tres0)
Call ffr(Rtest,Ztest,phitest,Ntest,dphi,nsteps,rout,zout,phiout,ierr_b,ilg)
call etime(tarray,tres)
print *,'Time: ',tres - tres0
!!$open(99,file='outline2.out')
!!$write(99,*) rout
!!$write(99,*) zout
!!$write(99,*) phiout
!!$close(99)
write(*,*) 'ierr_b',ierr_b
Write(*,*) 'i_last_good',ilg,nsteps+1




End program test_fieldline_follow

