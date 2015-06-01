!-----------------------------------------------------------------------------
!
! Tests fieldline following for 2d gfile eq + screened b-spline field
! -- JDL
!
Program test_fieldline_follow_screened
 
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
Use screening_module
Use fieldline_follow_mod
Use NSTX_routines_mod
Use bspline
Implicit None

! Local variables (scalar)
Real(rknd),Allocatable :: Btest(:,:),Rtest(:),Ztest(:),Phitest(:)
Real  :: tarray(2),tres
Real  :: tarray0(2),tres0
! Local variables (array)
Character(Len=100) :: gfilename

Integer(iknd),parameter :: nsteps  = 100
Integer :: i,ntest,ntorpts
Integer(iknd) :: ierr
Integer,Allocatable :: ierr_b(:),ilg(:)
Real(rknd) :: dphi, taper(6)
Real(rknd),Allocatable,Dimension(:,:) :: rout,zout,phiout,Bout
Character(Len=200) :: Afile_path
!- End of header -------------------------------------------------------------


gfilename = '../bfield_files/g135183.00433'
Call readg_g3d(gfilename)


Afile_path = '/home/jjl/Pavel/screened-8-19/'
Call setup_screening_vars(Afile_path)
!call etime(tarray0,tres0)
Call read_Afiles(Afile_path)
Call prepare_Afile_splinefits !(Afile_path)
call etime(tarray,tres)


Write(*,*) 'test fieldline following'
Ntest = 2
Allocate(Rtest(Ntest),Ztest(Ntest),Phitest(Ntest))
Allocate(rout(Ntest,nsteps+1))
Allocate(zout(Ntest,nsteps+1))
Allocate(phiout(Ntest,nsteps+1))
Allocate(ierr_b(Ntest),ilg(Ntest))
Allocate(Bout(Ntest,3))
Rtest = (/1.1d0,1.1d0/)
Ztest = (/0.d0,0.d0/)
phitest = (/0.0d0,0.0d0/)
dphi = -0.5d0*3.1415d0/180.d0

Write(*,*) 'test beval'
Call bfield_bspline(Rtest,phitest,Ztest,Ntest,Bout,ierr)
Write(*,*) Bout,'<--------'


ierr_b = 0
Write(*,*)
Write(*,*) 'test fieldline following  -- --'
rmp_method = 2
call etime(tarray0,tres0)
Call ffr(Rtest,Ztest,phitest,Ntest,dphi,nsteps,rout,zout,phiout,ierr_b,ilg)
call etime(tarray,tres)
print *,'Time: ',tres - tres0
open(99,file='outline3.out')
write(99,*) rout
write(99,*) zout
write(99,*) phiout
close(99)
write(*,*) 'ierr_b',ierr_b
Write(*,*) 'i_last_good',ilg,nsteps+1




End program test_fieldline_follow_screened

