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
Use NSTX_routines_mod
Use bspline
Implicit None

! Local variables (scalar)
Real(rknd),Allocatable :: Btest(:,:),Rtest(:),Ztest(:),Phitest(:)
Real  :: tarray(2),tres
Real  :: tarray0(2),tres0
! Local variables (array)
Character(Len=100) :: gfilename

Integer(iknd),parameter :: nsteps  = 10000
Integer :: i,ntest,ntorpts
Integer,Allocatable :: ierr_b(:),ilg(:)
Real(rknd) :: dphi, taper(6)
Real(rknd),Allocatable,Dimension(:,:) :: rout,zout,phiout

!- End of header -------------------------------------------------------------
Call Etime(tarray,tres0)

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


Write(*,*) 'building coils'
taper = (/1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0/)
taper = taper*1000.d0
ntorpts = 5
Allocate(rmp_coil(6*(2*ntorpts+1),3))
Allocate(rmp_current(6*(2*ntorpts+1)))
Call build_nstx_rwmcoils_jl(taper,ntorpts,rmp_coil,rmp_current,rmp_ncoil_pts)


ierr_b = 0
Write(*,*)
Write(*,*) 'test fieldline following  -- --'
rmp_method = 1
call etime(tarray0,tres0)
Call ffr(Rtest,Ztest,phitest,Ntest,dphi,nsteps,rout,zout,phiout,ierr_b,ilg)
call etime(tarray,tres)
print *,'Time: ',tres - tres0
!open(99,file='outline2.out')
!write(99,*) rout
!write(99,*) zout
!write(99,*) phiout
!close(99)
write(*,*) 'ierr_b',ierr_b
Write(*,*) 'i_last_good',ilg,nsteps+1



End program test_fieldline_follow

