!-----------------------------------------------------------------------------
!
! Tests fieldline following for 2d gfile eq, with intersection
! -- JDL
!
Program test_fieldline_follow_with_int
 
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
Use fieldline_follow_mod
Use math_geo_module
Use bspline
Implicit None

! Local variables (scalar)
Real(rknd),Allocatable :: Btest(:,:),Rtest(:),Ztest(:),Phitest(:)
Real  :: tarray(2),tres
! Local variables (array)
Character(Len=100) :: gfilename

Integer(iknd),parameter :: nsteps  = 10000
Integer(iknd) :: i,ntest,iline,jpt,npts,nline,inside
Integer,Allocatable :: ierr_b(:),ilg(:),ithit(:),hit_jpt(:)
Real(rknd) :: dphi
Real(rknd),Allocatable,Dimension(:,:) :: rout,zout,phiout

!- End of header -------------------------------------------------------------


!gfilename = '../bfield_files/g135183.00433'
gfilename = '../bfield_files/g140508.00403'

Call readg_g3d(gfilename)



Write(*,*) 'test fieldline following'
Ntest = 2
Allocate(Rtest(Ntest),Ztest(Ntest),Phitest(Ntest))
Allocate(rout(Ntest,nsteps+1))
Allocate(zout(Ntest,nsteps+1))
Allocate(phiout(Ntest,nsteps+1))
Allocate(ierr_b(Ntest),ilg(Ntest))
Rtest = (/1.08d0,1.2d0/)
Ztest = (/-1.3d0,0.d0/)
phitest = (/0.0d0,0.0d0/)
dphi = -0.5d0*3.1415d0/180.d0
Call follow_fieldlines_rzphi(Rtest,Ztest,phitest,Ntest,dphi,nsteps,rout,zout,phiout,ierr_b,ilg,1)
call etime(tarray,tres)
print *,'Time: ',tres
open(99,file='outline.out')
write(99,*) rout
write(99,*) zout
write(99,*) phiout
close(99)
write(*,*) 'ierr_b',ierr_b
Write(*,*) 'i_last_good',ilg,nsteps+1

nline = Ntest
npts = nsteps + 1
Allocate(ithit(nline))
Allocate(hit_jpt(nline))
hit_jpt = 0
ithit = 0
Do iline = 1,nline
  Do jpt = 1,npts
    inside = inside_poly(rout(iline,jpt),zout(iline,jpt),g_lim(1,:),g_lim(2,:),g_limitr)
    if (inside .ne. 1) Then 
      ithit(iline) = 1
      hit_jpt(iline) = jpt
      Exit
    Endif
  Enddo
Enddo
write(*,*) 'It hit?',ithit
write(*,*) 'hit jpt',hit_jpt

End program test_fieldline_follow_with_int

