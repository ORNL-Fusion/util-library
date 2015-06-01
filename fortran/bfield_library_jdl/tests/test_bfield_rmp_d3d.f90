!-----------------------------------------------------------------------------
!
! Tests bfield_geq_bicub
! -- JDL
!
Program test_bfield_rmp_d3d
 
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
Use bspline
Use bfield_module
Use diiid_routines_mod
Implicit None

! Local variables (scalar)
Real(rknd),Allocatable :: Btest(:,:),Rtest(:),Ztest(:)
Real  :: tarray(2),tres
! Local variables (array)
Character(Len=100) :: gfilename

Integer(iknd),parameter :: npts_line  = 10000
Integer :: i,ierr_b,Ntest

Real(rknd),dimension(npts_line+1) :: rout,zout,phiout

Real(rknd),allocatable :: rmp_coil(:,:),rmp_current(:)
Real(Rknd) :: taper(12)
Integer(iknd) :: ntorpts,ncoil_pts

Real(rknd), allocatable, Dimension(:) :: Bx,By,Bz,P_x,P_y,P_z,Phitest,Br,Bphi

!- End of header -------------------------------------------------------------


gfilename = '../bfield_files/g160884.04015_915'

Call readg_g3d(gfilename)

Write(*,*) 'Test call to bfield_geq_bicub'
Ntest = 2
Allocate(Rtest(Ntest),Ztest(Ntest),Phitest(Ntest),Btest(Ntest,3))
Rtest = (/1.8d0,1.8d0/)
Ztest = (/-0.5d0,0.5d0/)
Phitest = (/0.0d0,0.1d0/)

Write(*,*) '--------------------------------'
Write(*,*) '         TEST POINTS'
Write(*,*) 'R (m)',Rtest
Write(*,*) 'Z (m)',Ztest
Write(*,*) 'P (r)',Phitest
Write(*,*) '--------------------------------'

Call bfield_geq_bicub(Rtest,Ztest,Ntest,Btest,ierr_b)
Write(*,*) 'And the answer: ',Btest(1,:),ierr_b
Write(*,*) 'And the answer: ',Btest(2,:),ierr_b

Write(*,*) 'building coils'
taper = (/-3740.d0, 3863.d0,-3720.d0, 3855.d0,-3718.d0, 3858.d0, &
           3862.d0,-3791.d0, 3884.d0,-3854.d0, 3923.d0,-3847.d0/)
ntorpts = 6
Allocate(rmp_coil(12*(2*ntorpts+1),3))
Allocate(rmp_current(12*(2*ntorpts+1)))
Call build_d3d_icoils_jl(taper,ntorpts,rmp_coil,rmp_current,ncoil_pts)

!do i =1,156
!  write(*,*) rmp_coil(i,1:3)
!enddo
!stop

Write(*,*) 'Test call to bfield_bs_jdl'
Allocate(P_x(Ntest),P_y(Ntest),P_z(Ntest))
Allocate(Bx(Ntest),By(Ntest),Bz(Ntest))
Allocate(Br(Ntest),Bphi(Ntest))
P_x = Rtest*cos(Phitest)
P_y = Rtest*sin(Phitest)
P_z = Ztest
Call bfield_bs_jdl(P_x,P_y,P_z,Ntest,rmp_coil,rmp_current,ncoil_pts,Bx,By,Bz)

Write(*,*) 'Bx',Bx
Write(*,*) 'By',By
Write(*,*) 'Bz',Bz


Write(*,*) 'Test call to bfield_bs_cyl'
write(*,*) '  '
Call bfield_bs_cyl(Rtest,Phitest,Ztest,Ntest,rmp_coil,rmp_current,ncoil_pts,Br,Bphi,Bz)

Write(*,*) 'Br',Br
Write(*,*) 'Bphi',Bphi
Write(*,*) 'Bz',Bz

End program test_bfield_rmp_d3d

