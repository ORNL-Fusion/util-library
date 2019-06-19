!-----------------------------------------------------------------------------
!
! -- JDL
!
Program test_bfield_bgrid
Use kind_mod, Only: int32, real64
Use bgrid_module
Use phys_const, Only: pi
Implicit None

! Local variables (scalar)
Real(real64),Allocatable :: Bftest(:,:),Rtest(:),Ztest(:)
! Local variables (array)
Character(Len=100) :: fname, run_path

Integer :: Ntest,ierr,ifield_type

Integer(int32) :: i
Real(real64), allocatable, Dimension(:) :: Phitest

!- End of header -------------------------------------------------------------

Write(*,*) 'Test call to bfield_bgrid'
Ntest = 1
Allocate(Rtest(Ntest),Ztest(Ntest),Phitest(Ntest),Bftest(Ntest,3))
Rtest = (/6.1d0/)
Ztest = (/0.5d0/)
Phitest = (/0.1d0/)
!Phitest(1) = 2._real64*pi

Write(*,*) '--------------------------------'
Write(*,*) '         TEST POINTS'
Write(*,*) 'R (m)',Rtest
Write(*,*) 'Z (m)',Ztest
Write(*,*) 'P (r)',Phitest
Write(*,*) '--------------------------------'

fname = '/home/jjl/BGRID/Bgrid_vac_narrow_mirror_90x82x65'
Call open_bgrid_fields(fname)

Call bfield_bgrid(Rtest,Phitest,Ztest,Ntest,Bftest,ierr)




Write(*,*) '                   Br        Bz        Bphi, ierr'
Do i=1,Ntest
  Write(*,'(a,3f18.8,i8)') 'And the answer: ',Bftest(i,1),Bftest(i,2),Bftest(i,3),ierr
Enddo

Call close_bgrid_fields

End program test_bfield_bgrid

