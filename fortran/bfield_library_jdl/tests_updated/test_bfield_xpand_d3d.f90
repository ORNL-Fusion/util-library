!-----------------------------------------------------------------------------
!
! -- JDL
!
Program test_bfield_xpand_d3d
Use kind_mod, Only: int32, real64
Use xpand_module
Use phys_const, Only: pi
Implicit None

! Local variables (scalar)
Real(real64),Allocatable :: Bftest(:,:),Rtest(:),Ztest(:)
! Local variables (array)
Character(Len=100) :: gfilename, fname, run_path

Integer :: ierr_b,Ntest,ierr,ifield_type

Integer(int32) :: i
Real(real64), allocatable, Dimension(:) :: Bx,By,Bz,P_x,P_y,P_z,Phitest,Br,Bphi

!- End of header -------------------------------------------------------------

Write(*,*) 'Test call to bfield_xpand'
Ntest = 1
Allocate(Rtest(Ntest),Ztest(Ntest),Phitest(Ntest),Bftest(Ntest,3))
Rtest = (/1.8d0/)
Ztest = (/0.5d0/)
!Phitest = (/100.1d0/)
Phitest(1) = 2._real64*pi

Write(*,*) '--------------------------------'
Write(*,*) '         TEST POINTS'
Write(*,*) 'R (m)',Rtest
Write(*,*) 'Z (m)',Ztest
Write(*,*) 'P (r)',Phitest
Write(*,*) '--------------------------------'

fname = '/home/jjl/XPAND/164723/3059/xpand_164723_3059.dat'
Call read_xpand_field_file(fname)

ifield_type = 1

Call bfield_xpand(Rtest,Phitest,Ztest,Ntest,Bftest,ierr,ifield_type)




Write(*,*) '                   Br        Bz        Bphi'
Do i=1,Ntest
  Write(*,'(a,3f18.8,i0)') 'And the answer: ',Bftest(i,1),Bftest(i,2),Bftest(i,3),ierr_b
Enddo

Call close_xpand_fields

End program test_bfield_xpand_d3d

