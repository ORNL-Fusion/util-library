!-----------------------------------------------------------------------------
!
! Tests bfield_geq_bicub, bfield_bs_jdl, bfield_bs_cyl, build_d3d_icoils_jl
! -- JDL
!
Program test_bfield_rmp_d3d
  Use kind_mod, Only: int32, real64
Use g3d_module
Use rmpcoil_module
Use biotsavart_module, only: bfield_bs_jdl, bfield_bs_cyl
Use diiid_routines_mod
Implicit None

Real(real64),Allocatable :: Btest(:,:),Rtest(:),Ztest(:)
Character(Len=100) :: gfilename
Integer(int32) :: ierr_b,Ntest,ntorpts
Type(g_type) :: g
Type(coil_type) :: coil
Real(Real64) :: taper(12)
Real(real64), allocatable, Dimension(:) :: Bx,By,Bz,P_x,P_y,P_z,Phitest,Br,Bphi

!- End of header -------------------------------------------------------------


!gfilename = './g160884.03014_251'
gfilename = '/home/jjl/DIII-D/164723/g164723.03059_410'

Call readg_g3d(gfilename,g)

Write(*,*) 'Test call to bfield_geq_bicub'
Ntest = 1
Allocate(Rtest(Ntest),Ztest(Ntest),Phitest(Ntest),Btest(Ntest,3))
Rtest = (/2.2d0/)
Ztest = (/0.d0/)
Phitest = (/0.d0/)

!Rtest = (/1.8d0,1.8d0/)
!Ztest = (/-0.5d0,0.5d0/)
!Phitest = (/0.0d0,0.1d0/)

Write(*,*) '--------------------------------'
Write(*,*) '         TEST POINTS'
Write(*,*) 'R (m)',Rtest
Write(*,*) 'Z (m)',Ztest
Write(*,*) 'P (r)',Phitest
Write(*,*) '--------------------------------'

Call bfield_geq_bicub(g,Rtest,Ztest,Ntest,Btest,ierr_b)
Write(*,*) '                   Br        Bz        Bphi'
Write(*,'(a,3f10.4,i0)') 'And the answer: ',Btest(1,:),ierr_b
if (Ntest == 2) Then
  Write(*,'(a,3f10.4,i0)') 'And the answer: ',Btest(2,:),ierr_b
Endif

Write(*,*) 'building I coils'

!taper = (/-3740.d0, 3863.d0,-3720.d0, 3855.d0,-3718.d0, 3858.d0, &
!           3862.d0,-3791.d0, 3884.d0,-3854.d0, 3923.d0,-3847.d0/)

taper = (/-2900.d0, 2900.d0, -2900.d0, 2900.d0, -2900.d0, 2900.d0, &
          -2900.d0, 2900.d0, -2900.d0, 2900.d0, -2900.d0, 2900.d0/)


Write(*,*) 'Taper:',taper
ntorpts = 6
Call build_d3d_icoils_jl(coil,taper,ntorpts)
!ntorpts = 4
!Call build_d3d_ccoils_jl(coil,taper(1:6),ntorpts)

Write(*,*) 'Test call to bfield_bs_jdl'
Allocate(P_x(Ntest),P_y(Ntest),P_z(Ntest))
Allocate(Bx(Ntest),By(Ntest),Bz(Ntest))
Allocate(Br(Ntest),Bphi(Ntest))
P_x = Rtest*cos(Phitest)
P_y = Rtest*sin(Phitest)
P_z = Ztest
Call bfield_bs_jdl(P_x,P_y,P_z,Ntest,coil,Bx,By,Bz)

Write(*,*) 'Bx',Bx
Write(*,*) 'By',By
Write(*,*) 'Bz',Bz


Write(*,*) 'Test call to bfield_bs_cyl'
write(*,*) '  '
Call bfield_bs_cyl(Rtest,Phitest,Ztest,Ntest,coil,Br,Bphi,Bz)

Write(*,*) 'Br',Br
Write(*,*) 'Bphi',Bphi
Write(*,*) 'Bz',Bz


End program test_bfield_rmp_d3d

