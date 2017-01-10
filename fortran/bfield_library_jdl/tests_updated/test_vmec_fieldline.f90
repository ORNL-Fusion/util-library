Program test_vmec_bfield
  Use kind_mod
  Use setup_bfield_module
  Use fieldline_follow_mod, Only: follow_fieldlines_rzphi
  Use bfield
  Use phys_const
  Implicit None
  Character(Len=240) :: file_name
  Real(real64) :: rad(1),phi(1),zz(1),Bout(1,3)
  Real(real64) :: extcur(100)
  Integer(int32) :: ierr, METHOD, i

  Integer(int32) :: nstart_fl, nsteps, ntransits
  Real(real64), Allocatable :: r1d(:), z1d(:),phistart_arr(:), fl_r(:,:), fl_z(:,:), fl_p(:,:)
  Integer(int32), Allocatable :: ilg(:), fl_ierr(:)
  Real(real64) :: dphi_line, dphi_line_deg, rstart(1), zstart(1), phistart(1)

  file_name = "/home/jjl/VMEC_RUNS/W7X_2016/coils/coils.w7x"

  METHOD = 2

  EXTCUR(1) =  11923.d0 !Nonplanar 1-5
  EXTCUR(2) =  11799.d0
  EXTCUR(3) =  12048.d0
  EXTCUR(4) =  13290.d0
  EXTCUR(5) =  13414.d0
  EXTCUR(6) =   8197.d0 ! A
  EXTCUR(7) =  -2981.d0 ! B
  EXTCUR(8) =  -2515.d0 ! Sweep 1
  EXTCUR(9) =   2515.d0 ! Sweep 2 (repeated)
  EXTCUR(10) = -2515.d0
  EXTCUR(11) =  2515.d0
  EXTCUR(12) = -2515.d0
  EXTCUR(13) =  2515.d0
  EXTCUR(14) = -2515.d0
  EXTCUR(15) =  2515.d0
  EXTCUR(16) = -2515.d0
  EXTCUR(17) =  2515.d0
  EXTCUR(18) =  0.d0 ! Trim coils 1-5
  EXTCUR(19) =  0.d0
  EXTCUR(20) =  0.d0
  EXTCUR(21) =  0.d0
  EXTCUR(22) =  0.d0

  Write(*,*) Trim(adjustl(file_name))

  rad = 5.5
  phi = 3.14159d0/2.d0
  zz = 0.25

  vmec_coils_file = Trim(file_name)
  vmec_extcur_set(1:22) = EXTCUR(1:22)
  
  Call setup_bfield_vmec_coils
  Call calc_B_rzphi_general(bfield,rad,zz,phi,1,Bout(1,1),Bout(1,2),Bout(1,3))
  Write(*,*) Bout
  
  ntransits = 1
  dphi_line_deg = 0.5d0

  rstart = 5.9d0
  zstart = 0.1d0
  phistart = 0.d0
  
  dphi_line = -dphi_line_deg*pi/180.d0
  nsteps = Floor(ntransits*2.d0*pi/Abs(dphi_line))
  
  nstart_fl = 1
  Allocate(ilg(nstart_fl),fl_ierr(nstart_fl))
  Allocate(fl_r(nstart_fl,nsteps+1),fl_z(nstart_fl,nsteps+1),fl_p(nstart_fl,nsteps+1))
  fl_r = 0.d0; fl_z = 0.d0; fl_p = 0.d0

  Write(*,*) 'following fl from [r,z,phi_deg] = ',rstart,zstart,phistart*180.d0/pi
  Call follow_fieldlines_rzphi(bfield,rstart,zstart,phistart,nstart_fl, dphi_line,nsteps,fl_r,fl_z,fl_p,fl_ierr,ilg)

  Open(99,file='vmec_fieldline.out')
  Do i = 1,nsteps+1
    Write(99,*) fl_r(1,i),fl_z(1,i),fl_p(1,i)
  Enddo

End Program test_vmec_bfield
