Program test_vmec_bfield
  Use kind_mod
  Use vmec_routines_mod
  Use setup_bfield_module
  Use bfield
  Use math_geo_module, Only : rlinspace
  Use phys_const
  Implicit None
  
  Character(Len=240) :: file_name
  Real(real64), allocatable :: rad(:),phi(:),zz(:),Bout(:,:)
  Real(real64) :: extcur(100)
  Integer(int32) :: ierr, METHOD, npts, i

  
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

  npts = 1
  allocate(rad(npts),phi(npts),zz(npts),bout(npts,3))

  rad(:) = 6.1d0
  zz(:) = 0.1d0
  phi(1) = 2.5d0
!  phi = rlinspace(0.d0,2.d0*pi,npts)
  write(*,*) phi
  
  If (METHOD == 1) Then
    Call read_vmec_coils_file(file_name)
    vmec_extcur(1:22) = extcur(1:22)    
    Call bfield_vmec_coils(rad,phi,zz,npts,Bout,ierr)  
  Else
    vmec_coils_file = Trim(file_name)
    vmec_extcur_set(1:22) = EXTCUR(1:22)
    
    Call setup_bfield_vmec_coils
    Call calc_B_rzphi_general(bfield,rad,zz,phi,npts,Bout(:,1),Bout(:,2),Bout(:,3))
  Endif

  
  
  Open(99,file='vmec_bfield.out')
  Do i = 1,npts
    Write(99,*) Bout(i,1:3)
    if (npts .lt. 10) write(*,*) Bout(i,1:3)
  Enddo
End Program test_vmec_bfield
