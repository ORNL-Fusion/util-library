Program compare_m3dc1_and_coils

  Use kind_mod
  Use m3dc1_routines_mod
  Use g3d_module
  Use rmpcoil_module
  Use biotsavart_module, only: bfield_bs_jdl, bfield_bs_cyl
  Use diiid_routines_mod
  Implicit None
  Integer, Parameter :: Ntest = 1
  Real(real64) :: bm3d(Ntest,3), rg(Ntest), zg(Ntest), pg(Ntest),Bg(Ntest,3)
  Real(real64), Dimension(Ntest) :: Br,Bphi,Bz
  Integer(int32) :: ierr_b, num_pts, i, m3dc1_nsets
  Character(250) :: m3dc1_filenames(2)
  Character(Len=100) :: gfilename
  Integer(int32) :: ntorpts, nsteps
  Type(g_type) :: g
  Type(coil_type) :: coil
  Real(Real64) :: taper(12), dphi
  Real(Real64), Parameter :: pi = 3.14159d0


  m3dc1_nsets = 2  ! MUST EQUAL LENGTH OF FILENAMES AND SCALE_FACTORS!!
  m3dc1_filenames(1) = '/home/jjl/m3dc1/164723/03059/n=3/eb1_2f_iu/C1.h5'
  m3dc1_filenames(2) = '/home/jjl/m3dc1/164723/03059/n=3/eb1_2f_il/C1.h5'

  m3dc1_itime = 0 ! 0:vac 1 response
  m3dc1_factors(1) = 3.69d0
  m3dc1_factors(2) = 3.69d0
  m3dc1_phases_deg(1) = 30.d0
  m3dc1_phases_deg(2) = 30.d0
  
  m3dc1_toroidal_on_err = .false.
  m3dc1_field_type = 1  ! 1=perturbed_only
  Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))

  gfilename = '/home/jjl/DIII-D/164723/g164723.03059_410'   
  Call readg_g3d(gfilename,g)
  
  taper = (/-2900.d0, 2900.d0, -2900.d0, 2900.d0, -2900.d0, 2900.d0, &
            -2900.d0, 2900.d0, -2900.d0, 2900.d0, -2900.d0, 2900.d0/)
  ntorpts = 6
  Call build_d3d_icoils_jl(coil,taper,ntorpts)


  Open(99,file='compare_m3dc1_and_coils.out')
  
  nsteps = 100
  dphi = 2.d0*pi/nsteps  
  pg = 0.d0
  do i = 1,nsteps
    rg = 2.2d0
    zg = .5d0
    pg = pg + dphi
    
    Write(*,*)
    Write(*,*) '----------------------------------------'  
    Write(*,*) 'Evaluating B at [R,Z,phi] = ',rg,zg,pg*180.d0/pi
    
    !  Call bfield_m3dc1_2d(rg,zg,Ntest,bm3d,ierr_b)     
    !  write(*,*) 'Bout [Br,Bz,Bphi] = ',bm3d
    
    Call bfield_m3dc1(rg,pg,zg,Ntest,bm3d,ierr_b)     
    write(*,*) 'Bout m3d [Br,Bz,Bphi] = ',bm3d
    
    
!    Call bfield_geq_bicub(g,rg,zg,Ntest,Bg,ierr_b)
!    write(*,*) 'Bout [Br,Bz,Bphi] = ',bg        
    
    Call bfield_bs_cyl(rg,pg,zg,Ntest,coil,Br,Bphi,Bz)
    write(*,*) 'Bout rmp [Br,Bz,Bphi] = ',Br,Bz,Bphi
    Write(99,*) rg,zg,pg,bm3d,Br,Bz,Bphi
    Write(*,*) '----------------------------------------'
  Enddo
  Close(99)

End Program Compare_M3dc1_And_Coils

  
