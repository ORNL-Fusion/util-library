!-----------------------------------------------------------------------------
!+ Program to make poincare plots using jdl bfield library
!-----------------------------------------------------------------------------
Program poincare_driver_2d
!
! Author(s): J.D. Lore - 02/20/2014
!
! Modules used:
Use kind_mod, Only: int32, real64
Use rmpcoil_module, Only : rmp_coil, rmp_coil_current, rmp_ncoil_pts,bfield_bs_cyl,build_d3d_ccoils_jl,&
  build_d3d_icoils_jl
Use M3DC1_routines_mod, Only : prepare_m3dc1_fields, m3dc1_factors, m3dc1_itime, m3dc1_toroidal_on_err,bfield_m3dc1, &
     m3dc1_field_type
Use g3d_module, Only : readg_g3d, bfield_geq_bicub, get_psi_bicub
Use math_geo_module, Only : rlinspace
Use fieldline_follow_mod, Only: bfield_method, follow_fieldlines_rzphi
Use util_routines, Only: get_psin_2d
Use phys_const, Only: pi
Use ipec_module, Only: open_ipec_fields
Use xpand_module, Only: open_xpand_fields
Use util_routines, Only: num_lines_file
Use math_geo_module, Only: inside_poly
Implicit none

Logical, Parameter :: m3dc1_toroidal_off_grid = .true.
Integer(int32), parameter :: max_rmp_coils       = 12

Integer(int32) :: ntor_pts_coil, i, nstart_fl, nsteps, j, i1_lc, i2_lc, inside1, inside2
Integer(int32) :: iocheck, ind_poin, ierr, num_vv
Real(real64), Allocatable :: r1d(:), z1d(:),phistart_arr(:), fl_r(:,:), fl_z(:,:), fl_p(:,:), &
     psiout(:), psiNout(:), fl_r2(:,:), fl_z2(:,:), fl_p2(:,:), r1d_tmp(:), z1d_tmp(:), psiNmin(:,:), &
     rp(:), zp(:), lc1(:), lc2(:), lc(:)
Integer(int32), Allocatable :: ilg(:), fl_ierr(:), ilg2(:), fl_ierr2(:)
Real(real64) :: dphi_line, Adphirat
Real(Kind=4)  :: tarray(2),tres,tres0
Logical :: calc_psiN_min = .false., follow_both_ways = .false., found1, found2
Integer(int32), parameter :: max_m3dc1_files = 10
!---------------------------------------------------------------------------
! Namelist variables:
Real(real64) :: &
 rmp_current(max_rmp_coils) = 0.d0, &
 m3dc1_scale_factors(max_m3dc1_files) = 0.d0, &
 phistart_deg, rstart, rend, zstart, zend, dphi_line_deg

Character(Len=120) :: &
 rmp_type = 'none',  &
 gfile_name = 'none', &
 rmp_coil_type = 'none', &
 m3dc1_filenames(max_m3dc1_files), &
 ipec_run_path = 'none', &
 xpand_fname = 'none'

Integer(int32) :: &
 m3dc1_time = -1, &
 num_pts = 2, &
 ntransits = 1, &
 Nsym = 1, &
 m3dc1_nsets = 0, &
 ipec_field_eval_type = -1, &
 xpand_field_eval_type = -1

! Namelist files
Namelist / settings_nml / gfile_name, rmp_type, rmp_coil_type, m3dc1_filenames, &
     m3dc1_scale_factors, rmp_current, m3dc1_time, phistart_deg, rstart, rend, zstart, zend, &
     num_pts, ntransits, dphi_line_deg, Nsym, calc_psiN_min, follow_both_ways, m3dc1_nsets, &
     ipec_run_path, ipec_field_eval_type, xpand_fname, xpand_field_eval_type

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

Call Etime(tarray,tres0)
Write(*,'(/1a)') "-------------------------------------------------------------"
Write(*,'(a)') " Starting Poincare driver"

!
! 1) Read namelist
!
! Read namelist file 
Open(99,file="poincare_settings.nml",status="old",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening namelist file'
  Stop 'Exiting: I/O Error in poincare_driver.f90'
Endif
Read(99,nml=settings_nml)
Close(99)

!
! 2) Magnetic equilibrium data
!
! Setup rmp field
Select Case (rmp_type)
  Case ('g3d')
    Write(*,'(a)') '-----> BFIELD METHOD IS G3D'
    bfield_method = 0
    Call readg_g3d(gfile_name)
  Case ('g3d+rmpcoil')
    Write(*,'(a)') '-----> BFIELD METHOD IS G3D+RMPCOIL'
    bfield_method = 1
    Call readg_g3d(gfile_name)
    ! load rmp coil    
    Select Case (rmp_coil_type)
      Case ('d3d_ccoils')
        Write(*,'(a)') '-----> Rmp coils are DIII-D C coils'
        Write(*,'(a,6f12.3)') ' Coil currents: ', rmp_current(1:6)
        ntor_pts_coil = 4
        rmp_ncoil_pts = 6*(2*ntor_pts_coil+1)
        Allocate(rmp_coil(rmp_ncoil_pts,3))
        Allocate(rmp_coil_current(rmp_ncoil_pts))
        Call build_d3d_ccoils_jl(rmp_current(1:6),ntor_pts_coil,rmp_coil,rmp_coil_current,rmp_ncoil_pts)
      Case ('d3d_icoils')
        Write(*,'(a)') '----->  Rmp coils are DIII-D I coils'
        Write(*,'(a,12f12.3)') ' Coil currents: ', rmp_current(1:12)
        ntor_pts_coil = 6
        rmp_ncoil_pts = 12*(2*ntor_pts_coil+1)
        Allocate(rmp_coil(rmp_ncoil_pts,3))
        Allocate(rmp_coil_current(rmp_ncoil_pts))
        Call build_d3d_icoils_jl(rmp_current(1:12),ntor_pts_coil,rmp_coil,rmp_coil_current,rmp_ncoil_pts)
      Case Default
        Write(*,*) 'Error!  Unknown rmp_coil_type'
        Write(*,*) 'Current options are:'
        Write(*,*) '''d3d_ccoils'''        
        Write(*,*) '''d3d_icoils'''        
        Stop
      End Select
  Case ('g3d+m3dc1')
    Write(*,'(a)') '-----> BFIELD METHOD IS G3D+M3DC1'
    bfield_method = 3
    Call readg_g3d(gfile_name)
    ! Setup field    
    If (m3dc1_time .eq. -1) Then
      Write(*,*) 'Error: Bfield type of M3DC1 is set but m3dc1_time is not. Exiting.'
      Stop
    Endif
    m3dc1_itime = m3dc1_time           ! Variables in module
    m3dc1_factors = m3dc1_scale_factors
    m3dc1_toroidal_on_err = m3dc1_toroidal_off_grid
    m3dc1_field_type = 1
    Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
    Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
    Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)
    If (m3dc1_toroidal_on_err) Write(*,'(a)') '---------> M3DC1 fields will be set to B=Bt=1 off grid!'
  Case ('m3dc1_full_field')
    Write(*,'(a)') '-----> BFIELD METHOD IS M3DC1_FULL_FIELD'
    bfield_method = 4
    ! Setup field    
    If (m3dc1_time .eq. -1) Then
      Write(*,*) 'Error: Bfield type of M3DC1 is set but m3dc1_time is not. Exiting.'
      Stop
    Endif
    m3dc1_itime = m3dc1_time           ! Variables in module
    m3dc1_factors = m3dc1_scale_factors
    m3dc1_toroidal_on_err = m3dc1_toroidal_off_grid
    m3dc1_field_type = 0
    Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
    Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
    Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)
    If (m3dc1_toroidal_on_err) Write(*,'(a)') '---------> M3DC1 fields will be set to B=Bt=1 off grid!'
  Case ('m3dc1_as')
    Write(*,'(a)') '-----> BFIELD METHOD IS m3dc1_as'
    bfield_method = 5
    ! Setup field    
    If (m3dc1_time .eq. -1) Then
      Write(*,*) 'Error: Bfield type of M3DC1 is set but m3dc1_time is not. Exiting.'
      Stop
    Endif
    m3dc1_itime = m3dc1_time           ! Variables in module
    m3dc1_factors = m3dc1_scale_factors
    m3dc1_toroidal_on_err = m3dc1_toroidal_off_grid
    m3dc1_field_type = 0
    Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
    Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
    Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)
    If (m3dc1_toroidal_on_err) Write(*,'(a)') '---------> M3DC1 fields will be set to B=Bt=1 off grid!'
  Case ('ipec')
    Write(*,'(a)') '-----> BFIELD METHOD IS IPEC'
    Write(*,'(a,i0)') '-----> ipec_field_eval_type is: ',ipec_field_eval_type
    If (ipec_field_eval_type .eq. 0) Then
      Write(*,'(a)') '-----> Evaluating IPEC fields as EQUILIBRIUM ONLY!'
      bfield_method = 7
    Elseif (ipec_field_eval_type .eq. 1) Then
      Write(*,'(a)') '-----> Evaluating IPEC fields as EQ + VACUUM!'
      bfield_method = 8
    Elseif (ipec_field_eval_type .eq. 2) Then
      Write(*,'(a)') '-----> Evaluating IPEC fields as EQ + PERT!'
      bfield_method = 9
    Else
      Stop "Did not recognize ipec_field_eval_type"
    Endif
    Call readg_g3d(gfile_name)
    Call open_ipec_fields(ipec_run_path)
  Case ('xpand')
    Write(*,'(a)') '-----> BFIELD METHOD IS XPAND'
    Write(*,'(a,i0)') '-----> xpand_field_eval_type is: ',xpand_field_eval_type
    If (xpand_field_eval_type .eq. 0) Then
      Write(*,'(a)') '-----> Evaluating XPAND fields as PERTURBED!'
      bfield_method = 10
    Elseif (xpand_field_eval_type .eq. 1) Then
      Write(*,'(a)') '-----> Evaluating XPAND fields as VACUUM!'
      bfield_method = 11
    Else
      Stop "Did not recognize xpand_field_eval_type"
    Endif
    Call readg_g3d(gfile_name)
    Call open_xpand_fields(xpand_fname)    
  Case Default
    Write(*,*) 'Unknown rmp_type in poincare_driver!'
    Write(*,*) 'Current options are:'
    Write(*,*) '''g3d'''
    Write(*,*) '''g3d+rmpcoil'''
    Write(*,*) '''g3d+m3dc1'''
    Write(*,*) '''m3dc1_full_field'''
    Write(*,*) '''m3dc1_as'''
    Write(*,*) '''ipec'''
    Write(*,*) '''xpand'''
    Stop      
End Select


! 
! 3) Follow fls
! 


Allocate(r1d(num_pts**2),z1d(num_pts**2),z1d_tmp(num_pts),r1d_tmp(num_pts))

r1d_tmp = rlinspace(rstart,rend,num_pts)
z1d_tmp = rlinspace(zstart,zend,num_pts)

Do i = 1,num_pts
  r1d((i-1)*num_pts+1:i*num_pts) = r1d_tmp
  z1d((i-1)*num_pts+1:i*num_pts) = z1d_tmp(i)
Enddo

dphi_line = dphi_line_deg*pi/180.d0
nsteps = Floor(ntransits*2.d0*pi/Abs(dphi_line))

Adphirat = Abs(360.d0/dphi_line_deg/Nsym)
If (Adphirat - Real(Nint(Adphirat)) > 1.d-8) Then
  Write(*,*) 'Error!: 2*pi/Nsym must be an integer multiple of dphi_line_deg'
  Write(*,*) 'Exiting'
  Stop
Endif
ind_poin = Nint(Adphirat)

Write(*,'(/1a,i0,a,i0,a)') 'Following ',num_pts,' fls for ',ntransits,' transits.'
Write(*,'(a,f12.3,a)') 'Toroidal step size is ',dphi_line_deg,' degrees'
Write(*,'(a,i0)') 'Poincare plot step size index is ',ind_poin
Write(*,'(a,f12.3)') 'Starting fl at phi = ',phistart_deg
Write(*,'(a,i0)') 'Number of steps = ',nsteps

nstart_fl = num_pts**2
Allocate(ilg(nstart_fl),fl_ierr(nstart_fl),phistart_arr(nstart_fl))
Allocate(fl_r(nstart_fl,nsteps+1),fl_z(nstart_fl,nsteps+1),fl_p(nstart_fl,nsteps+1))
fl_r = 0.d0; fl_z = 0.d0; fl_p = 0.d0
phistart_arr = phistart_deg*pi/180.d0

Call follow_fieldlines_rzphi(r1d,z1d,phistart_arr,nstart_fl, dphi_line,nsteps,fl_r,fl_z,fl_p,fl_ierr,ilg)

Allocate(ilg2(nstart_fl),fl_ierr2(nstart_fl))
Allocate(fl_r2(nstart_fl,nsteps+1),fl_z2(nstart_fl,nsteps+1),fl_p2(nstart_fl,nsteps+1))
fl_r2 = 0.d0; fl_z2 = 0.d0; fl_p2 = 0.d0 
Call follow_fieldlines_rzphi(r1d,z1d,phistart_arr,nstart_fl,-dphi_line,nsteps,fl_r2,fl_z2,fl_p2,fl_ierr2,ilg2)
  

!
! Calculate min psi_N
!

Allocate(psiNmin(2,nstart_fl))

If (calc_psiN_min) Then
  Write(*,*) 'Calculating minimum psi_N'
  Allocate(psiout(nsteps+1),psiNout(nsteps+1))

  Do i = 1,nstart_fl
    
    !Call get_psi_bicub(fl_r(i,:),fl_z(i,:),nsteps+1,psiout,psiNout,ierr)
    psiNout = get_psiN_2d(fl_r(i,:),fl_z(i,:),nsteps+1,ierr)
    
    Where (psiNout < 1.e-3) psiNout = 1000000.d0
    psiNmin(1,i) = Minval(psiNout)
  Enddo

  Deallocate(psiout,psiNout)
Endif

If (follow_both_ways) Then
  ! Calculate min psi_N
  If (calc_psiN_min) Then
    Write(*,*) 'Calculating minimum psi_N'
    Allocate(psiout(nsteps+1),psiNout(nsteps+1))    
    
    Do i = 1,nstart_fl
      psiNout = get_psiN_2d(fl_r2(i,:),fl_z2(i,:),nsteps+1,ierr)
!      Call get_psi_bicub(fl_r2(i,:),fl_z2(i,:),nsteps+1,psiout,psiNout,ierr)      
      Where (psiNout < 1.e-3) psiNout = 1000000.d0
      psiNmin(2,i) = Minval(psiNout)
    Enddo

    Deallocate(psiout,psiNout)
  Endif
Endif


! read vvfile (in mm)
num_vv = num_lines_file('vvfile.ogr')
Open(98,file='vvfile.ogr',status='old',iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening vvfile.ogr file'
  Stop 'Exiting: I/O Error in poincare_driver.f90'
Endif


Allocate(rp(num_vv),zp(num_vv))
Do i = 1,num_vv
  Read(98,*) rp(i),zp(i)
Enddo
rp = rp/1.d3
zp = zp/1.d3
Close(98)

! Calculate L_c
Allocate(lc(nstart_fl),lc1(nstart_fl),lc2(nstart_fl))
lc = 0.d0
lc1 = 0.d0
lc2 = 0.d0
Do i = 1,nstart_fl
  write(*,*) 'Lc checking line ',i,' of ',nstart_fl
  found1 = .false.
  found2 = .false.
  i1_lc = 0
  i2_lc = 0
  Do j = 1,ilg(i)
    if (.NOT. found1) Then
      inside1 = inside_poly(fl_r(i,j),fl_z(i,j),rp,zp,num_vv)
!      if (j == 1) Then
!        Write(*,*) 'Starts inside?',inside1
!        write(*,*) fl_r(i,j),fl_z(i,j)
!      Endif
      
      if (inside1 .eq. 0) Then
        found1 = .true.
        lc1(i) = Sum( Sqrt( &
             fl_r(i,1:j-1)**2 + &
             fl_r(i,2:j)**2 + &        
             (fl_z(i,2:j) - fl_z(i,1:j-1))**2 - &
             2.d0*fl_r(i,1:j-1)*fl_r(i,2:j)*cos(fl_p(i,2:j) - fl_p(i,1:j-1)) ) )             
      Endif
    endif
  Enddo
  Do j = 1,ilg2(i)
    if (.NOT. found2) Then
      inside2 = inside_poly(fl_r2(i,j),fl_z2(i,j),rp,zp,num_vv)
!      if (j == 1) Write(*,*) 'Starts inside (rev)?',inside2
      if (inside2 .eq. 0) Then
        found2 = .true.
        lc2(i) = Sum( Sqrt( &
             fl_r2(i,1:j-1)**2 + &
             fl_r2(i,2:j)**2 + &        
             (fl_z2(i,2:j) - fl_z2(i,1:j-1))**2 - &
             2.d0*fl_r2(i,1:j-1)*fl_r2(i,2:j)*cos(fl_p2(i,2:j) - fl_p2(i,1:j-1)) ) )                     
      Endif
    Endif
  Enddo
  lc(i) = lc1(i) + lc2(i)
Enddo


Open(99,file="poincare_output_2d.out",status="unknown",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening output file'
  Stop 'Exiting: I/O Error in poincare_driver.f90'
Endif
Write(99,*) phistart_deg,nstart_fl,nsteps/ind_poin+1

Do i = 1,nstart_fl 
!  Write(99,'(i0,4f18.12)') i, fl_r(i,1), fl_z(i,1), Minval(psiNmin(1:2,i)), lc(i)
  Write(99,*) i, fl_r(i,1), fl_z(i,1), Minval(psiNmin(1:2,i)), lc(i)  
Enddo
Close(99)
Deallocate(r1d,z1d)



Deallocate(ilg,fl_ierr,fl_r,fl_z,fl_p,phistart_arr)
Deallocate(ilg2,fl_ierr2,fl_r2,fl_z2,fl_p2)



Call Etime(tarray,tres)
Write(*,*) ' Poincare_driver took ',tres-tres0,' seconds'


End Program poincare_driver_2d





