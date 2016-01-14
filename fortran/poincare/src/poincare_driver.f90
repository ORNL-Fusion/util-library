!-----------------------------------------------------------------------------
!+ Program to make poincare plots using jdl bfield library
!-----------------------------------------------------------------------------
Program poincare_driver
!
! Author(s): J.D. Lore - 02/20/2014
!
! Modules used:
Use kind_mod, Only: int32, real64
Use rmpcoil_module, Only : rmp_coil, rmp_coil_current, rmp_ncoil_pts,bfield_bs_cyl,build_d3d_ccoils_jl,&
  build_d3d_icoils_jl
Use M3DC1_routines_mod, Only : prepare_m3dc1_fields, m3dc1_factor, m3dc1_itime, m3dc1_toroidal_on_err,bfield_m3dc1, &
     m3dc1_field_type
Use g3d_module, Only : readg_g3d, bfield_geq_bicub, get_psi_bicub
Use math_geo_module, Only : rlinspace
Use fieldline_follow_mod, Only: bfield_method, follow_fieldlines_rzphi
Use util_routines, Only: get_psin_2d
Use phys_const, Only: pi
Implicit none

Logical, Parameter :: m3dc1_toroidal_off_grid = .true.
Integer(int32), parameter :: max_rmp_coils       = 12

Integer(int32) :: ntor_pts_coil, i, nstart_fl, nsteps
Integer(int32) :: iocheck, itest, ierr_b, ind_poin, ierr
Real(real64), Allocatable :: r1d(:), z1d(:),phistart_arr(:), fl_r(:,:), fl_z(:,:), fl_p(:,:), &
     psiout(:), psiNout(:), fl_r2(:,:), fl_z2(:,:), fl_p2(:,:)
Integer(int32), Allocatable :: ilg(:), fl_ierr(:), ilg2(:), fl_ierr2(:)
Real(real64) :: dphi_line, Adphirat
character(10) :: junk
Real(Kind=4)  :: tarray(2),tres,tres0
Logical :: calc_psiN_min = .false., follow_both_ways = .false.
!---------------------------------------------------------------------------
! Namelist variables:
Real(real64) :: &
 rmp_current(max_rmp_coils) = 0.d0, &
 m3dc1_scale_factor = 0.d0, &
 phistart_deg, rstart, rend, zstart, zend, dphi_line_deg

Character(Len=120) :: &
 rmp_type = 'none',  &
 gfile_name = 'none', &
 rmp_coil_type = 'none', &
 m3dc1_filename = 'none.none'

Integer(int32) :: &
 m3dc1_time = -1, &
 num_pts = 2, &
 ntransits = 1, &
 Nsym = 1

! Namelist files
Namelist / settings_nml / gfile_name, rmp_type, rmp_coil_type, m3dc1_filename, &
     m3dc1_scale_factor, rmp_current, m3dc1_time, phistart_deg, rstart, rend, zstart, zend, &
     num_pts, ntransits, dphi_line_deg, Nsym, calc_psiN_min, follow_both_ways

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
        Call build_d3d_ccoils_jl(rmp_current(1:12),ntor_pts_coil,rmp_coil,rmp_coil_current,rmp_ncoil_pts)
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
    m3dc1_factor = m3dc1_scale_factor
    m3dc1_toroidal_on_err = m3dc1_toroidal_off_grid
    m3dc1_field_type = 1
    Call prepare_m3dc1_fields(m3dc1_filename)
    Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
    Write(*,'(a,f12.4)') '---------> M3DC1 scale factor: ',m3dc1_factor    
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
    m3dc1_factor = m3dc1_scale_factor
    m3dc1_toroidal_on_err = m3dc1_toroidal_off_grid
    m3dc1_field_type = 0
    Call prepare_m3dc1_fields(m3dc1_filename)
    Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
    Write(*,'(a,f12.4)') '---------> M3DC1 scale factor: ',m3dc1_factor
    If (m3dc1_toroidal_on_err) Write(*,'(a)') '---------> M3DC1 fields will be set to B=Bt=1 off grid!'
  Case Default
    Write(*,*) 'Unknown rmp_type in poincare_driver!'
    Write(*,*) 'Current options are:'
    Write(*,*) '''g3d'''
    Write(*,*) '''g3d+rmpcoil'''
    Write(*,*) '''g3d+m3dc1'''
    Write(*,*) '''m3dc1_full_field'''
    Stop      
End Select


! 
! 3) Follow fls
! 


Allocate(r1d(num_pts),z1d(num_pts))
r1d = rlinspace(rstart,rend,num_pts)
z1d = rlinspace(zstart,zend,num_pts)

dphi_line = dphi_line_deg*pi/180.d0
nsteps = Floor(ntransits*2.d0*pi/Abs(dphi_line))

!If (bfield_method .eq. 0 ) Then
!  ind_poin = 1
!Else 
  Adphirat = Abs(360.d0/dphi_line_deg/Nsym)
  If (Adphirat - Real(Nint(Adphirat)) > 1.d-8) Then
    Write(*,*) 'Error!: 2*pi/Nsym must be an integer multiple of dphi_line_deg'
    Write(*,*) 'Exiting'
    Stop
  Endif
  ind_poin = Nint(Adphirat)
!Endif

Write(*,'(/1a,i0,a,i0,a)') 'Following ',num_pts,' fls for ',ntransits,' transits.'
Write(*,'(a,f12.3,a)') 'Toroidal step size is ',dphi_line_deg,' degrees'
Write(*,'(a,i0)') 'Poincare plot step size index is ',ind_poin
Write(*,'(a,f12.3)') 'Starting fl at phi = ',phistart_deg
Write(*,'(a,i0)') 'Number of steps = ',nsteps

nstart_fl = num_pts
Allocate(ilg(nstart_fl),fl_ierr(nstart_fl),phistart_arr(nstart_fl))
Allocate(fl_r(nstart_fl,nsteps+1),fl_z(nstart_fl,nsteps+1),fl_p(nstart_fl,nsteps+1))
fl_r = 0.d0; fl_z = 0.d0; fl_p = 0.d0
phistart_arr = phistart_deg*pi/180.d0

Call follow_fieldlines_rzphi(r1d,z1d,phistart_arr,nstart_fl, dphi_line,nsteps,fl_r,fl_z,fl_p,fl_ierr,ilg)

If (follow_both_ways) Then
  Allocate(ilg2(nstart_fl),fl_ierr2(nstart_fl))
  Allocate(fl_r2(nstart_fl,nsteps+1),fl_z2(nstart_fl,nsteps+1),fl_p2(nstart_fl,nsteps+1))
  fl_r2 = 0.d0; fl_z2 = 0.d0; fl_p2 = 0.d0 
  Call follow_fieldlines_rzphi(r1d,z1d,phistart_arr,nstart_fl,-dphi_line,nsteps,fl_r2,fl_z2,fl_p2,fl_ierr2,ilg2)
Endif
  
Open(99,file="poincare_output.out",status="unknown",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening output file'
  Stop 'Exiting: I/O Error in poincare_driver.f90'
Endif
Write(99,*) phistart_deg,nstart_fl,nsteps/ind_poin+1
!Write(*,*) phistart_deg,nstart_fl,nsteps/ind_poin+1

Do i = 1,nstart_fl
  Write(99,*) i
  Write(99,'(6f18.12)') fl_r(i,1:nsteps+1:ind_poin)
  Write(99,'(6f18.12)') fl_z(i,1:nsteps+1:ind_poin)
Enddo
Close(99)
Deallocate(r1d,z1d)

If (follow_both_ways) Then
  Open(99,file="poincare_output2.out",status="unknown",form="formatted",iostat=iocheck)
  If ( iocheck /= 0 ) Then
    Write(*,*) 'Error opening output file'
    Stop 'Exiting: I/O Error in poincare_driver.f90'
  Endif
  Write(99,*) phistart_deg,nstart_fl,nsteps/ind_poin+1
  
  Do i = 1,nstart_fl
    Write(99,*) i
    Write(99,'(6f18.12)') fl_r2(i,1:nsteps+1:ind_poin)
    Write(99,'(6f18.12)') fl_z2(i,1:nsteps+1:ind_poin)
  Enddo
  Close(99)
Endif

!
! Calculate min psi_N

!
If (calc_psiN_min) Then
  Write(*,*) 'Calculating minimum psi_N'
  Allocate(psiout(nsteps+1),psiNout(nsteps+1))

  Open(99,file="psiN_min_output.out",status="unknown",form="formatted",iostat=iocheck)
  If ( iocheck /= 0 ) Then
    Write(*,*) 'Error opening output file'
    Stop 'Exiting: I/O Error in poincare_driver.f90 (2)'
  Endif
  Write(99,*) nstart_fl

  Do i = 1,nstart_fl
    
    !Call get_psi_bicub(fl_r(i,:),fl_z(i,:),nsteps+1,psiout,psiNout,ierr)
    psiNout = get_psiN_2d(fl_r(i,:),fl_z(i,:),nsteps+1,ierr)
    
    Where (psiNout < 1.e-3) psiNout = 1000000.d0
    Write(99,*) Minval(psiNout)
  Enddo

  Close(99)
  Deallocate(psiout,psiNout)
Endif

If (follow_both_ways) Then
  ! Calculate min psi_N
  If (calc_psiN_min) Then
    Write(*,*) 'Calculating minimum psi_N'
    Allocate(psiout(nsteps+1),psiNout(nsteps+1))
    
    Open(99,file="psiN_min_output2.out",status="unknown",form="formatted",iostat=iocheck)
    If ( iocheck /= 0 ) Then
      Write(*,*) 'Error opening output file'
      Stop 'Exiting: I/O Error in poincare_driver.f90 (2)'
    Endif
    Write(99,*) nstart_fl
    
    Do i = 1,nstart_fl
      psiNout = get_psiN_2d(fl_r2(i,:),fl_z2(i,:),nsteps+1,ierr)
!      Call get_psi_bicub(fl_r2(i,:),fl_z2(i,:),nsteps+1,psiout,psiNout,ierr)      
      Where (psiNout < 1.e-3) psiNout = 1000000.d0
      Write(99,*) Minval(psiNout)
    Enddo
    
    Close(99)
    Deallocate(psiout,psiNout)
  Endif
Endif


Deallocate(ilg,fl_ierr,fl_r,fl_z,fl_p,phistart_arr)

If (follow_both_ways) Then
  Deallocate(ilg2,fl_ierr2,fl_r2,fl_z2,fl_p2)
Endif

Call Etime(tarray,tres)
Write(*,*) ' Poincare_driver took ',tres-tres0,' seconds'


End Program poincare_driver





