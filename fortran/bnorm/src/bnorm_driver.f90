!-----------------------------------------------------------------------------
!+ Program to make bnorm plots
!-----------------------------------------------------------------------------
Program bnorm_driver
!
! Author(s): J.D. Lore - 10/13/2016
!
! Modules used:
Use kind_mod, Only: int32, real64
Use rmpcoil_module, Only : build_d3d_ccoils_jl, build_d3d_icoils_jl
Use M3DC1_routines_mod, Only : prepare_m3dc1_fields, m3dc1_factors, m3dc1_itime, m3dc1_toroidal_on_err,bfield_m3dc1, &
     m3dc1_field_type
Use g3d_module, Only : readg_g3d, get_psi_derivs_bicub
Use bnorm_routines, Only : get_pest_coords
Use math_geo_module, Only : rlinspace
Use phys_const, Only: pi
Use ipec_module, Only: open_ipec_fields
Use xpand_module, Only: open_xpand_fields
Use bfield, Only : bfield_type, coil_type, g_type, calc_b_rzphi_general, set_bfield_pert_only, reset_bfield
Implicit none

Logical, Parameter :: m3dc1_toroidal_off_grid = .true.
Integer(int32), parameter :: max_rmp_coils       = 12

Integer(int32) :: i, k, mind
Integer(int32) :: ierr,iocheck
Real(Kind=4)  :: tarray(2),tres,tres0
Integer(int32), parameter :: max_m3dc1_files = 10
Type(bfield_type) :: bfield
Type(g_type) :: g
Type(coil_type) :: coil
Integer(int32), Allocatable :: marr(:)
Real(real64) :: dphi, dtheta
Real(real64), Allocatable :: pnwant(:), phi(:), theta(:)
Real(real64), Allocatable :: rpest(:,:), zpest(:,:), jpest(:,:),dpsidr(:),dpsidz(:),psiout(:)
Real(real64), Allocatable :: rn(:),zn(:),br(:),bz(:),bphi(:),bnorm(:), area(:),phiarr(:),alpha_mn(:)
Real(real64), Allocatable :: br_c(:,:), br_s(:,:), Br_mn(:,:)
!---------------------------------------------------------------------------
! Namelist variables:
Real(real64) :: &
 rmp_current(max_rmp_coils) = 0.d0, &
 m3dc1_scale_factors(max_m3dc1_files) = 0.d0, &
 pnmin = 0.d0, &
 pnmax = 0.d0 

Character(Len=120) :: &
 rmp_type                         = 'none',  &
 gfile_name                       = 'none', &
 rmp_coil_type                    = 'none', &
 m3dc1_filenames(max_m3dc1_files) = 'none', &
 ipec_run_path                    = 'none', &
 xpand_fname                      = 'none'

Integer(int32) :: &
 m3dc1_time = -1, &
 m3dc1_nsets = 0, &
 ipec_field_eval_type = -1, &
 xpand_field_eval_type = -1, &
 ipec_itype = 0, &
 nres = 0, &
 ntheta = 0, &
 nphi = 0, &
 mmax = 0, &
 numpn = 0

! Namelist files
Namelist / settings_nml / &
     ! driver stuff
     nres, ntheta, nphi, mmax, pnmin, pnmax, numpn, &
     ! Bfield stuff
     rmp_type, &
     gfile_name, &
     rmp_coil_type, rmp_current, &
     m3dc1_filenames, m3dc1_scale_factors, m3dc1_time, m3dc1_nsets, &
     ipec_run_path, ipec_field_eval_type, ipec_itype, &
     xpand_fname, xpand_field_eval_type

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

Call Etime(tarray,tres0)
Write(*,'(/1a)') "-------------------------------------------------------------"
Write(*,'(a)') " Starting bnorm driver"

!
! 1) Read namelist
!
! Read namelist file 
Open(99,file="bnorm_settings.nml",status="old",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening namelist file'
  Stop 'Exiting: I/O Error in bnorm_driver.f90'
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
    bfield%method    = 0
    bfield%method_2d = 0
    bfield%method_pert = -1
    Call readg_g3d(gfile_name,g)
    bfield%g = g
  Case ('g3d+rmpcoil')
    Write(*,'(a)') '-----> BFIELD METHOD IS G3D+RMPCOIL'
    bfield%method = 1
    bfield%method_2d = 0
    bfield%method_pert = 6
    Call readg_g3d(gfile_name,g)
    bfield%g = g
    ! load rmp coil    
    Select Case (rmp_coil_type)
      Case ('d3d_ccoils')
        Write(*,'(a)') '-----> Rmp coils are DIII-D C coils'
        Write(*,'(a,6f12.3)') ' Coil currents: ', rmp_current(1:6)
        Call build_d3d_ccoils_jl(coil,rmp_current(1:6))
      Case ('d3d_icoils')
        Write(*,'(a)') '----->  Rmp coils are DIII-D I coils'
        Write(*,'(a,12f12.3)') ' Coil currents: ', rmp_current(1:12)
        Call build_d3d_icoils_jl(coil,rmp_current(1:12))
      Case Default
        Write(*,*) 'Error!  Unknown rmp_coil_type'
        Write(*,*) 'Current options are:'
        Write(*,*) '''d3d_ccoils'''        
        Write(*,*) '''d3d_icoils'''        
        Stop
      End Select
      bfield%coil = coil
  Case ('g3d+m3dc1')
    Write(*,'(a)') '-----> BFIELD METHOD IS G3D+M3DC1'
    bfield%method    = 3
    bfield%method_2d = 0
    bfield%method_pert = 4
    Call readg_g3d(gfile_name,g)
    bfield%g = g
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
    bfield%method    = 4
    bfield%method_2d = 5
    bfield%method_pert = -1
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
    Stop "This method not valid for bnorm calculation"
  Case ('m3dc1_as')
    Write(*,'(a)') '-----> BFIELD METHOD IS m3dc1_as'
    bfield%method    = 5
    bfield%method_2d = 5
    bfield%method_pert = -1
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
      bfield%method = 7
      bfield%method_2d = 0
      bfield%method_pert = -1
    Elseif (ipec_field_eval_type .eq. 1) Then
      Write(*,'(a)') '-----> Evaluating IPEC fields as EQ + VACUUM!'
      bfield%method = 8
      bfield%method_2d = 0
      bfield%method_pert = -1
    Elseif (ipec_field_eval_type .eq. 2) Then
      Write(*,'(a)') '-----> Evaluating IPEC fields as EQ + PERT!'
      bfield%method = 9
      bfield%method_2d = 0
      bfield%method_pert = -1
    Else
      Stop "Did not recognize ipec_field_eval_type"
    Endif
    Call readg_g3d(gfile_name,g)
    bfield%g = g
    Call open_ipec_fields(ipec_run_path,ipec_itype)
  Case ('xpand')
    Write(*,'(a)') '-----> BFIELD METHOD IS XPAND'
    Write(*,'(a,i0)') '-----> xpand_field_eval_type is: ',xpand_field_eval_type
    If (xpand_field_eval_type .eq. 0) Then
      Write(*,'(a)') '-----> Evaluating XPAND fields as PERTURBED!'
      bfield%method = 10
      bfield%method_2d = 0
      bfield%method_pert = -1
    Elseif (xpand_field_eval_type .eq. 1) Then
      Write(*,'(a)') '-----> Evaluating XPAND fields as VACUUM!'
      bfield%method = 11
      bfield%method_2d = 0
      bfield%method_pert = -1
    Else
      Stop "Did not recognize xpand_field_eval_type"
    Endif
    Call readg_g3d(gfile_name,g)
    bfield%g = g
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


Write(*,'(/1a)') "-------------------------------------------------------------"
Write(*,*) 'Setting up bnorm'
Write(*,*) 'Nres = ',nres

Allocate(pnwant(numpn))
pnwant = rlinspace(pnmin,pnmax,numpn)
Allocate(phi(nphi),theta(ntheta))
phi = rlinspace(0.d0,2.d0*pi*Real(nphi-1,real64)/Real(nphi,real64),nphi)
theta = rlinspace(0.d0,2.d0*pi*Real(ntheta-1,real64)/Real(ntheta,real64),ntheta)
dphi = phi(2) - phi(1)
dtheta = theta(2) - theta(1)


Allocate(rpest(numpn,ntheta))
Allocate(zpest(numpn,ntheta))
Allocate(jpest(numpn,ntheta))

! Calculate pest coords for each surface, then Br_mn
Call get_pest_coords(bfield%g,pnwant,ntheta,rpest,zpest,jpest)

Allocate(dpsidr(ntheta),dpsidz(ntheta),psiout(ntheta))
Allocate(rn(ntheta),zn(ntheta))
Allocate(br(ntheta),bz(ntheta),bphi(ntheta))
Allocate(phiarr(ntheta))
Allocate(bnorm(ntheta))
Allocate(area(numpn))
Allocate(alpha_mn(ntheta))
Allocate(br_c(numpn,2*mmax+1),br_s(numpn,2*mmax+1),Br_mn(numpn,2*mmax+1))
Allocate(marr(2*mmax+1))
Do mind = 0,2*mmax
  marr(mind+1) = mind - mmax
Enddo

Call set_bfield_pert_only(bfield)
Do i = 1,numpn
  Call get_psi_derivs_bicub(bfield%g,rpest(i,:),zpest(i,:),ntheta,&
       psiout,dpsidr,dpsidz,ierr)

  rn = dpsidr/Sqrt(dpsidr**2 + dpsidz**2) ! Unit vector in grad(psi)
  zn = dpsidz/Sqrt(dpsidr**2 + dpsidz**2)

  area(i) = Sum(jpest(i,:))*2.d0*pi*dtheta
  
  Do k = 1,nphi
    phiarr(:) = phi(k)
    Call calc_B_rzphi_general(bfield,rpest(i,:),zpest(i,:),phiarr,ntheta,br,bz,bphi)
    bnorm = br*rn + bz*zn
    Do mind = 1,2*mmax + 1
      alpha_mn = nres*phi(k) - marr(mind)*theta
      br_c(i,mind) = br_c(i,mind) + 2.d0*dtheta*dphi*Sum(jpest(i,:)*bnorm*cos(alpha_mn))/area(i) !A.15
      br_s(i,mind) = br_s(i,mind) + 2.d0*dtheta*dphi*Sum(jpest(i,:)*bnorm*sin(alpha_mn))/area(i)
    Enddo
  Enddo
Enddo
Br_mn = Sqrt(br_c**2 + br_s**2) ! After A.15
Call reset_bfield(bfield)

Open(99,file="Brmn.out",status="unknown",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening output file'
  Stop 'Exiting: I/O Error in bnorm_driver.f90'
Endif
Write(99,*) Br_mn
Close(99)

Open(99,file="marr.out",status="unknown",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening output file'
  Stop 'Exiting: I/O Error in bnorm_driver.f90'
Endif
Write(99,*) marr
Close(99)

Open(99,file="pn.out",status="unknown",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening output file'
  Stop 'Exiting: I/O Error in bnorm_driver.f90'
Endif
Write(99,*) pnwant
Close(99)

Deallocate(Br_mn)
Deallocate(marr)
Deallocate(br_c,br_s)
Deallocate(alpha_mn)
Deallocate(dpsidr,dpsidz)
Deallocate(psiout)

Deallocate(phiarr,rn,zn)
Deallocate(br,bz,bphi)
Deallocate(area)
Deallocate(bnorm)


! CLEANUP


Deallocate(rpest,zpest,jpest)
Deallocate(pnwant,phi,theta)
Call Etime(tarray,tres)
Write(*,*) ' bnorm_driver took ',tres-tres0,' seconds'


End Program bnorm_driver





