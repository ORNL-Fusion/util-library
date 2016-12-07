!-----------------------------------------------------------------------------
!+ Program to make bnorm plots
!-----------------------------------------------------------------------------
Program bnorm_driver
!
! Author(s): J.D. Lore - 10/13/2016
!
! Modules used:
Use kind_mod, Only: int32, real64
Use g3d_module, Only : get_psi_derivs_bicub
Use bnorm_routines, Only : get_pest_coords
Use math_geo_module, Only : rlinspace
Use phys_const, Only: pi
Use bfield, Only : calc_b_rzphi_general, set_bfield_pert_only, reset_bfield
Use bnorm_bfield_mod
Implicit none


Integer(int32) :: i, k, mind
Integer(int32) :: ierr,iocheck
Real(real64) :: tstart, tend

Integer(int32), Allocatable :: marr(:)
Real(real64) :: dphi, dtheta
Real(real64), Allocatable :: pnwant(:), phi(:), theta(:)
Real(real64), Allocatable :: rpest(:,:), zpest(:,:), jpest(:,:),dpsidr(:),dpsidz(:),psiout(:)
Real(real64), Allocatable :: rn(:),zn(:),br(:),bz(:),bphi(:),bnorm(:), area(:),phiarr(:),alpha_mn(:)
Real(real64), Allocatable :: br_c(:,:), br_s(:,:), Br_mn(:,:)
!---------------------------------------------------------------------------
! Namelist variables:
Real(real64) :: &
 pnmin = 0.d0, &
 pnmax = 0.d0 

Integer(int32) :: &
 nres = 0, &
 ntheta = 0, &
 nphi = 0, &
 mmax = 0, &
 numpn = 0

! Namelist files
Namelist / settings_nml / &
     nres, ntheta, nphi, mmax, pnmin, pnmax, numpn
     
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

Call cpu_time(tstart)
Write(*,'(/1a)') "-------------------------------------------------------------"
Write(*,'(a)') " Starting bnorm driver"

Open(99,file="bnorm_settings.nml",status="old",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening namelist file'
  Stop 'Exiting: I/O Error in bnorm_driver.f90'
Endif
Read(99,nml=settings_nml)
Close(99)

Write(*,'(/1a)') "-------------------------------------------------------------"
Write(*,*) 'Setting up bnorm'
Write(*,*) 'Nres = ',nres

Call setup_bfield_bnorm

Allocate(pnwant(numpn))
pnwant = rlinspace(pnmin,pnmax,numpn)

Allocate(rpest(numpn,ntheta))
Allocate(zpest(numpn,ntheta))
Allocate(jpest(numpn,ntheta))

! Calculate pest coords for each surface, then Br_mn
Write(*,*) 'Calculating PEST coordinates'
Call get_pest_coords(bfield%g,pnwant,ntheta,rpest,zpest,jpest)
Call cpu_time(tend)
Write(*,*) 'Calculating PEST coordinates took ',tend-tstart,' seconds'

Allocate(phi(nphi),theta(ntheta))
phi = rlinspace(0.d0,2.d0*pi*Real(nphi-1,real64)/Real(nphi,real64),nphi)
theta = rlinspace(0.d0,2.d0*pi*Real(ntheta-1,real64)/Real(ntheta,real64),ntheta)
dphi = phi(2) - phi(1)
dtheta = theta(2) - theta(1)

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

br_c = 0.d0
br_s = 0.d0
Call set_bfield_pert_only(bfield)
Do i = 1,numpn
  Write(*,*) 'Working on surf ',i,' of ',numpn
  Call get_psi_derivs_bicub(bfield%g,rpest(i,:),zpest(i,:),ntheta,psiout,dpsidr,dpsidz,ierr)

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
Call reset_bfield(bfield)

Br_mn = Sqrt(br_c**2 + br_s**2) ! After A.15

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
Call cpu_time(tend)
Write(*,*) ' bnorm_driver took ',tend-tstart,' seconds'


End Program bnorm_driver





