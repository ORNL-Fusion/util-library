!-----------------------------------------------------------------------------
!
!   Routines/modules related to XPAND
!
!
!   Contains:
!-----------------------------------------------------------------------------
Module xpand_module
Use kind_mod, Only: real64, int32
Implicit None

Real(real64), Allocatable, Save :: xpand_r(:), xpand_z(:), xpand_phi(:)
Real(real64), Allocatable, Dimension(:,:,:), Save :: &
     xpand_vac_br, xpand_vac_bz, xpand_vac_bphi,&
     xpand_pert_br,xpand_pert_bz, xpand_pert_bphi
Integer(int32), Save :: xpand_nr, xpand_nz, xpand_nphi

Contains

!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine read_xpand_field_file(fname)
! Author(s): J.D. Lore
  Use kind_mod, Only: int32
  Use phys_const, Only: pi
Implicit None
Character(len=*), Intent(In) :: fname
Integer(int32) :: iocheck, nphi,ind,nr,nz,field_size
Character(len=256) :: sjunk
Real(real64), Allocatable :: data(:,:)
!- End of header -------------------------------------------------------------

Write(*,'(a,a)') ' Reading ', Trim(fname)
open(99,file=fname,IOSTAT=iocheck,status='old')
If (iocheck /= 0) Then
  Write(*,*) 'Error opening file: ',fname
  Stop 'Exiting: I/O error in function read_xpand_field_file'
Endif

Read(99,'(A)') sjunk  ! label
Write(*,*) sjunk
Read(99,'(A)') sjunk  ! NR =
ind = Scan(sjunk,"=")
!write(*,*) 'THIS: |',sjunk(ind+1:ind+6),"|"
read(sjunk(ind+1:ind+6),*) nr
Read(99,'(A)') sjunk  ! Np =
ind = Scan(sjunk,"=")
read(sjunk(ind+1:ind+6),*) nphi
Read(99,'(A)') sjunk  ! NZ =
ind = Scan(sjunk,"=")
read(sjunk(ind+1:ind+6),*) nz

write(*,*) 'nr   = ',nr
write(*,*) 'nphi = ',nphi
write(*,*) 'nz   = ',nz
Read(99,'(A)') sjunk  !label


field_size = 10
Allocate(data(10,nr*nz*nphi))
Read(99,*) data
Close(99)

nphi = nphi + 1 ! Will add first slice to 2*pi

Allocate(xpand_r(nr))
Allocate(xpand_z(nz))
Allocate(xpand_phi(nphi))
Allocate(xpand_pert_br(nr,nz,nphi))
Allocate(xpand_pert_bz(nr,nz,nphi))
Allocate(xpand_pert_bphi(nr,nz,nphi))
Allocate(xpand_vac_br(nr,nz,nphi))
Allocate(xpand_vac_bz(nr,nz,nphi))
Allocate(xpand_vac_bphi(nr,nz,nphi))

xpand_r = data(1,1:nr)
xpand_z = data(3,1:nz*nr:nr)
xpand_phi(1:nphi-1) = data(2,1:nz*nr*(nphi-1):nr*nz)
xpand_phi(nphi) = 2._real64*pi

xpand_pert_br(1:nr,1:nz,1:nphi-1)     = Reshape(data(4,1:nr*nz*(nphi-1)),(/nr,nz,nphi-1/))
xpand_pert_br(:,:,nphi) = xpand_pert_br(:,:,1)
xpand_pert_bphi(1:nr,1:nz,1:nphi-1)   = Reshape(data(5,1:nr*nz*(nphi-1)),(/nr,nz,nphi-1/))
xpand_pert_bphi(:,:,nphi) = xpand_pert_bphi(:,:,1)
xpand_pert_bz(1:nr,1:nz,1:nphi-1)     = Reshape(data(6,1:nr*nz*(nphi-1)),(/nr,nz,nphi-1/))
xpand_pert_bz(:,:,nphi) = xpand_pert_bz(:,:,1)

xpand_vac_br(1:nr,1:nz,1:nphi-1)     = Reshape(data(8,1:nr*nz*(nphi-1)),(/nr,nz,nphi-1/))
xpand_vac_br(:,:,nphi) = xpand_vac_br(:,:,1)
xpand_vac_bphi(1:nr,1:nz,1:nphi-1)   = Reshape(data(9,1:nr*nz*(nphi-1)),(/nr,nz,nphi-1/))
xpand_vac_bphi(:,:,nphi) = xpand_vac_bphi(:,:,1)
xpand_vac_bz(1:nr,1:nz,1:nphi-1)     = Reshape(data(10,1:nr*nz*(nphi-1)),(/nr,nz,nphi-1/))
xpand_vac_bz(:,:,nphi) = xpand_vac_bz(:,:,1)

Deallocate(data)

xpand_nr = nr
xpand_nz = nz
xpand_nphi = nphi

End Subroutine read_xpand_field_file

Subroutine open_xpand_fields(fname)
  Implicit None
  Character(Len=*), Intent(In) :: fname
  Call read_xpand_field_file(fname)
End Subroutine open_xpand_fields

Subroutine close_xpand_fields
  Implicit None
  Deallocate(xpand_r,xpand_z,xpand_phi)
  Deallocate(xpand_pert_br,xpand_pert_bz,xpand_pert_bphi)
  Deallocate(xpand_vac_br,xpand_vac_bz,xpand_vac_bphi)
End Subroutine close_xpand_fields

!-----------------------------------------------------------------------------
!+ Evaluate B(r,z) using Equilibrium only XPAND fields
!-----------------------------------------------------------------------------
Subroutine bfield_xpand(r,phi,z,Npts,Bout,ierr,ifield_type)
! ifield_type : 0=pert, 1=vac
  ! Output:
!   Bout = (:,[Br,Bz,Bt])
  Use kind_mod, Only: int32, real64
  Use phys_const, Only: pi
Implicit None
Real(Real64), Intent(In), Dimension(Npts) :: r, z, phi
Integer(int32), Intent(In) :: Npts, ifield_type
Real(real64), Intent(Out), Dimension(Npts,3) :: Bout
Integer(int32), Intent(Out) :: ierr
! Local variables
Real(real64) :: phi_tmp, phi_fac, dphi_grid
Real(real64) :: dz_grid, dr_grid, dr1, dr2, dz1, dz2, QQ1(2,2), QQ2(2,2)
Integer(int32) :: i, ir, iz, iphi
ierr = 0

Do i=1,Npts

  Bout(i,1:3) = 0._real64
  If (r(i) .lt. xpand_r(1) .OR. r(i) .gt. xpand_r(xpand_nr-1) &
       .OR. z(i) .lt. xpand_z(1) .OR. z(i) .gt. xpand_z(xpand_nz-1)) Then
    ierr = 1
    Cycle
  Endif
  
  ir = Floor((r(i) - xpand_r(1))/(xpand_r(2)-xpand_r(1))) + 1
  iz = Floor((z(i) - xpand_z(1))/(xpand_z(2)-xpand_z(1))) + 1
  phi_tmp = phi(i)
  Do While (phi_tmp .lt. 0._real64)
    phi_tmp = phi_tmp + 2._real64*pi
  Enddo
  phi_tmp = Mod(phi_tmp,2._real64*pi)
  iphi = Floor((phi_tmp - xpand_phi(1))/(xpand_phi(2)-xpand_phi(1))) + 1
  If (iphi .lt. 1 .OR. iphi .gt. xpand_nphi - 1) Then
    Write(*,*) 'iphi out of range... should not happen',iphi
    Write(*,*) (iphi .lt. 1)
    Write(*,*) xpand_nphi -1
    Stop "giving up"
  Endif
  
  dr_grid = xpand_r(ir+1) - xpand_r(ir)
  dz_grid = xpand_z(iz+1) - xpand_z(iz)
  dphi_grid = xpand_phi(iphi+1) - xpand_phi(iphi)

  phi_fac = (phi_tmp - xpand_phi(iphi))/dphi_grid
  
  dr2 = xpand_r(ir+1) - r(i)
  dr1 = dr_grid - dr2
  dz2 = xpand_z(iz+1) - z(i)
  dz1 = dz_grid - dz2

  If (ifield_type .eq. 0) Then
    QQ1 = xpand_pert_br(ir:ir+1,iz:iz+1,iphi)
    QQ2 = xpand_pert_br(ir:ir+1,iz:iz+1,iphi+1)
    Bout(i,1) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
         phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
    QQ1 = xpand_pert_bz(ir:ir+1,iz:iz+1,iphi)
    QQ2 = xpand_pert_bz(ir:ir+1,iz:iz+1,iphi+1)
    Bout(i,2) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
         phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
    QQ1 = xpand_pert_bphi(ir:ir+1,iz:iz+1,iphi)
    QQ2 = xpand_pert_bphi(ir:ir+1,iz:iz+1,iphi+1)
    Bout(i,3) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
         phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
  Elseif (ifield_type .eq. 1) Then
    Write(*,*) 'This only works with gfile!'
    Stop
    QQ1 = xpand_vac_br(ir:ir+1,iz:iz+1,iphi)
    QQ2 = xpand_vac_br(ir:ir+1,iz:iz+1,iphi+1)
    Bout(i,1) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
         phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
    QQ1 = xpand_vac_bz(ir:ir+1,iz:iz+1,iphi)
    QQ2 = xpand_vac_bz(ir:ir+1,iz:iz+1,iphi+1)
    Bout(i,2) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
         phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
    QQ1 = xpand_vac_bphi(ir:ir+1,iz:iz+1,iphi)
    QQ2 = xpand_vac_bphi(ir:ir+1,iz:iz+1,iphi+1)
    Bout(i,3) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
         phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
  Else
    Stop "Did not recognize ifield_type in bfield_xpand"
  Endif
Enddo

End Subroutine bfield_xpand

End Module xpand_module
