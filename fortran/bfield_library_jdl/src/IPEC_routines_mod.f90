!-----------------------------------------------------------------------------
!
!   Routines/modules related to IPEC
!
!
!   Contains:
!     Subroutine read_ipec_field_file
!     Subroutine bfield_ipec
!     Subroutine open_ipec_fields
!     Subroutine close_ipec_fields
!-----------------------------------------------------------------------------
Module ipec_module
Use kind_mod, Only: real64, int32
Implicit None
!Integer(int32), Save :: ipec_field_eval_type = -1

Real(real64), Allocatable, Save :: ipec_r(:), ipec_z(:)
Real(real64), Allocatable, Dimension(:,:), Save :: &
     ipec_eq_br, ipec_eq_bz, ipec_eq_bphi, &
     ipec_vac_rbr, ipec_vac_ibr, ipec_vac_rbz, ipec_vac_ibz, ipec_vac_rbphi, ipec_vac_ibphi, &
     ipec_pert_rbr, ipec_pert_ibr, ipec_pert_rbz, ipec_pert_ibz, ipec_pert_rbphi, ipec_pert_ibphi     
Integer(int32), Save :: ipec_nr, ipec_nz, ipec_vac_n, ipec_pert_n

Contains

!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine read_ipec_field_file(fname,ifield)
!
! Description:
! Input:
!
! Output:
!
! Calls:
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!   
! Author(s): J.D. Lore 1/27/2014
Use kind_mod, Only: int32
Implicit None
Character(len=*), Intent(In) :: fname
Integer(int32), Intent(In) :: ifield
Integer(int32) :: iocheck, n,ind,nr,nz,field_size
Character(len=100) :: sjunk
Real(real64), Allocatable :: data(:,:)
!- End of header -------------------------------------------------------------

Write(*,'(a,a)') ' Reading ', Trim(fname)
open(99,file=fname,IOSTAT=iocheck,STATUS="old")
If (iocheck /= 0) Then
  Write(*,*) 'Error opening file: ',fname
  Stop 'Exiting: I/O error in function read_ipec_field_file'
Endif

Read(99,'(A)') sjunk  ! name  
Read(99,'(A)') sjunk  !version
Read(99,'(A)') sjunk  !blank
Read(99,'(A)') sjunk  ! n =
ind = Scan(sjunk,"=")
!write(*,*) 'THIS: |',sjunk(ind+1:ind+6),"|"
read(sjunk(ind+1:ind+6),*) n
write(*,*) 'n = ',n
Read(99,'(A)') sjunk  ! nr, nz
ind = Scan(sjunk,"=",.true.)
read(sjunk(ind+1:ind+6),*) nr
ind = Scan(sjunk,"=")
read(sjunk(ind+1:ind+6),*) nz
write(*,*) 'nr = ',nr,',nz = ',nz
Read(99,'(A)') sjunk  !blank
Read(99,'(A)') sjunk  !label

field_size = 6
Select Case (ifield)
Case (0)

  ipec_nr = nr
  ipec_nz = nz
  
  field_size = 3
  Allocate(data(3+field_size,nr*nz))
  Read(99,*) data
  Allocate(ipec_r(nr))
  Allocate(ipec_z(nz))
  Allocate(ipec_eq_br(nr,nz))
  Allocate(ipec_eq_bz(nr,nz))
  Allocate(ipec_eq_bphi(nr,nz))

  ipec_r = data(2,1:nr*nz:nz)
  ipec_z = data(3,1:nz)

  ipec_eq_br   = Transpose(Reshape(data(4,1:nr*nz),(/nr,nz/)))
  ipec_eq_bz   = Transpose(Reshape(data(5,1:nr*nz),(/nr,nz/)))
  ipec_eq_bphi = Transpose(Reshape(data(6,1:nr*nz),(/nr,nz/)))
Case (1)
  ipec_vac_n = n
  field_size = 6
  Allocate(data(3+field_size,nr*nz))
  Read(99,*) data
  Allocate(ipec_vac_rbr(nr,nz))
  Allocate(ipec_vac_rbz(nr,nz))
  Allocate(ipec_vac_rbphi(nr,nz))
  Allocate(ipec_vac_ibr(nr,nz))
  Allocate(ipec_vac_ibz(nr,nz))
  Allocate(ipec_vac_ibphi(nr,nz))

  ipec_vac_rbr   = Transpose(Reshape(data(4,1:nr*nz),(/nr,nz/)))
  ipec_vac_ibr   = Transpose(Reshape(data(5,1:nr*nz),(/nr,nz/)))
  ipec_vac_rbz   = Transpose(Reshape(data(6,1:nr*nz),(/nr,nz/)))
  ipec_vac_ibz   = Transpose(Reshape(data(7,1:nr*nz),(/nr,nz/)))
  ipec_vac_rbphi = Transpose(Reshape(data(8,1:nr*nz),(/nr,nz/)))
  ipec_vac_ibphi = Transpose(Reshape(data(9,1:nr*nz),(/nr,nz/)))
Case (2)
  ipec_pert_n = n
  field_size = 6
  Allocate(data(3+field_size,nr*nz))
  Read(99,*) data
  Allocate(ipec_pert_rbr(nr,nz))
  Allocate(ipec_pert_rbz(nr,nz))
  Allocate(ipec_pert_rbphi(nr,nz))
  Allocate(ipec_pert_ibr(nr,nz))
  Allocate(ipec_pert_ibz(nr,nz))
  Allocate(ipec_pert_ibphi(nr,nz))

  ipec_pert_rbr   = Transpose(Reshape(data(4,1:nr*nz),(/nr,nz/)))
  ipec_pert_ibr   = Transpose(Reshape(data(5,1:nr*nz),(/nr,nz/)))
  ipec_pert_rbz   = Transpose(Reshape(data(6,1:nr*nz),(/nr,nz/)))
  ipec_pert_ibz   = Transpose(Reshape(data(7,1:nr*nz),(/nr,nz/)))
  ipec_pert_rbphi = Transpose(Reshape(data(8,1:nr*nz),(/nr,nz/)))
  ipec_pert_ibphi = Transpose(Reshape(data(9,1:nr*nz),(/nr,nz/)))  
Case Default
  Stop "Did not recognize case in read_ipec_field_file"
End Select

Deallocate(data)
Close(99)


End Subroutine read_ipec_field_file

Subroutine open_ipec_fields(run_path,itype)
  Use kind_mod, Only: int32
  Implicit None
  Character(Len=*), Intent(In) :: run_path
  Integer(int32), Intent(In) :: itype ! 1 = ipec, 2 = gpec

  Select Case (itype)
    Case (1) 
      Call read_ipec_field_file(Trim(run_path)//'/ipec_eqbrzphi_n3.out',0)
      Call read_ipec_field_file(Trim(run_path)//'/ipec_cbrzphi_n3.out',1)
      Call read_ipec_field_file(Trim(run_path)//'/ipec_brzphi_n3.out',2)
    Case (2)
      Call read_ipec_field_file(Trim(run_path)//'/gpec_eqbrzphi_n3.out',0)
      Call read_ipec_field_file(Trim(run_path)//'/gpec_cbrzphi_n3.out',1)
      Call read_ipec_field_file(Trim(run_path)//'/gpec_brzphi_n3.out',2)      
    Case Default
      Write(*,*) "Did not recognize itype in open_ipec_fields",itype
      Stop "quitting from open_ipec_fields"
    End Select

End Subroutine open_ipec_fields

Subroutine close_ipec_fields
  Implicit None
  Deallocate(ipec_r,ipec_z)
  Deallocate(ipec_eq_br,ipec_eq_bz,ipec_eq_bphi)
  Deallocate(ipec_vac_rbr,ipec_vac_ibr,ipec_vac_rbz,ipec_vac_ibz,ipec_vac_ibphi,ipec_vac_rbphi)
  Deallocate(ipec_pert_rbr,ipec_pert_ibr,ipec_pert_rbz,ipec_pert_ibz,ipec_pert_ibphi,ipec_pert_rbphi)  
End Subroutine close_ipec_fields

!-----------------------------------------------------------------------------
!+ Evaluate B(r,z) using Equilibrium only IPEC fields
!-----------------------------------------------------------------------------
Subroutine bfield_ipec(r,phi,z,Npts,Bout,ierr,ifield_type)
  ! ifield_type : 0=eq, 1=eq+vac, 2=eq+pert, 3=vac_only, 4=pert_only
  ! Output:
  !   Bout = (:,[Br,Bz,Bt])
Use kind_mod, Only: int32, real64
Implicit None
Real(Real64), Intent(In), Dimension(Npts) :: r, z, phi
Integer(int32), Intent(In) :: Npts, ifield_type
Real(real64), Intent(Out), Dimension(Npts,3) :: Bout
Integer(int32), Intent(Out) :: ierr
! Local variables
Real(real64) :: dz_grid, dr_grid, dr1, dr2, dz1, dz2, QQ(2,2)
Real(real64) :: breal, bimg, cosphi, sinphi
Integer(int32) :: i, ir, iz
ierr = 0
Bout(:,:) = 0._real64
Do i=1,Npts
  If (r(i) .lt. ipec_r(1) .OR. r(i) .gt. ipec_r(ipec_nr-1) &
       .OR. z(i) .lt. ipec_z(1) .OR. z(i) .gt. ipec_z(ipec_nz-1)) Then
    ierr = 1
    Bout(i,1:3) = 0._real64
    Cycle
  Endif
  ir = Floor((r(i) - ipec_r(1))/(ipec_r(2)-ipec_r(1))) + 1
  iz = Floor((z(i) - ipec_z(1))/(ipec_z(2)-ipec_z(1))) + 1
  
  dr_grid = ipec_r(ir+1) - ipec_r(ir)
  dz_grid = ipec_z(iz+1) - ipec_z(iz)
  
  dr2 = ipec_r(ir+1) - r(i)
  dr1 = dr_grid - dr2
  dz2 = ipec_z(iz+1) - z(i)
  dz1 = dz_grid - dz2

  ! 2D COMPONENT -- skip for perturbation part only
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  If (ifield_type .lt. 3) Then
    QQ = ipec_eq_br(ir:ir+1,iz:iz+1)
    Bout(i,1) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    QQ = ipec_eq_bz(ir:ir+1,iz:iz+1)
    Bout(i,2) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    QQ = ipec_eq_bphi(ir:ir+1,iz:iz+1)
    Bout(i,3) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
  Endif

  
  ! 3D COMPONENT -- skip for AS only
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  If (ifield_type .eq. 0) Then
    Cycle
  Elseif (ifield_type .eq. 1 .OR. ifield_type .eq. 3) Then
    cosphi = cos(ipec_vac_n*phi(i))
    sinphi = sin(ipec_vac_n*phi(i))
    
    QQ = ipec_vac_rbr(ir:ir+1,iz:iz+1)
    breal = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    QQ = ipec_vac_ibr(ir:ir+1,iz:iz+1)
    bimg = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)        
    Bout(i,1) = Bout(i,1) + breal*cosphi + bimg*sinphi

    QQ = ipec_vac_rbz(ir:ir+1,iz:iz+1)
    breal = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    QQ = ipec_vac_ibz(ir:ir+1,iz:iz+1)
    bimg = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    Bout(i,2) = Bout(i,2) + breal*cosphi + bimg*sinphi
        
    QQ = ipec_vac_rbphi(ir:ir+1,iz:iz+1)
    breal = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    QQ = ipec_vac_ibphi(ir:ir+1,iz:iz+1)
    bimg = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    Bout(i,3) = Bout(i,3) + breal*cosphi + bimg*sinphi
  Elseif (ifield_type .eq. 2 .OR. ifield_type .eq. 4) Then
    cosphi = cos(ipec_pert_n*phi(i))
    sinphi = sin(ipec_pert_n*phi(i))
    
    QQ = ipec_pert_rbr(ir:ir+1,iz:iz+1)
    breal = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    QQ = ipec_pert_ibr(ir:ir+1,iz:iz+1)
    bimg = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)        
    Bout(i,1) = Bout(i,1) + breal*cosphi + bimg*sinphi

    QQ = ipec_pert_rbz(ir:ir+1,iz:iz+1)
    breal = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    QQ = ipec_pert_ibz(ir:ir+1,iz:iz+1)
    bimg = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    Bout(i,2) = Bout(i,2) + breal*cosphi + bimg*sinphi
        
    QQ = ipec_pert_rbphi(ir:ir+1,iz:iz+1)
    breal = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    QQ = ipec_pert_ibphi(ir:ir+1,iz:iz+1)
    bimg = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid)
    Bout(i,3) = Bout(i,3) + breal*cosphi + bimg*sinphi    
  Else
    Write(*,*) "Did not recognize ifield_type in bfield_ipec ",ifield_type
    Stop "Exiting"
  Endif
Enddo

End Subroutine bfield_ipec

End Module ipec_module
