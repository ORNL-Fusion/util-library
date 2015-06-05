!-----------------------------------------------------------------------------
!
!   Routines/modules related to screening
!   
!-----------------------------------------------------------------------------
Module screening_module
Use kind_mod
Implicit None

!QQ make private

Real(rknd), Dimension(:,:,:), Allocatable :: Ar, Aphi, Az
Real(rknd), Dimension(:,:,:), Allocatable :: Arcoeff, Aphicoeff, Azcoeff
Real(rknd), Dimension(:), Allocatable :: Rvals, Zvals, Phivals
Real(rknd), Dimension(:), Allocatable :: Rknot, Zknot, Phiknot
Real(rknd) :: R_min, R_max, Z_min, Z_max
Integer(iknd) nr, nz, nphi

Integer(iknd), Parameter :: spline_ord = 5_iknd

Contains

!-----------------------------------------------------------------------------
!+ returns Bcyl from B-spline of vector potential
!-----------------------------------------------------------------------------
Subroutine bfield_bspline(rvec,pvec,zvec,Npts,Bout,ierr)
! Description: 
!  
! Output:
!   Bout = (Br,Bz,Bt)
!
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/14/2011  Original Code.  JL
! 
! Author(s): J. Lore 04/14/2011 
!
! Modules used:
Use kind_mod                  ! Import rknd, iknd specifications
Use bspline
Use phys_const, Only: pi
Implicit None

! Input/output                      !See above for descriptions
Integer(iknd),Intent(in) :: Npts
Real(rknd),Dimension(Npts),Intent(in)  :: rvec,pvec,zvec
Real(rknd),Dimension(Npts,3),Intent(out) :: Bout(Npts,3)
Integer(iknd),Intent(out) :: ierr

! Local variables
Integer(iknd) :: ii
Real(rknd) :: rr,zz,pp
Real(rknd) :: dAr_dz,dAr_dp,dAz_dr,dAz_dp,dAp_dr,dAp_dz,Aphi_local
Real(rknd) :: Br,Bz,Bphi

!- End of header -------------------------------------------------------------

If (.not. Allocated(Rknot) ) Then
  Write(*,*) 'B-spline VARIABLES NOT ALLOCATED, EXITING FROM bfield_bspline!'
  Stop
Endif

ierr = 0
Do ii = 1,Npts 
  
  rr = rvec(ii)
  zz = zvec(ii)
  pp = pvec(ii)
  Do While (pp .lt. 0._rknd)
    pp = pp + 2._rknd*pi
  End Do
  pp = Mod(pp,2._rknd*pi)

  ! Check for points off grid
  If ( (rr .lt. R_min) .or. (rr .gt. R_max) ) Then
    Write(*,*) 'bfield_bspline: Point off grid in R: R = ',rr,'. [Rmin,Rmax] = [',R_min,',',R_max,']'
    ierr = 1
    Bout = 0.d0
    return
  Endif
  If ( (zz .lt. Z_min) .or. (zz .gt. Z_max) ) Then
    Write(*,*) 'bfield_bspline: Point off grid in Z: Z = ',zz,'. [Zmin,Zmax] = [',Z_min,',',Z_max,']'
    ierr = 1
    Bout = 0.d0
    return
  Endif


  dAr_dz = dbs3dr(0,1,0,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Arcoeff)
  dAr_dp = dbs3dr(0,0,1,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Arcoeff)
  dAz_dr = dbs3dr(1,0,0,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Azcoeff)
  dAz_dp = dbs3dr(0,0,1,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Azcoeff)
  dAp_dr = dbs3dr(1,0,0,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Aphicoeff)
  dAp_dz = dbs3dr(0,1,0,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Aphicoeff)
  !      Ar     = dbs3dr(0,0,0,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Arcoeff)
  !      Az     = dbs3dr(0,0,0,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Azcoeff)
  Aphi_local   = dbs3dr(0,0,0,rr,zz,pp,spline_ord,spline_ord,spline_ord,Rknot,Zknot,Phiknot,nr,nz,nphi,Aphicoeff)

  Br   = 1.d0/rr * dAz_dp - dAp_dz
  Bphi = dAr_dz - dAz_dr
  Bz   = 1.d0/rr * Aphi_local + dAp_dr - 1.d0/rr * dAr_dp

  Bout(ii,1) = Br
  Bout(ii,2) = Bz
  Bout(ii,3) = Bphi
Enddo

End Subroutine bfield_bspline


!-----------------------------------------------------------------------------
!+ Allocates memory for screening variables
!-----------------------------------------------------------------------------
Subroutine setup_screening_vars(Afile_path)
Use kind_mod
Implicit None
Character(Len=200), Intent(in) :: Afile_path
Character(Len=200) :: fname
Integer(iknd) :: iocheck

Write(*,*) 'Reading A_in2.dat'
fname = Afile_path(1:len_trim(Afile_path))//'A_in2.dat'
Open(99,FILE=fname,STATUS="old",IOSTAT=iocheck)
! Check for success
If ( iocheck /= 0 ) then
  Write(*,*) 'Error opening file: ', fname
  Stop 'Exiting: I/O Error in subroutine read_spline_data'
Endif
Read(99,*) nr
Read(99,*) nz
Read(99,*) nphi
Read(99,*) R_min
Read(99,*) R_max
Read(99,*) Z_min
Read(99,*) Z_max
Close(99)
nphi = nphi + 1  ! 2pi point will be added
Write (*,*) 'Grid resolution: nr=',nr,'nz=',nz,'nphi=', nphi
Write (*,*) 'R range: ',R_min, R_max
Write (*,*) 'Z range: ',Z_min, Z_max
Write (*,*)

Allocate(Ar(nr,nz,nphi))
Allocate(Az(nr,nz,nphi))
Allocate(Aphi(nr,nz,nphi))

Allocate(Rvals(nr))
Allocate(Zvals(nz))
Allocate(phivals(nphi))

! B-Spline Data
Allocate(Rknot(nr+spline_ord), Zknot(nz+spline_ord), Phiknot(nphi+spline_ord))
Allocate(Arcoeff(nr, nz, nphi)) 
Allocate(Azcoeff(nr, nz, nphi))
Allocate(Aphicoeff(nr, nz, nphi))







End Subroutine setup_screening_vars



!-----------------------------------------------------------------------------
!+ Reads B-spline data for vector potential
!-----------------------------------------------------------------------------
Subroutine read_spline_data(spline_data_path)
Use kind_mod
Implicit None
Character(Len=200), Intent(in) :: spline_data_path
Character(Len=200) :: fname
Integer(iknd) :: i,j,k, iocheck

Real(rknd) :: R_min_tmp, R_max_tmp, Z_min_tmp, Z_max_tmp
Integer(iknd) :: spline_ord_tmp, nr_tmp, nz_tmp, nphi_tmp

fname = spline_data_path(1:len_trim(spline_data_path))//'spline_data.dat'
Write(*,*) 'Reading B-spline coefficients from file:',fname
Open(99,FILE=fname,STATUS="old",IOSTAT=iocheck)
! Check for success
If ( iocheck /= 0 ) then
  Write(*,*) 'Error opening file: ', fname
  Stop 'Exiting: I/O Error in subroutine read_spline_data'
Endif
Read(99,*) R_min_tmp, R_max_tmp, Z_min_tmp, Z_max_tmp
Read(99,*) nr_tmp, nz_tmp, nphi_tmp
Read(99,*) spline_ord_tmp
Read(99,*) Rknot
Read(99,*) Zknot
Read(99,*) Phiknot
Write(*,*) 'Reading Arcoeff'
Do i = 1,nr
  Do j = 1,nz
    Do k = 1,nphi
      Read(99,*) Arcoeff(i,j,k)
    Enddo
  Enddo
Enddo
Write(*,*) 'Reading Azcoeff'
Do i = 1,nr
  Do j = 1,nz
    Do k = 1,nphi
      Read(99,*) Azcoeff(i,j,k)
    Enddo
  Enddo
Enddo
Write(*,*) 'Read Aphicoeff'
Do i = 1,nr
  Do j = 1,nz
    Do k = 1,nphi
      Read(99,*) Aphicoeff(i,j,k)
    Enddo
  Enddo
Enddo

Close(99)

Write(*,*) '---> Done reading B-spline coefficients'



End Subroutine read_spline_data


!-----------------------------------------------------------------------------
!+ Reads Ar,Az,Afi files
!-----------------------------------------------------------------------------
Subroutine read_Afiles(Afile_path)
Use kind_mod
Use math_geo_module
Use phys_const, Only: pi
Implicit None
Character(Len=200), Intent(in) :: Afile_path
Character(Len=200) :: fname
Integer(iknd) :: i,j,k

Write(*,*) 'Reading Afile data from:',Afile_path



Rvals = rlinspace(R_min,R_max,nr)
Zvals = rlinspace(Z_min,Z_max,nz)
phivals = rlinspace(0._rknd,2._rknd*pi,nphi)

Write(*,*) 'Reading Ar.dat'
fname = Afile_path(1:len_trim(Afile_path))//'Ar.dat'
Open(99,file=fname)
Do i = 1,nr
  Do j = 1,nz
    Do k = 1,nphi - 1 
      Read(99,*) Ar(i,j,k)
    Enddo
  Enddo
Enddo
Close(99)

Write(*,*) 'Reading Az.dat'
fname = Afile_path(1:len_trim(Afile_path))//'Az.dat'
Open(99,file=fname)
Do i = 1,nr
  Do j = 1,nz
    Do k = 1,nphi - 1
      Read(99,*) Az(i,j,k)
    Enddo
  Enddo
Enddo
Close(99)


Write(*,*) 'Reading Afi.dat'
fname = Afile_path(1:len_trim(Afile_path))//'Afi.dat'
Open(99,file=fname)
Do i = 1,nr
  Do j = 1,nz
    Do k = 1,nphi - 1
      Read(99,*) Aphi(i,j,k)
    Enddo
  Enddo
Enddo
Close(99)

! Make periodic
Write(*,*) ' Adding 2pi point to A arrays'
Ar(:,:,nphi) = Ar(:,:,1)
Az(:,:,nphi) = Az(:,:,1)
Aphi(:,:,nphi) = Aphi(:,:,1)
End Subroutine read_Afiles



!-----------------------------------------------------------------------------
!+ Performs B-spline fits to vector potentials
!  -- Note that Ar, Az, Aphi are deallocated after this call!!
!-----------------------------------------------------------------------------
Subroutine prepare_Afile_splinefits !(Afile_path)

Use kind_mod
Use bspline
Implicit None
!Character(Len=200), Intent(in) :: Afile_path
!Character(Len=200) :: fname
!Integer(iknd) :: i,j,k


Write(*,*) 'Preparing spline fits to Afile data'

Call dbsnak(nr,Rvals,spline_ord,Rknot)
Call dbsnak(nz,Zvals,spline_ord,Zknot)
Call dbsnak(nphi,phivals,spline_ord,Phiknot)

Write(*,*) 'Preparing B-spline coefficients for Ar'
Call dbs3in(nr,Rvals,nz,Zvals,nphi,Phivals,Ar,nr,nz,spline_ord,spline_ord,spline_ord, &
     Rknot,Zknot,Phiknot,Arcoeff)
Deallocate(Ar)
Write(*,*) 'Preparing B-spline coefficients for Az'
Call dbs3in(nr,Rvals,nz,Zvals,nphi,Phivals,Az,nr,nz,spline_ord,spline_ord,spline_ord, &
     Rknot,Zknot,Phiknot,Azcoeff)
Deallocate(Az)
Write(*,*) 'Preparing B-spline coefficients for Aphi'
Call dbs3in(nr,Rvals,nz,Zvals,nphi,Phivals,Aphi,nr,nz,spline_ord,spline_ord,spline_ord, &
     Rknot,Zknot,Phiknot,Aphicoeff)
Deallocate(Aphi)

!!$Write(*,*) 'Saving B-spline coefficients'
!!$fname = Afile_path(1:len_trim(Afile_path))//'spline_data.dat'
!!$Open(99,file=fname)
!!$Write(99,*) R_min, R_max, Z_min, Z_max
!!$Write(99,*) nr, nz, nphi
!!$Write(99,*) spline_ord
!!$Write(99,*) Rknot
!!$Write(99,*) Zknot
!!$Write(99,*) Phiknot
!!$Write(*,*) 'Saving Arcoeff'
!!$Do i = 1,nr
!!$  Do j = 1,nz
!!$    Do k = 1,nphi
!!$      Write(99,*) Arcoeff(i,j,k)
!!$    Enddo
!!$  Enddo
!!$Enddo
!!$Write(*,*) 'Saving Azcoeff'
!!$Do i = 1,nr
!!$  Do j = 1,nz
!!$    Do k = 1,nphi
!!$      Write(99,*) Azcoeff(i,j,k)
!!$    Enddo
!!$  Enddo
!!$Enddo
!!$Write(*,*) 'Saving Aphicoeff'
!!$Do i = 1,nr
!!$  Do j = 1,nz
!!$    Do k = 1,nphi
!!$      Write(99,*) Aphicoeff(i,j,k)
!!$    Enddo
!!$  Enddo
!!$Enddo
!!$Close(99)


Write(*,*) 'Note: Ar, Az, and Apfi have been deallocated!'
End Subroutine prepare_Afile_splinefits


End Module screening_module
