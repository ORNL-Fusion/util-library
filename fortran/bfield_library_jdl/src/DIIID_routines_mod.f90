!-----------------------------------------------------------------------------
!
!   Routines/modules related to DIII-D
!
!   Contains:
!     Subroutine build_d3d_ccoils_jl
!     Subroutine build_d3d_icoils_jl
!   
!-----------------------------------------------------------------------------
Module DIIID_routines_mod
  Use coil_typedef, Only : coil_type, allocate_coiltype
  Implicit None
  Private
  Public :: build_d3d_ccoils_jl
  Public :: build_d3d_icoils_jl
Contains

!-----------------------------------------------------------------------------
!+ Builds the d3d c-coil set as series of current filaments
!-----------------------------------------------------------------------------
Subroutine build_d3d_ccoils_jl(coil,taper,ntorpts)
!
! Set up so TAPER can be taken directly from curntic in diiidsup.in
!
! Description:
!  Builds the coil geometry and sets the current (in Amp-turns) for 
!  the 6 D3D C coils. A single loop is returned for each entire 
!  coil set, with one turn per coil. Output can be used in bfield_bs_jdl, bfield_bs_cyl, etc.
!
!  Coils 1:6 --> C79, C139, C199, C259, C319, C19
!    Coil toroidal angles are in RH cylindrical coordinates, so negative
!    of machine geographical toroidal angle.  <-- Following M. Schaffer routines
!
! Input:
!  taper : real64(6) : The coil currents in Amps.  Number of turns (4) is applied in the routine
!  ntorpts : int32 : number of toroidal points used to generate the coils.  (4 was used in Shaffer's routine)
!
! Output:
!  coil : real64(ncoil_pts,3) : Array of coil filament beginning and end points,
!                             in Cartesian coordinates (m), second index is (1,2,3) = x,y,z
!  current : real64(ncoil_pts) : The current (in Amp-turns) corresponding to the current filaments. 
!                              The filament defined from coil(i,:) to coil(i+1,:) has a current
!                              set by current(i). The last value of array current is not used.  
!  ncoil_pts : int32 : number of coil points
!
! Calls:
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!   
! Author(s): J.D. Lore
!  -- Built on code by J.M. Canik and M. Schaffer
Use kind_mod, Only: real64, int32
Use phys_const, Only: pi  
Implicit None

Integer, Parameter :: nc = 6
Real(real64), Intent(in) :: taper(nc)
Integer(int32), Intent(in), Optional :: ntorpts
Type(coil_type), Intent(Out) :: coil


Integer(int32) :: npts,nturn,i,j,myntorpts
Real(real64) :: phicens(nc),phiext,R(2),Z(2),taper2(nc)
Real(real64), Allocatable :: phi(:)
!- End of header -------------------------------------------------------------

If (Present(ntorpts)) Then
  myntorpts = ntorpts
Else
  myntorpts = 4
Endif

Call allocate_coiltype(nc,myntorpts,coil)

! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------
!             COIL DATA
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------
! Coils are approximately rectangular.
phicens = (/-79.d0,-139.d0,-199.d0,-259.d0,-319.d0,-19.d0/)  ! Center of each coil
phicens = phicens*pi/180.d0
phiext = 58.d0*pi/180.d0               ! Toroidal extent of each coil
R = (/3.23d0,3.23d0/)                  ! Major radius of coil
Z = (/0.8d0,-0.8d0/)                   ! Vertical extent of coil.
nturn = 4
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------

taper2 = real(nturn,real64)*taper

Allocate(phi(myntorpts))
npts = coil%ncoil_pts/coil%num_coils

Do i = 0,5
  Do j = 1,myntorpts
    phi(j) = phiext*real(j-1,real64)/real(myntorpts-1,real64) + phicens(i+1) - phiext/2.d0
  Enddo
  coil%coilxyz(i*npts+1:i*npts+myntorpts,1) = R(1)*cos(phi)
  coil%coilxyz(i*npts+1:i*npts+myntorpts,2) = R(1)*sin(phi)
  coil%coilxyz(i*npts+1:i*npts+myntorpts,3) = Z(1)
  coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,1) = R(2)*cos(phi(myntorpts:1:-1))
  coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,2) = R(2)*sin(phi(myntorpts:1:-1))
  coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,3) = Z(2)
  coil%coilxyz((i+1)*npts,:) = coil%coilxyz(i*npts+1,:)
  coil%current(i*npts+1:(i+1)*npts) = taper2(i+1)
  coil%current((i+1)*npts) = 0.0d0   ! Sticks connecting the coils have no current.
Enddo
Deallocate(phi)
End Subroutine build_d3d_ccoils_jl


!-----------------------------------------------------------------------------
!+ Builds the d3d i-coil set as series of current filaments
!-----------------------------------------------------------------------------
Subroutine build_d3d_icoils_jl(coil,taper,ntorpts)
!
! Set up so TAPER can be taken directly from curntic in diiidsup.in
!
! Description:
!  Builds the coil geometry and sets the current (in Amp-turns) for 
!  the 12 D3D I coils. A single loop is returned for each entire 
!  coil set, with one turn per coil. Output can be used in bfield_bs_jdl, bfield_bs_cyl, etc.
!
!  Coils 1:12 --> IU30, IU90, IU150, IU210, IU270, IU330, IL30, IL90, IL150, IL210, IL270, IL330
!    Coil toroidal angles are in RH cylindrical coordinates, so negative
!    of machine geographical toroidal angle.  <-- Following M. Shaffer routines
!
! Input:
!  taper : real64(12) : The coil currents in Amps.  Number of turns (1) is applied in the routine
!  ntorpts : int32 : number of toroidal points used to generate the coils.  (6 was used in Schaffer's routine)
!
! Output:
!  coil : real64(ncoil_pts,3) : Array of coil filament beginning and end points,
!                             in Cartesian coordinates (m), second index is (1,2,3) = x,y,z
!  current : real64(ncoil_pts) : The current (in Amp-turns) corresponding to the current filaments. 
!                              The filament defined from coil(i,:) to coil(i+1,:) has a current
!                              set by current(i). The last value of array current is not used.  
!  ncoil_pts : int32 : number of coil points
!
! Calls:
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!   
! Author(s): J.D. Lore
!  -- Built on code by J.M. Canik and M. Schaffer
Use kind_mod, Only: real64, int32
Use phys_const, Only: pi  
Implicit None
Integer(int32), Parameter :: nc = 12
Real(real64), Intent(in) :: taper(nc)
Integer(int32), Intent(in), Optional :: ntorpts
Type(coil_type), Intent(Out) :: coil

Integer(int32) :: npts,nturn,i,j,myntorpts
Real(real64) :: phicens(nc),phiext,R(2),Z(2),taper2(nc)
Real(real64), Allocatable :: phi(:)
!- End of header -------------------------------------------------------------

If (Present(ntorpts)) Then
  myntorpts = ntorpts
Else
  myntorpts = 6
Endif

Call allocate_coiltype(nc,myntorpts,coil)

! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------
!             COIL DATA
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------
! Coils are approximately rectangular.
phicens = (/-32.7d0,-87.3d0,-152.7d0,-207.3d0,-272.7d0,-327.3d0, &
            -32.7d0,-87.3d0,-152.7d0,-207.3d0,-272.7d0,-327.3d0/) ! Center of each coil
phicens = phicens*pi/180.d0
phiext = 51.72d0*pi/180.d0               ! Toroidal extent of each coil
! 2003
!R = (/2.184d0,2.394d0/)       ! Major radius of coil  --- Upper coil, lower just invert and flip Z
!Z = (/1.012d0,0.504d0/)       ! Vertical extent of coil.
! 2006 'revised'
R = (/2.164d0,2.373d0/)        ! Major radius of coil  --- Upper coil, lower just invert and flip Z
Z = (/1.016d0,0.504d0/)        ! Vertical extent of coil.
nturn = 1
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------

taper2 = -1.d0*Real(nturn,real64)*taper   ! -1 accounts for DIII-D convention

Allocate(phi(myntorpts))
npts = coil%ncoil_pts/coil%num_coils

! UPPER COILS
Do i = 0,5
  Do j = 1,myntorpts
    phi(j) = phiext*Real(j-1,real64)/Real(myntorpts-1,real64) + phicens(i+1) - phiext/2.d0
  Enddo
  coil%coilxyz(i*npts+1:i*npts+myntorpts,1) = R(1)*cos(phi)
  coil%coilxyz(i*npts+1:i*npts+myntorpts,2) = R(1)*sin(phi)
  coil%coilxyz(i*npts+1:i*npts+myntorpts,3) = Z(1)
  coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,1) = R(2)*cos(phi(myntorpts:1:-1))
  coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,2) = R(2)*sin(phi(myntorpts:1:-1))
  coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,3) = Z(2)
  coil%coilxyz((i+1)*npts,:) = coil%coilxyz(i*npts+1,:)
  coil%current(i*npts+1:(i+1)*npts) = taper2(i+1)
  coil%current((i+1)*npts) = 0.0d0   ! Sticks connecting the coils have no current.
Enddo
! LOWER COILS
Do i = 6,11
  Do j = 1,myntorpts
    phi(j) = phiext*real(j-1,real64)/real(myntorpts-1,real64) + phicens(i+1) - phiext/2.d0
  Enddo
    coil%coilxyz(i*npts+1:i*npts+myntorpts,1) = R(2)*cos(phi)
    coil%coilxyz(i*npts+1:i*npts+myntorpts,2) = R(2)*sin(phi)
    coil%coilxyz(i*npts+1:i*npts+myntorpts,3) = -1.d0*Z(2)
    coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,1) = R(1)*cos(phi(myntorpts:1:-1))
    coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,2) = R(1)*sin(phi(myntorpts:1:-1))
    coil%coilxyz(i*npts+myntorpts+1:i*npts+2*myntorpts,3) = -1.d0*Z(1)
    coil%coilxyz((i+1)*npts,:) = coil%coilxyz(i*npts+1,:)
    coil%current(i*npts+1:(i+1)*npts) = taper2(i+1)
    coil%current((i+1)*npts) = 0.0d0   ! Sticks connecting the coils have no current.
Enddo
Deallocate(phi)

End Subroutine build_d3d_icoils_jl


End Module DIIID_routines_mod
