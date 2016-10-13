!-----------------------------------------------------------------------------
!
!   Routines/modules related to NSTX
!   
!-----------------------------------------------------------------------------
Module NSTX_routines_mod
  Use coil_typedef, Only : coil_type, allocate_coiltype
  Implicit None
  Private
  Public :: build_nstx_rwmcoils_jl
Contains

Subroutine build_nstx_rwmcoils_jl(coil,taper,ntorpts)
! Builds the coil geometry and sets the current (in Amp-turns) for 
! the 6 NSTX RWM coils. A single loop is returned for each entire 
! coil set, with one turn per RWM coil.
!
! Inputs:
!	taper(1:6) -- Contains the coil current in Amps for RWM#1-6. This
!                 value is multiplied by 2 turns/coil. If length(taper) == 1
!                 then this is the current for RWM1 and an n=3 phase is used.
!                 This should have the sign of the irwm# MDSplus signals.
!	ntorpts -- Number of toroidal cuts used to define each coil. Default 5.
!
! Output:
!	rwm.coil(npts,1:3) -- Coil set coordinates in meters. The second index corresponds
!	                      to the x,y, or z coordinate. 
!   rwm.current(npts)  -- Current in each filament. The filaments connecting the RWM 
!					      coils have 0 current, and the last value is irrelevant 
!					     (number of filaments = number of points - 1).
!
!
! Notes:
!	1) Should return geometry aligned to 'NSTX physics' toroidal angle.
!	   0 degrees is at the center of RWM coil 1, and phi increases CCW when
!      viewing NSTX from above, co-Ip, and opposite the direction of increasing
!      RWM coil number.
!	2) Viewed from outside the machine, each 'real' RWM coil current flows CCW for 
!	   a positive value of MDSplus signal "irwm#", thus Br is positive for 
!      positive "irwm#" signal. My coils have the opposite sense, so the taper values
!      are multiplied by -1 below to give the proper field orientation.
! 
! JDL
Use kind_mod, Only: real64, int32
Use phys_const, Only : pi
Implicit None

Integer, Parameter :: nc = 6
Real(real64), Intent(in) :: taper(nc)
Integer(int32), Intent(in), Optional :: ntorpts
Type(coil_type), Intent(out) :: coil

Integer(int32) :: npts,nturn,i,j,myntorpts
Real(real64) :: phicens(nc),phiext,R(2),Z(2),taper2(nc)
Real(real64), Allocatable :: phi(:)
!- End of header -------------------------------------------------------------

If (Present(ntorpts)) Then
  myntorpts = ntorpts
Else
  myntorpts = 5
Endif

Call allocate_coiltype(nc,myntorpts,coil)

! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------
!             COIL DATA
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------
! Coils are approximately rectangular.
phicens = (/0.d0,300.d0,240.d0,180.d0,120.d0,60.d0/)  ! Center of each coil
phicens = phicens*pi/180.d0
phiext = 56.d0*pi/180.d0               ! Toroidal extent of each coil
R = (/1.76d0,1.76d0/)                  ! Major radius of coil
Z = (/0.4826d0,-0.4826d0/)             ! Vertical extent of coil.
nturn = 2
! -------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------

taper2 = -Real(nturn,real64)*taper;  ! Minus sign accounts for the handedness of the coils as I have 
                       ! defined them relative to the 'real' orientation described above.
		       ! I.e., the returned coils are CW as viewed from outside NSTX.

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
  coil%current((i+1)*npts) = 0.0   ! Sticks connecting the coils have no current.
Enddo
Deallocate(phi)
End Subroutine build_nstx_rwmcoils_jl

End Module NSTX_routines_mod
