!-----------------------------------------------------------------------------
!
!   Routines/modules related to NSTX
!   
!-----------------------------------------------------------------------------
Module NSTX_routines_mod
Use kind_mod
Implicit None
Contains

Subroutine build_nstx_rwmcoils_jl(taper,ntorpts,coil,current,ncoil_pts)
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
Use kind_mod
Use phys_const, Only : pi
Implicit None

Real(rknd), Intent(in) :: taper(6)
Integer(iknd), Intent(in) :: ntorpts
Integer(iknd), Intent(out) :: ncoil_pts  
Real(rknd),Intent(out) :: coil(6*(2*ntorpts+1),3),current(6*(2*ntorpts+1))

Integer(iknd) :: npts,nturn,i,j
Real(rknd) :: phicens(6),phiext,R(2),Z(2),phi(ntorpts),taper2(6)


!if (size(taper) .eq. 1)
!    taper = taper*[1,-1,1,-1,1,-1];
!end
!if nargin < 2 
!    ntorpts = 5;
!end
npts = 2*ntorpts + 1
ncoil_pts = 6*npts

! Coils are approximately rectangular.
phicens = (/0.d0,300.d0,240.d0,180.d0,120.d0,60.d0/)  ! Center of each coil
phicens = phicens*pi/180.d0
phiext = 56.d0*pi/180.d0               ! Toroidal extent of each coil
R = (/1.76d0,1.76d0/)                          ! Major radius of coil
Z = (/0.4826d0,-0.4826d0/)                     ! Vertical extent of coil.
nturn = 2

coil = 0.d0
current = 0.d0

taper2 = -Real(nturn,rknd)*taper;  ! Minus sign accounts for the handedness of the coils as I have 
                       ! defined them relative to the 'real' orientation described above.
		       ! I.e., the returned coils are CW as viewed from outside NSTX.

Do i = 0,5
  Do j = 1,ntorpts
    phi(j) = phiext*real(j-1,rknd)/real(ntorpts-1,rknd) + phicens(i+1) - phiext/2.d0
  Enddo

    coil(i*npts+1:i*npts+ntorpts,1) = R(1)*cos(phi)
    coil(i*npts+1:i*npts+ntorpts,2) = R(1)*sin(phi)
    coil(i*npts+1:i*npts+ntorpts,3) = Z(1)
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,1) = R(2)*cos(phi(ntorpts:1:-1))
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,2) = R(2)*sin(phi(ntorpts:1:-1))
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,3) = Z(2)
    coil((i+1)*npts,:) = coil(i*npts+1,:)
    current(i*npts+1:(i+1)*npts) = taper2(i+1)
    current((i+1)*npts) = 0.0   ! Sticks connecting the coils have no current.
Enddo

End Subroutine build_nstx_rwmcoils_jl

End Module NSTX_routines_mod
