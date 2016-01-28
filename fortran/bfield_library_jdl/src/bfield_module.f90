!-----------------------------------------------------------------------------
!
!   Routines/modules related to bfield calculations
!
!   Contains:
!     Subroutine bfield_bs_cyl
!     Subroutine bfield_bs_jdl
!   
!-----------------------------------------------------------------------------
Module bfield_module
  Implicit None
  Private
  Public :: bfield_bs_cyl
  Public :: bfield_bs_jdl
Contains

!-----------------------------------------------------------------------------
!+ Performs Biot-Savart in cylindrical coordinates
!-----------------------------------------------------------------------------
Subroutine bfield_bs_cyl(P_r,P_phi,P_z,npts,coil,current,ncoil_pts,Br,Bphi,Bz)
!
! Description:
!  Performs Biot-Savart evaluation at a set of points given a series of 
!  current filaments. Points are in cylindrical coordinates (m,radians),
!  output in Tesla.
!
! Input:
!  P_r, P_phi, P_z : real64(1:npts) : Evaluation points (m, radians, m)
!  npts : int32 : number of evaluation points
!  coil : real64(ncoil_pts,3) : Array of coil filament beginning and end points,
!                             in Cartesian coordinates (m), second index is (1,2,3) = x,y,z
!  current : real64(ncoil_pts) : The current (in Amp-turns) corresponding to the current filaments. 
!                              The filament defined from coil(i,:) to coil(i+1,:) has a current
!                              set by current(i). The last value of array current is not used.  
!  ncoil_pts : int32 : number of coil points
!
! Output:
!  Br,Bphi,Bz : real64(npts) : Field components (Tesla)
!
! Calls:
!  bfield_bs_jdl
!
! Author(s): J.D. Lore
Use kind_mod, Only: real64, int32
Implicit None
!Input/output
Real(real64), Intent(in), Dimension(npts) :: P_r, P_phi, P_z
Real(real64), Intent(in), Dimension(ncoil_pts) :: current
Real(real64), Intent(in), Dimension(ncoil_pts,3) :: coil
Integer(int32), Intent(in) :: npts, ncoil_pts
Real(real64), Intent(out),Dimension(npts) :: Br, Bphi, Bz
! Local Variables
Real(real64), Dimension(npts) :: Bx, By, cp, sp 
Real(real64), Dimension(npts) :: P_x, P_y
!- End of header -------------------------------------------------------------

cp = cos(P_phi)
sp = sin(P_phi)

P_x = P_r*cp
P_y = P_r*sp

Call bfield_bs_jdl(P_x,P_y,P_z,npts,coil,current,ncoil_pts,Bx,By,Bz)

Br = Bx*cp + By*sp
Bphi = -Bx*sp + By*cp

End Subroutine bfield_bs_cyl

!-----------------------------------------------------------------------------
!+ Performs Biot-Savart in Cartesian coordinates
!-----------------------------------------------------------------------------
Subroutine bfield_bs_jdl(P_x,P_y,P_z,npts,coil,current,ncoil_pts,Bx,By,Bz)
!
! Description:
!  Performs Biot-Savart evaluation at a set of points given a series of 
!  current filaments. Points are in Cartesian coordinates (m,radians),
!  output in Tesla.
!
! Input:
!  P_r, P_y, P_z : real64(1:npts) : Evaluation points (m)
!  npts : int32 : number of evaluation points
!  coil : real64(ncoil_pts,3) : Array of coil filament beginning and end points,
!                             in Cartesian coordinates (m), second index is (1,2,3) = x,y,z
!  current : real64(ncoil_pts) : The current (in Amp-turns) corresponding to the current filaments. 
!                              The filament defined from coil(i,:) to coil(i+1,:) has a current
!                              set by current(i). The last value of array current is not used.  
!  ncoil_pts : int32 : number of coil points
!
! Output:
!  Bx,By,Bz : real64(npts) : Field components (Tesla)
!
! Calls:
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!
! Author(s): J.D. Lore 
Use kind_mod, Only: real64, int32
Implicit None
Real(real64), Intent(in), Dimension(npts) :: P_x, P_y, P_z
Real(real64), Intent(in), Dimension(ncoil_pts) :: current
Real(real64), Intent(in), Dimension(ncoil_pts,3) :: coil
Integer(int32), Intent(in) :: npts, ncoil_pts
Real(real64), Intent(out),Dimension(npts) :: Bx, By, Bz

Real(real64), Dimension(ncoil_pts-1) :: Ri_x, Ri_y, Ri_z, Ri_mag
Real(real64), Dimension(ncoil_pts-1) :: Rf_x, Rf_y, Rf_z, Rf_mag
Real(real64), Dimension(ncoil_pts-1) :: RicrRf_x, RicrRf_y, RicrRf_z
Real(real64), Dimension(ncoil_pts-1) :: RidotRf, back, Ival
Integer(int32) :: i, nfil
!- End of header -------------------------------------------------------------

nfil = ncoil_pts - 1
Ival=current(1:nfil)*1.e-7_real64 ! this is mu0*I/4pi
Do i = 1,npts

    Ri_x = P_x(i) - coil(1:nfil,1)     !R_i is from beg of filament to point x
    Ri_y = P_y(i) - coil(1:nfil,2)
    Ri_z = P_z(i) - coil(1:nfil,3)
    Ri_mag = sqrt(Ri_x*Ri_x + Ri_y*Ri_y + Ri_z*Ri_z)
    
    Rf_x(1:nfil-1) = Ri_x(2:nfil)
    Rf_x(nfil) = P_x(i) - coil(nfil+1,1)                    !R_f is from end of filament to x
    Rf_y(1:nfil-1) = Ri_y(2:nfil)
    Rf_y(nfil) = P_y(i) - coil(nfil+1,2)                    !it is thus the same as Ri except one point before
    Rf_z(1:nfil-1) = Ri_z(2:nfil)
    Rf_z(nfil) = P_z(i) - coil(nfil+1,3)
    Rf_mag(1:nfil-1) = Ri_mag(2:nfil)
    Rf_mag(nfil) = sqrt(Rf_x(nfil)*Rf_x(nfil) + Rf_y(nfil)*Rf_y(nfil) + Rf_z(nfil)*Rf_z(nfil))
    
    RicrRf_x=Ri_y*Rf_z - Ri_z*Rf_y      !cross product: RixRf
    RicrRf_y=Ri_z*Rf_x - Ri_x*Rf_z
    RicrRf_z=Ri_x*Rf_y - Ri_y*Rf_x
    
    RidotRf = Ri_x*Rf_x + Ri_y*Rf_y + Ri_z*Rf_z
    
    back = Ival(1:nfil)*(Ri_mag+Rf_mag)/(Ri_mag*Rf_mag*(Ri_mag*Rf_mag+RidotRf))
    
    Bx(i) = sum(RicrRf_x*back)
    By(i) = sum(RicrRf_y*back)
    Bz(i) = sum(RicrRf_z*back)
    
Enddo

End Subroutine Bfield_bs_jdl

End Module bfield_module
