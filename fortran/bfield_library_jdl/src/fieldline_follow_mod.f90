!-----------------------------------------------------------------------------
!
!   Routines/modules related to fieldline following
!
!   Contains:
!     Subroutine follow_fieldlines_rzphi_AS
!     Subroutine follow_fieldlines_rzphi
!     Subroutine fl_derivs_fun
!     Subroutine rk45_fixed_step_integrate
!     Subroutine rk4_core
!   
!   Set bfield_method to control the fieldline deriviative calls.
!   Appropriate loading must be done before calls to any fieldline following routine 
!   (e.g., gfile loading, rmp coil generation)
!
!    bfield_method == 
!                     0 -- gfile field only
!                     1 -- gfile + rmp coil field
!                     2 -- gfile + Pavel's screened fields with bspline interpolation (not fully implemented!)
!                     3 -- M3DC1 fields (not fully implemented!)
! 
!-----------------------------------------------------------------------------
Module fieldline_follow_mod
Implicit None
Integer :: bfield_method ! Used to select fl derivs
Contains

!-----------------------------------------------------------------------------
!+ Follows fieldlines in cylindrical coords --> Assuming an axisymmetric field!
!-----------------------------------------------------------------------------
Subroutine follow_fieldlines_rzphi_AS(rstart,zstart,phistart,Npts,dphi,nsteps,r,z,phi,ierr,i_last_good)
!
! Description: 
!  Follows fieldlines by integrating along toroidal angle in cylindrical coordinates. This version
!  assumes deriviatives are not a function of phi, and is therefore a little faster than subroutine
!  follow_fieldlines_rzphi. bfield_method is ignored for this routine (only gfile fields used)
!
!  Integration method is hard-coded fixed-step size RK4.
!
! Input:
!  rstart, zstart, phistart : rknd(Npts) : Launching points for fieldlines (m,m,radians)
!  Npts : iknd : Number of starting points
!  dphi : rknd : Integration step size (radians)
!  nsteps : iknd : Number of integration steps to take
! Output:
!  r,z,phi : rknd(Npts,nsteps+1) : Field line trajectories (m,m,radians)
!  ierr : iknd(Npts) : Error flag for each fl. 0 indicates no error, 1 indicates an error from bfield_geq_bicub
!  i_last_good : iknd(Npts) : If an error occured for the fl, this array gives the last 'good' index, i.e., before the 
!                             error occured
!
! Calls:
!  Subroutine bfield_geq_bicub
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Ported from Canik's g3d.  JDL
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod                ! Import rknd, iknd specifications
Use g3d_module, Only : bfield_geq_bicub
Implicit None

! Input/output                !See above for descriptions
Integer(iknd),Intent(in) :: Npts, nsteps
Real(rknd),Intent(in),Dimension(Npts) :: & 
  rstart(Npts),zstart(Npts),phistart(Npts)
Real(rknd),Intent(in) :: dphi
Real(rknd),Intent(out),Dimension(Npts,nsteps+1) :: &
  r,z,phi
Integer(iknd),Intent(out) :: ierr(Npts), i_last_good(Npts)

! Local variables
Real(rknd) :: br,bz,bphi
Real(Rknd) :: k1r,k2r,k3r,k4r,k1z,k2z,k3z,k4z
Integer(iknd) :: ii,ipt,ierr_b
Real(rknd) :: Bval(1,3),R1(1), Z1(1)
!- End of header -------------------------------------------------------------

r(:,1)   = rstart(:)
z(:,1)   = zstart(:)
phi(:,1) = phistart(:)

Do ipt = 1,Npts 
  ierr(ipt) = 0
  i_last_good(ipt) = nsteps + 1
  Do ii=2,nsteps+1
    phi(ipt,ii) = phi(ipt,ii-1) + dphi

    R1 = r(ipt,ii-1)
    Z1 = z(ipt,ii-1)
    Call bfield_geq_bicub(R1,Z1,1,Bval,ierr_b)     
    br = Bval(1,1)
    bz = Bval(1,2)
    bphi = Bval(1,3)
    If (ierr_b .ne. 0) Then
      ierr(ipt) = 1
      i_last_good(ipt) = ii - 1
      Exit
    Endif

    k1r = dphi*r(ipt,ii-1)*br/bphi
    k1z = dphi*r(ipt,ii-1)*bz/bphi
    
    R1 = r(ipt,ii-1) + 0.5_rknd*k1r
    Z1 = z(ipt,ii-1) + 0.5_rknd*k1z
    Call bfield_geq_bicub(R1,Z1,1,Bval,ierr_b)     
    br = Bval(1,1)
    bz = Bval(1,2)
    bphi = Bval(1,3)
    If (ierr_b .ne. 0) Then
      ierr(ipt) = 1
      i_last_good(ipt) = ii - 1
      Exit
    Endif

    k2r = dphi*(r(ipt,ii-1)+0.5_rknd*k1r)*br/bphi
    k2z = dphi*(r(ipt,ii-1)+0.5_rknd*k1r)*bz/bphi

    R1 = r(ipt,ii-1) + 0.5_rknd*k2r
    Z1 = z(ipt,ii-1) + 0.5_rknd*k2z
    Call bfield_geq_bicub(R1,Z1,1,Bval,ierr_b)     
    br = Bval(1,1)
    bz = Bval(1,2)
    bphi = Bval(1,3)
    If (ierr_b .ne. 0) Then
      ierr(ipt) = 1
      i_last_good(ipt) = ii - 1
      Exit
    Endif

    k3r = dphi*(r(ipt,ii-1)+0.5_rknd*k2r)*br/bphi
    k3z = dphi*(r(ipt,ii-1)+0.5_rknd*k2r)*bz/bphi


    R1 = r(ipt,ii-1) + k3r
    Z1 = z(ipt,ii-1) + k3z
    Call bfield_geq_bicub(R1,Z1,1,Bval,ierr_b)     
    br = Bval(1,1)
    bz = Bval(1,2)
    bphi = Bval(1,3)  
    If (ierr_b .ne. 0) Then
      ierr(ipt) = 1
      i_last_good(ipt) = ii - 1
      Exit
    Endif
    k4r = dphi*(r(ipt,ii-1)+k3r)*br/bphi
    k4z = dphi*(r(ipt,ii-1)+k3r)*bz/bphi

    r(ipt,ii) = r(ipt,ii-1) + k1r/6._rknd + k2r/3._rknd + k3r/3._rknd + k4r/6._rknd
    z(ipt,ii) = z(ipt,ii-1) + k1z/6._rknd + k2z/3._rknd + k3z/3._rknd + k4z/6._rknd

  Enddo
Enddo


EndSubroutine follow_fieldlines_rzphi_AS

!-----------------------------------------------------------------------------
!+ Follows fieldlines in cylindrical coords with diffusion
!-----------------------------------------------------------------------------
Subroutine follow_fieldlines_rzphi_diffuse(rstart,zstart,phistart,Npts,dphi,nsteps,r,z,phi,ierr,i_last_good,dmag)
!
! Description: 
!  Follows fieldlines by integrating along toroidal angle in cylindrical coordinates. 
!
!  Integration method is RK45 with fixed step size
!
! Input:
!  rstart, zstart, phistart : rknd(Npts) : Launching points for fieldlines (m,m,radians)
!  Npts : iknd : Number of starting points
!  dphi : rknd : Integration step size (radians)
!  nsteps : iknd : Number of integration steps to take
!  dmag : magnetic diffusivity (m^2/m)
! Output:
!  r,z,phi : rknd(Npts,nsteps+1) : Field line trajectories (m,m,radians)
!  ierr : iknd(Npts) : Error flag for each fl. 0 indicates no error, 1 indicates an error from bfield_geq_bicub
!  i_last_good : iknd(Npts) : If an error occured for the fl, this array gives the last 'good' index, i.e., before the 
!                             error occured
!
! Calls:
!  Subroutine rk45_fixed_step_integrate
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod                ! Import rknd, iknd specifications
Implicit None
! Input/output                      !See above for descriptions
Integer(iknd),Intent(in) :: Npts, nsteps
Real(rknd),Intent(in),Dimension(Npts) :: rstart(Npts),zstart(Npts),phistart(Npts)
Real(rknd),Intent(in) :: dphi, dmag

Real(rknd),Intent(out),Dimension(Npts,nsteps+1) :: r,z,phi
Integer(iknd),Intent(out) :: ierr(Npts), i_last_good(Npts)
! Local variables
Integer(iknd), Parameter :: n = 2
Real(rknd) :: y(n),x,dx,xout(nsteps+1),yout(n,nsteps+1)
Integer(iknd) :: ierr_rk45, i_last_good_rk45, ipt
!- End of header -------------------------------------------------------------

r = 0._rknd
z = 0._rknd
phi = 0._rknd
dx = dphi
Do ipt = 1,Npts 
  y(1) = rstart(ipt)
  y(2) = zstart(ipt)
  x = phistart(ipt)
  Call rk45_fixed_step_integrate_diffuse(y,n,x,dx,nsteps,fl_derivs_fun,yout,xout,ierr_rk45,i_last_good_rk45,dmag)
  r(ipt,1:nsteps+1)   = yout(1,1:nsteps+1)
  z(ipt,1:nsteps+1)   = yout(2,1:nsteps+1)
  phi(ipt,1:nsteps+1) = xout
  ierr(ipt) = ierr_rk45
  i_last_good(ipt) = i_last_good_rk45
Enddo

EndSubroutine follow_fieldlines_rzphi_diffuse


!-----------------------------------------------------------------------------
!+ Follows fieldlines in cylindrical coords
!-----------------------------------------------------------------------------
Subroutine follow_fieldlines_rzphi(rstart,zstart,phistart,Npts,dphi,nsteps,r,z,phi,ierr,i_last_good)
!
! Description: 
!  Follows fieldlines by integrating along toroidal angle in cylindrical coordinates. 
!
!  Integration method is RK45 with fixed step size
!
! Input:
!  rstart, zstart, phistart : rknd(Npts) : Launching points for fieldlines (m,m,radians)
!  Npts : iknd : Number of starting points
!  dphi : rknd : Integration step size (radians)
!  nsteps : iknd : Number of integration steps to take
! Output:
!  r,z,phi : rknd(Npts,nsteps+1) : Field line trajectories (m,m,radians)
!  ierr : iknd(Npts) : Error flag for each fl. 0 indicates no error, 1 indicates an error from bfield_geq_bicub
!  i_last_good : iknd(Npts) : If an error occured for the fl, this array gives the last 'good' index, i.e., before the 
!                             error occured
!
! Calls:
!  Subroutine rk45_fixed_step_integrate
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod                ! Import rknd, iknd specifications
Implicit None
! Input/output                      !See above for descriptions
Integer(iknd),Intent(in) :: Npts, nsteps
Real(rknd),Intent(in),Dimension(Npts) :: rstart(Npts),zstart(Npts),phistart(Npts)
Real(rknd),Intent(in) :: dphi

Real(rknd),Intent(out),Dimension(Npts,nsteps+1) :: r,z,phi
Integer(iknd),Intent(out) :: ierr(Npts), i_last_good(Npts)
! Local variables
Integer(iknd), Parameter :: n = 2
Real(rknd) :: y(n),x,dx,xout(nsteps+1),yout(n,nsteps+1)
Integer(iknd) :: ierr_rk45, i_last_good_rk45, ipt
!- End of header -------------------------------------------------------------

r = 0._rknd
z = 0._rknd
phi = 0._rknd
dx = dphi
Do ipt = 1,Npts 
  y(1) = rstart(ipt)
  y(2) = zstart(ipt)
  x = phistart(ipt)
  Call rk45_fixed_step_integrate(y,n,x,dx,nsteps,fl_derivs_fun,yout,xout,ierr_rk45,i_last_good_rk45)
  r(ipt,1:nsteps+1)   = yout(1,1:nsteps+1)
  z(ipt,1:nsteps+1)   = yout(2,1:nsteps+1)
  phi(ipt,1:nsteps+1) = xout
  ierr(ipt) = ierr_rk45
  i_last_good(ipt) = i_last_good_rk45
Enddo

EndSubroutine follow_fieldlines_rzphi

!-----------------------------------------------------------------------------
!+ Routines field line equation derivatives
!-----------------------------------------------------------------------------
Subroutine fl_derivs_fun(n,phi,RZ,df,ierr)
!
! Description: 
!  Evaluates field line deriviatives (based on bfield_method). Should be easy to generalize to different
!  integration methods, but right now it is hard-coded to assume two simultaneous equations evaluated at 1 pt.
!
! Input:
!  n : iknd : Number of equations
!  phi : rknd : Integration variable (radians)
!  RZ : rknd(n) : Evaluation points (solution vector)  RZ(1) = R, RZ(2) = Z in meters
! 
! Output:
!  df : rknd(n) : derivative evaluation
!  ierr : iknd : Error flag (0 = no error)
!
! Calls:
!  Subroutine bfield_geq_bicub
!  Subroutine bfield_bs_cyl
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod
Use g3d_module, Only : bfield_geq_bicub
Use rmpcoil_module, Only : rmp_coil, rmp_coil_current, rmp_ncoil_pts
Use screening_module, Only : bfield_bspline
Use bfield_module, Only : bfield_bs_cyl
Use m3dc1_routines_mod, Only : bfield_m3dc1
Implicit None
Real(rknd), Intent(In) :: phi
Integer(iknd), Intent(In) :: n
Integer(iknd), Intent(Out) :: ierr
Real(rknd), Intent(In), Dimension(n) :: RZ
Real(rknd), Intent(Out), Dimension(n) :: df

Integer(iknd),Parameter :: Npts = 1
Real(rknd) :: bval(Npts,3), phi_tmp(Npts), bval_screened(Npts,3), bval_tmp(Npts,3)
Integer(iknd) :: ierr_b, ierr_rmp
Real(rknd) :: Bz, Br, Bphi, Br_rmp(Npts), Bphi_rmp(Npts), Bz_rmp(Npts)
!- End of header -------------------------------------------------------------

bval = 0._rknd
ierr_b = 0
ierr_rmp = 0
If (bfield_method == 0) Then     ! gfile field only
  Call bfield_geq_bicub(RZ(1),RZ(2),Npts,bval,ierr_b)     
  Br   = bval(1,1)
  Bz   = bval(1,2)
  Bphi = bval(1,3)
Elseif (bfield_method == 1) Then ! g + rmp coils
  Call bfield_geq_bicub(RZ(1),RZ(2),Npts,bval,ierr_b)     
  phi_tmp(1) = phi
!  If (Allocated(rmp_coil)) Then
    Call bfield_bs_cyl(RZ(1),phi_tmp,RZ(2),Npts,rmp_coil,rmp_coil_current,rmp_ncoil_pts,Br_rmp,Bphi_rmp,Bz_rmp)
!  Else
!    Write(*,*) 'RMP VARIABLES NOT ALLOCATED, EXITING FROM fl_derivs_fun!'
!    Stop
!  Endif
  ierr_rmp = 0
  Br   = bval(1,1) + Br_rmp(1)
  Bz   = bval(1,2) + Bz_rmp(1)
  Bphi = bval(1,3) + Bphi_rmp(1)
Elseif (bfield_method == 2) Then ! g + screening B-spline
  Call bfield_geq_bicub(RZ(1),RZ(2),Npts,bval,ierr_b)     
  phi_tmp(1) = phi
  Call bfield_bspline(RZ(1),phi_tmp,RZ(2),Npts,bval_screened,ierr_rmp)
  Br   = bval(1,1) + bval_screened(1,1)
  Bz   = bval(1,2) + bval_screened(1,2)
  Bphi = bval(1,3) + bval_screened(1,3)
Elseif (bfield_method == 3) Then  ! m3dc1 
  phi_tmp(1) = phi
  bval_tmp =0.d0
  Call bfield_geq_bicub(RZ(1),RZ(2),Npts,bval,ierr_b)     
  Call bfield_m3dc1(RZ(1),phi_tmp(1),RZ(2),Npts,bval_tmp,ierr_rmp)
  bval = bval + bval_tmp
  Br   = bval(1,1)
  Bz   = bval(1,2)
  Bphi = bval(1,3)
  ierr_rmp = ierr_b + ierr_rmp
Else
  Write(*,*) 'Unknown bfield_method in fl_derivs_fun'
  stop
Endif

If ((ierr_b .ne. 0) .or. (ierr_rmp .ne. 0)) Then
  ierr = 1
  df = 0._rknd
  Return
Else
  ierr = 0
Endif
df(1) = RZ(1)*Br/Bphi
df(2) = RZ(1)*Bz/Bphi

End Subroutine fl_derivs_fun


!-----------------------------------------------------------------------------
!+ Main routine for RK45 fixed step integration
!-----------------------------------------------------------------------------
Subroutine rk45_fixed_step_integrate(y0,n,x0,dx,nsteps,odefun,yout,xout,ierr,i_last_good)
!
! Description: 
!  Should be a general implementation of RK45 fixed step integration
!
! Input:
!  y0     : rknd(n) : initial values
!  n      : iknd    : number of initial values
!  x0     : rknd    : Location of initial values
!  dx     : rknd    : Step size
!  nsteps : iknd    : Number of integration steps
!  odefun : External function that evaluates derivatives
! 
! Output:
!  yout : rknd(n,nsteps+1) : Solution
!  xout : rknd(nsteps+1)   : Solution evaluation locations
!  ierr : iknd             : Error flag (0 = no error)
!  i_last_good : iknd      : Index of last good evaluation before error
!
! Calls:
!  Subroutine odefun
!  Subroutine rk4_core
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod
Implicit None
Real(rknd), Intent(In), Dimension(n) :: y0
Integer(iknd), Intent(In) :: n, nsteps
Real(rknd), Intent(In) :: x0, dx
Real(rknd), Intent(Out), Dimension(n,nsteps+1) :: yout
Real(rknd), Intent(Out), Dimension(nsteps+1) :: xout
Integer(iknd), Intent(Out) :: ierr, i_last_good

Integer(iknd) :: i, ierr_odefun, ierr_rk4core
Real(rknd), Dimension(n) :: y, dydx, ytmp
Real(rknd) :: x
External odefun
!- End of header -------------------------------------------------------------

! Store initial point
yout(:,:) = 0._rknd
xout(:) = 0._rknd
yout(:,1) = y0
xout(1) = x0

y = y0
x = x0
ierr = 0
i_last_good = nsteps+1
Do i=1,nsteps
  Call odefun(n,x,y,dydx,ierr_odefun)
  If (ierr_odefun == 1) Then
    ierr = 1
    i_last_good = i
    Return
  Endif
  Call rk4_core(y,dydx,n,x,dx,odefun,ytmp,ierr_rk4core)
  If (ierr_rk4core == 1) Then
    ierr = 1
    i_last_good = i
    Return
  Endif

  x = x + dx

  yout(:,i+1) = ytmp
  xout(i+1) = x
  y = ytmp
Enddo

End Subroutine rk45_fixed_step_integrate

!-----------------------------------------------------------------------------
!+ Main routine for RK45 fixed step integration with diffusion
!-----------------------------------------------------------------------------
Subroutine rk45_fixed_step_integrate_diffuse(y0,n,x0,dx,nsteps,odefun,yout,xout,ierr,i_last_good,dmag)
!
! Description: 
!  Should be a general implementation of RK45 fixed step integration
!
! Input:
!  y0     : rknd(n) : initial values
!  n      : iknd    : number of initial values
!  x0     : rknd    : Location of initial value
!  dx     : rknd    : Step size
!  dmag   : rknd    : Magnetic diffusivity (m^2/m)
!  nsteps : iknd    : Number of integration steps
!  odefun : External function that evaluates derivatives
! 
! Output:
!  yout : rknd(n,nsteps+1) : Solution
!  xout : rknd(nsteps+1)   : Solution evaluation locations
!  ierr : iknd             : Error flag (0 = no error)
!  i_last_good : iknd      : Index of last good evaluation before error
!
! Calls:
!  Subroutine odefun
!  Subroutine rk4_core
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod
Use g3d_module, Only : bfield_geq_bicub
Use gfile_var_pass, Only : g_rmaxis, g_zmaxis
Use rmpcoil_module, Only : rmp_coil, rmp_coil_current, rmp_ncoil_pts
Use screening_module, Only : bfield_bspline
Use bfield_module, Only : bfield_bs_cyl
Use m3dc1_routines_mod, Only : bfield_m3dc1
Implicit None


Real(rknd), Parameter :: diff_mag = 0.4d0
Real(rknd), Parameter :: nfac_diff = 3.d0
!Real(rknd), Parameter :: sigtheta = 20.d0
Real(rknd), Parameter :: sigtheta = 20.d0*3.1415d0/180.d0
Real(rknd)  :: theta, phi_factor, pol_factor

Real(rknd), Parameter :: pi = 3.141592653589793238462643383279502_rknd
Real(rknd), Intent(In), Dimension(n) :: y0
Integer(iknd), Intent(In) :: n, nsteps
Real(rknd), Intent(In) :: x0, dx, dmag
Real(rknd), Intent(Out), Dimension(n,nsteps+1) :: yout
Real(rknd), Intent(Out), Dimension(nsteps+1) :: xout
Integer(iknd), Intent(Out) :: ierr, i_last_good
Logical, parameter :: verbose = .false.
Integer(iknd) :: i, ierr_odefun, ierr_rk4core
Real(rknd), Dimension(n) :: y, dydx, ytmp
Real(rknd) :: x, RZ(2), perpdir1(3), perpdir2(3), alpha, dca, dsa, delta_x, dL
Real(rknd) :: bval(1,3), phi_tmp(1), bval_screened(1,3), bval_tmp(1,3), phi, rnum
Integer(iknd) :: ierr_b, ierr_rmp
Real(rknd) :: Bz, Br, Bphi, Br_rmp(1), Bphi_rmp(1), Bz_rmp(1)
External odefun
!- End of header -------------------------------------------------------------

! Store initial point
yout(:,:) = 0._rknd
xout(:) = 0._rknd
yout(:,1) = y0
xout(1) = x0

y = y0
x = x0
ierr = 0
i_last_good = nsteps+1
Do i=1,nsteps
  Call odefun(n,x,y,dydx,ierr_odefun)
  If (ierr_odefun == 1) Then
    ierr = 1
    i_last_good = i
    Return
  Endif
  Call rk4_core(y,dydx,n,x,dx,odefun,ytmp,ierr_rk4core)
  If (ierr_rk4core == 1) Then
    ierr = 1
    i_last_good = i
    Return
  Endif


  x = x + dx

  !
  ! Diffuse
  ! 

  RZ(1) = ytmp(1)
  RZ(2) = ytmp(2)
  phi = x
  bval = 0._rknd
  ierr_b = 0
  ierr_rmp = 0
  If (bfield_method == 0) Then     ! gfile field only
    Call bfield_geq_bicub(RZ(1),RZ(2),1,bval,ierr_b,verbose)     !-- turn off verbose
    Br   = bval(1,1)
    Bz   = bval(1,2)
    Bphi = bval(1,3)
  Elseif (bfield_method == 1) Then ! g + rmp coils
    Call bfield_geq_bicub(RZ(1),RZ(2),1,bval,ierr_b,verbose)     
    phi_tmp(1) = phi
    Call bfield_bs_cyl(RZ(1),phi_tmp,RZ(2),1,rmp_coil,rmp_coil_current,rmp_ncoil_pts,Br_rmp,Bphi_rmp,Bz_rmp)
    ierr_rmp = 0
    Br   = bval(1,1) + Br_rmp(1)
    Bz   = bval(1,2) + Bz_rmp(1)
    Bphi = bval(1,3) + Bphi_rmp(1)
  Elseif (bfield_method == 2) Then ! g + screening B-spline
    Call bfield_geq_bicub(RZ(1),RZ(2),1,bval,ierr_b,verbose)     
    phi_tmp(1) = phi
    Call bfield_bspline(RZ(1),phi_tmp,RZ(2),1,bval_screened,ierr_rmp)
    Br   = bval(1,1) + bval_screened(1,1)
    Bz   = bval(1,2) + bval_screened(1,2)
    Bphi = bval(1,3) + bval_screened(1,3)
  Elseif (bfield_method == 3) Then  ! m3dc1 
    phi_tmp(1) = phi
    bval_tmp =0.d0
    Call bfield_geq_bicub(RZ(1),RZ(2),1,bval,ierr_b,verbose)     
    Call bfield_m3dc1(RZ(1),phi_tmp(1),RZ(2),1,bval_tmp,ierr_rmp)
    bval = bval + bval_tmp
    Br   = bval(1,1)
    Bz   = bval(1,2)
    Bphi = bval(1,3)
    ierr_rmp = ierr_b + ierr_rmp
  Else
    Write(*,*) 'Unknown bfield_method in rk45_fixed_step_integrate_diffuse'
    stop
  Endif
  If ((ierr_b .ne. 0) .or. (ierr_rmp .ne. 0)) Then
    ierr = 1
    i_last_good = i
    Return
  Endif
  
  dL = sqrt(ytmp(1)*ytmp(1) + y(1)*y(1) - 2._rknd*ytmp(1)*y(1)*cos(dx) & 
       + ytmp(2)*ytmp(2) + y(2)*y(2) - 2._rknd*ytmp(2)*y(2))
  
  ! B cross z^hat
  perpdir1(1) = Bphi   !r 
  perpdir1(2) = -Br    !phi
  perpdir1(3) = 0.d0   !z
  perpdir1 = perpdir1/Sqrt(bphi*bphi + br*br)
  
  ! B cross r^hat
  perpdir2(1) = 0.d0
  perpdir2(2) = bz
  perpdir2(3) = -bphi
  perpdir2 = perpdir2/Sqrt(bz*bz + bphi*bphi)

  Call Random_number(rnum)
  
  alpha = 2.d0*pi*(-1.d0 + 2.d0*rnum) ! -2pi to 2pi kick    
  dca = Cos(alpha)
  dsa = Sin(alpha)
  
  phi_factor = 1 + diff_mag*cos(nfac_diff*x)

!  theta = atan2(ytmp(1)-g_rmaxis,ytmp(2)-g_zmaxis)*180.d0/3.14159d0
!  pol_factor = (1/(2.506628274631*sigtheta))*exp(-theta**2/(2*sigtheta**2))*250.d0

  theta = atan2(ytmp(2)-g_zmaxis,ytmp(1)-g_rmaxis)
  pol_factor = (1/(2.506628274631*sigtheta))*exp(-theta**2/(2*sigtheta**2))

  
  delta_x = Sqrt(dmag*dL*phi_factor*pol_factor)
  
  
  ytmp(1)   = ytmp(1)   + delta_x*(dca*perpdir1(1) + dsa*perpdir2(1))
  x         = x         + delta_x*(dca*perpdir1(2) + dsa*perpdir2(2))
  ytmp(2)   = ytmp(2)   + delta_x*(dca*perpdir1(3) + dsa*perpdir2(3))

!  write(*,*) delta_x


  yout(:,i+1) = ytmp
  xout(i+1) = x
  y = ytmp


Enddo





End Subroutine rk45_fixed_step_integrate_diffuse


!-----------------------------------------------------------------------------
!+ RK4 stepper
!-----------------------------------------------------------------------------
Subroutine rk4_core(y,dydx,n,x,dx,odefun,yout,ierr)
! Advance y(x) to yout=y(x+dx) given dydx(x) using 
! the RK4 method. y, dydx are vectors of length n, 
! function odefun(n,x,y,dxdy) returns dydx. 
! JDL 5/2012
Use kind_mod
Implicit None
Real(rknd), Intent(In), Dimension(n) :: y, dydx
Real(rknd), Intent(Out), Dimension(n) :: yout
Real(rknd), Intent(In) :: x, dx
Integer(iknd), Intent(In) :: n
Integer(iknd), Intent(Out) :: ierr

Real(rknd), Parameter :: TWO = 2._rknd
Real(rknd), Parameter :: SIX = 6._rknd
Integer(iknd) :: ierr_odefun
Real(rknd), Dimension(n) :: d1,d2,d3,d4

External odefun
!- End of header -------------------------------------------------------------

! First step (uses supplied derivatives)
d1 = dx*dydx
ierr = 0
! Second step
Call odefun(n,x+dx/TWO,y+d1/TWO,dydx,ierr_odefun)
If (ierr_odefun == 1) Then
  ierr = 1
  yout = 0._rknd
  Return
Endif
d2 = dx*dydx

! Third step
Call odefun(n,x+dx/TWO,y+d2/TWO,dydx,ierr_odefun)
If (ierr_odefun == 1) Then
  ierr = 1
  yout = 0._rknd
  Return
Endif
d3 = dx*dydx

! Fourth step
Call odefun(n,x+dx,y+d3,dydx,ierr_odefun)
If (ierr_odefun == 1) Then
  ierr = 1
  yout = 0._rknd
  Return
Endif
d4 = dx*dydx

yout = y + (d1 + TWO*d2 + TWO*d3 + d4)/SIX

Return
End Subroutine rk4_core


End Module fieldline_follow_mod
