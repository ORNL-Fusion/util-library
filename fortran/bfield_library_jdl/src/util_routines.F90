!-----------------------------------------------------------------------------
!
!   Routines/modules with utility functions that are (ideally) not specific
!    to one field method
!
!   Contains:
!     Module util_routines
!       Subroutine calc_sep
!       Subroutine find_xpt_jdl
!-----------------------------------------------------------------------------
Module util_routines
Implicit None
Contains

!-----------------------------------------------------------------------------
!+ Calculates approximate curve(s) for the sepatrix
!-----------------------------------------------------------------------------
Subroutine calc_sep(second_sep,fname_out)
!
! Description:
!  Calculates the separatrix curve(s) 
!
!  bfield_method must be set first, and a gfile must also have been read
!
! Input:
! Output:
!
! Calls:
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!            6/4/15    Port from matlab
! Author(s): J.D. Lore 
Use kind_mod
Use fieldline_follow_mod, Only: follow_fieldlines_rzphi
Use gfile_var_pass, Only: g_ip_sign,g_lim,g_limitr,g_r,g_z,g_mw,g_mh
Use math_geo_module, Only : inside_poly
Use phys_const, Only: pi
Implicit None
!Input/output
Logical, Intent(In) :: second_sep
Character(Len=120),Intent(In) :: fname_out
! Local Variables
Real(rknd) :: rx,zx,rx2,zx2
Integer(iknd) :: idir,maxtries,ii,i,nlim,nbox,in1,in2,nsep1,nsep0,nseps,isep,icount
Real(rknd), Allocatable, Dimension(:) :: rsep0,zsep0,rsep0_tmp,zsep0_tmp,lim_r,lim_z
Real(rknd), Allocatable, Dimension(:) :: rsep1,zsep1
Real(rknd), Allocatable, Dimension(:) ::  box_r,box_z
! Fl folllowing vars
Real(rknd), Dimension(1) :: rstart,zstart,phistart
Integer(iknd), Dimension(1) :: fl_ierr,ilg
Real(rknd), Allocatable, Dimension(:,:) :: fl_r,fl_z,fl_p
Real(rknd) :: dphi_fl
Integer(iknd) :: nsteps_fl
! Local Parameters
Real(rknd), Parameter :: dx = 1.d-3
Integer(rknd), Parameter :: nsep_div = 1000  ! total number of sep points from fl follow divided by this to output
!- End of header -------------------------------------------------------------

Open(99,file=fname_out)
If (second_sep) Then
  nseps = 3
Else
  nseps = 1
Endif
Write(99,*) nseps

! Set up boundaries to search
! 1) gfile lim
Allocate(lim_r(g_limitr))
Allocate(lim_z(g_limitr))
! Throw away zero points 
nlim = 0
Do i = 1,g_limitr
  If (g_lim(1,i) > 1.d-4) Then
    nlim = nlim + 1
    lim_r(nlim) = g_lim(1,i)
    lim_z(nlim) = g_lim(2,i)
  Endif
Enddo
If (nlim == 0) Then
  Write(*,*) 'Error: could not create limiter bdry in calc_sep'
  Stop
Endif
! 2) box based on gfile domain size
nbox = 4
Allocate(box_r(nbox),box_z(nbox))
box_r = (/g_r(1),g_r(g_mw),g_r(g_mw),g_r(1)/)
box_z = (/g_z(1),g_z(1),g_z(g_mh),g_z(g_mh)/)


Call find_xpt_jdl(.true.,.true.,1.d-8,.false.,rx,zx,rx2,zx2)

dphi_fl = g_ip_sign*0.1d0*pi/180.d0
nsteps_fl = Floor(0.3d0*pi/Abs(dphi_fl))  ! Factor of 0.x to stop from leaving grid
write(*,*) dphi_fl,nsteps_fl
maxtries = 1000
Allocate(fl_r(1,nsteps_fl+1),fl_z(1,nsteps_fl+1),fl_p(1,nsteps_fl+1))
fl_r = 0.d0; fl_z = 0.d0; fl_p = 0.d0
phistart = 0.d0

Do isep=1,nseps
  Write(*,*) 'Working on sep ',isep,'of ',nseps
  Do idir=1,-1,-2
    Write(*,*) 'Working on direction ',idir
    If (isep == 1) Then
      rstart = rx-dx
      zstart = zx
    Elseif (isep == 2) Then
      rstart = rx2
      zstart = zx2-dx
    Else
      rstart = rx2
      zstart = zx2+dx      
    Endif
    
    Do ii = 1,maxtries
      Call follow_fieldlines_rzphi(rstart,zstart,phistart,1,dphi_fl*idir,nsteps_fl,fl_r,fl_z,fl_p,fl_ierr,ilg)
      Allocate(rsep0(nsteps_fl*ii),zsep0(nsteps_fl*ii))
      If ( ii .gt. 1) Then
        rsep0 = [rsep0_tmp,fl_r(1,2:nsteps_fl+1)]
        zsep0 = [zsep0_tmp,fl_z(1,2:nsteps_fl+1)]
        Deallocate(rsep0_tmp,zsep0_tmp)
      Else
        rsep0 = fl_r(1,2:nsteps_fl+1)
        zsep0 = fl_z(1,2:nsteps_fl+1)
      Endif
      rstart = fl_r(1,nsteps_fl+1)
      zstart = fl_z(1,nsteps_fl+1)
      
      in1 = inside_poly(rstart(1),zstart(1),lim_r,lim_z,nlim)
      in2 = inside_poly(rstart(1),zstart(1),box_r,box_z,nbox)    

      If ( in1 + in2 .NE. 2 ) Exit
      If ( ii == maxtries ) Then
        Write(*,*) 'Increase maxtries in calc_sep'
        Stop
      Endif
    
      Allocate(rsep0_tmp(nsteps_fl*ii),zsep0_tmp(nsteps_fl*ii))
      rsep0_tmp = rsep0
      zsep0_tmp = zsep0
      Deallocate(rsep0,zsep0)
    Enddo  ! ii to maxtries
    
    If (idir == 1) Then
      nsep1 = nsteps_fl*ii
      Allocate(rsep1(nsep1),zsep1(nsep1))
      rsep1 = rsep0
      zsep1 = zsep0
      Deallocate(rsep0,zsep0)
    Endif
  Enddo
  ! reverse second direction
  nsep0 = nsteps_fl*ii
  rsep0 = rsep0(nsep0:1:-1)
  zsep0 = zsep0(nsep0:1:-1)

  ! drop last point
  Write(99,*) (nsep0-1)/((nsep0+nsep1)/nsep_div)+nsep1/((nsep0+nsep1)/nsep_div)+2
  Do i = 1,nsep0-1,(nsep0+nsep1)/nsep_div
    Write(99,*) rsep0(i),zsep0(i)
  Enddo
  Do i = 1,nsep1,(nsep0+nsep1)/nsep_div
    Write(99,*) rsep1(i),zsep1(i)    
  Enddo
  Deallocate(rsep1,zsep1)
Enddo

Close(99)
Deallocate(fl_r,fl_z,fl_p)
Deallocate(lim_r,lim_z,box_r,box_z)


End Subroutine calc_sep







!-----------------------------------------------------------------------------
!+ Find x-point(s) from gfile or m3dc1 field
!-----------------------------------------------------------------------------
Subroutine find_xpt_jdl(second,refine,tol,quiet,rx,zx,rx2,zx2,phi_eval_deg,dx)
  ! tol is magnitude of bp at xpt
  ! phi_eval_deg only used for m3dc1 fields
  ! Both cases require a gfile field for initial guess!!!
Use kind_mod
Use gfile_var_pass, Only: g_bdry, g_r, g_z, g_mh, g_mw
Use math_geo_module, Only: rlinspace
Use fieldline_follow_mod, Only: bfield_method
#ifdef HAVE_M3DC1
Use m3dc1_routines_mod, Only: bfield_m3dc1
#endif
Use g3d_module, Only: bfield_geq_bicub
Use phys_const, Only: pi
Implicit None
Logical, Intent(in) :: second, refine, quiet
Real(rknd), Intent(in) :: tol
Real(rknd), Intent(in), Optional :: phi_eval_deg, dx
Real(rknd), Intent(out) :: rx, zx, rx2, zx2

Integer(iknd), Parameter :: niter_max = 15
Integer(iknd), Parameter :: n1 = 100  ! grid dimension (square)

Real(rknd), Allocatable :: rtmp(:), ztmp(:),Bout(:,:), Bout_tmp(:,:), phi_tmp(:)
Real(rknd) :: bp(n1,n1), rg(n1,n1), zg(n1,n1), rt(n1), zt(n1)
Real(rknd) :: bpx, err, de, bpx2, dx1_grid, dx2_grid, my_phi_eval
Integer(iknd) :: icount, i, npts_bdry, ierr, ix,ixjx(2), niter
! Local parameters               
!- End of header -------------------------------------------------------------
my_phi_eval = 0.d0
If (present(phi_eval_deg)) Then
  my_phi_eval = phi_eval_deg*pi/180.d0
Endif

! Always start by guessing xpt from Bpmin in gfile bdry
npts_bdry = Size(g_bdry,2)
Allocate(rtmp(npts_bdry))
Allocate(ztmp(npts_bdry))

! Throw away zero points (should not be any here anyway)
icount = 0
Do i = 1,npts_bdry
  If (g_bdry(1,i) > 1.d-4) Then
    icount = icount + 1
    rtmp(icount) = g_bdry(1,i)
    ztmp(icount) = g_bdry(2,i)
  Endif
Enddo
If (icount == 0) Then
  Write(*,*) 'Error 1 in find_xpt_jdl'
  Stop
Endif

! set search area based on grid size
dx1_grid = (g_r(g_mw) - g_r(1))*.15  ! These seem to work ok, really should make sure 
dx2_grid = (g_z(g_mh) - g_z(1))*.15  ! boundary not exceeded by these guesses
If (present(dx)) Then
  dx1_grid=dx
  dx2_grid=dx
Endif

Allocate(Bout(icount,3))

! Get minimum poloidal field on boundary (initial guess, still just gfile here)
Call bfield_geq_bicub(rtmp(1:icount),ztmp(1:icount),icount,Bout,ierr)
Bout(:,1) = Sqrt( Bout(:,1)**2 + Bout(:,2)**2 )
bpx = Minval(Bout(:,1))
ix  = Minloc(Bout(:,1),1)
rx  = rtmp(ix)
zx  = ztmp(ix)
If ( .NOT. quiet) Then
  Write(*,'(a,e12.3,2f12.3)') ' Very rough first X-point. [Bp,R,Z] = ',bpx,rx,zx
Endif
Deallocate(rtmp,ztmp,Bout)

! Evaluate on grid around initial guess to find min Bp
Allocate(ztmp(n1),Bout(n1,3))
Allocate(Bout_tmp(n1,3),phi_tmp(n1))
If (refine) Then
  err = bpx
  de = dx1_grid
  bp = 0.d0
  rg = 0.d0
  zg = 0.d0
  niter = 1
  Do While ((err .gt. tol) .AND. (niter .lt. niter_max))
    rt = rlinspace(rx-0.5d0*de,rx+0.5d0*de,n1)
    zt = rlinspace(zx-0.5d0*de,zx+0.5d0*de,n1)

    Do i = 1,n1 
      ztmp(:) = zt(i)            
      If (bfield_method == 0 ) Then      ! g only
        Call bfield_geq_bicub(rt,ztmp,n1,Bout,ierr)
#ifdef HAVE_M3DC1        
      Elseif ( bfield_method ==  3) Then ! g + m3dc1 
        phi_tmp(:) = my_phi_eval
        Call bfield_geq_bicub(rt,ztmp,n1,Bout,ierr)     
        Call bfield_m3dc1(rt,phi_tmp,ztmp,n1,Bout_tmp,ierr)
        Bout = Bout + Bout_tmp
      Elseif (bfield_method == 4) Then  ! m3dc1 total field
        phi_tmp(:) = my_phi_eval
        Call bfield_m3dc1(rt,phi_tmp,ztmp,n1,Bout,ierr)
#endif        
      Else
        Write(*,*) 'Bad value for bfield_method in find_xpt_jdl'
        Stop
      Endif
      
      rg(i,:) = rt
      zg(i,:) = ztmp
      bp(i,:) = Sqrt(Bout(:,1)**2 + Bout(:,2)**2)
    Enddo
    
    ixjx = Minloc(bp)
    bpx = bp(ixjx(1),ixjx(2))
    de = Sqrt( (rg(ixjx(1),ixjx(2))-rx)**2 + (zg(ixjx(1),ixjx(2))-zx)**2 )
    err = bpx

    rx = rg(ixjx(1),ixjx(2))
    zx = zg(ixjx(1),ixjx(2))

    niter = niter + 1
    If ( niter >= niter_max ) Then
        Write(*,*) 'Warning!! : niter_max exceeded for 1st x-point'
    Endif
  Enddo
  If ( .NOT. quiet ) Then
    Write(*,'(a,e12.3,2f12.3)') ' 1st X-point. [Bp,R,Z] = ',bpx,rx,zx
  Endif
Endif
Deallocate(Ztmp,Bout,Bout_tmp,phi_tmp)

! Find second xpoint
rx2 = 0.d0
zx2 = 0.d0
If (second) Then      
  ! Guess that config is approximately up-down symmetric
  rx2 = rx
  zx2 = -zx
  Allocate(Bout(1,3))
  Call bfield_geq_bicub((/rx2/),(/zx2/),1,Bout,ierr)
  bpx2 = Sqrt( Bout(1,1)**2 + Bout(1,2)**2 )
  err = bpx2
  Deallocate(Bout)
  Allocate(ztmp(n1),Bout(n1,3))
  Allocate(Bout_tmp(n1,3),phi_tmp(n1))  
  If (refine) Then
    de = dx2_grid
    bp = 0.d0
    rg = 0.d0
    zg = 0.d0
    niter = 1

    Do While ((err .gt. tol) .AND. (niter .lt. niter_max))
      rt = rlinspace(rx2-0.5d0*de,rx2+0.5d0*de,n1)
      zt = rlinspace(zx2-0.5d0*de,zx2+0.5d0*de,n1)
      
      Do i = 1,n1 
        ztmp(:) = zt(i)
        If (bfield_method == 0 ) Then      ! g only
          Call bfield_geq_bicub(rt,ztmp,n1,Bout,ierr)
#ifdef HAVE_M3DC1
        Elseif ( bfield_method ==  3) Then ! g + m3dc1 
          phi_tmp(:) = my_phi_eval
          Call bfield_geq_bicub(rt,ztmp,n1,Bout,ierr)     
          Call bfield_m3dc1(rt,phi_tmp,ztmp,n1,Bout_tmp,ierr)
          Bout = Bout + Bout_tmp
        Elseif (bfield_method == 4) Then  ! m3dc1 total field
          phi_tmp(:) = my_phi_eval
          Call bfield_m3dc1(rt,phi_tmp,ztmp,n1,Bout,ierr)
#endif          
        Else
          Write(*,*) 'Bad value for bfield_method in find_xpt_jdl'
          Stop
        Endif
        
        rg(i,:) = rt
        zg(i,:) = ztmp
        bp(i,:) = Sqrt(Bout(:,1)**2 + Bout(:,2)**2)
      Enddo
      
      ixjx = Minloc(bp)
      bpx2 = bp(ixjx(1),ixjx(2))
      de = Sqrt( (rg(ixjx(1),ixjx(2))-rx2)**2 + (zg(ixjx(1),ixjx(2))-zx2)**2 )
      err = bpx2
      
      rx2 = rg(ixjx(1),ixjx(2))
      zx2 = zg(ixjx(1),ixjx(2))

      niter = niter + 1
      If ( niter .ge. niter_max ) Then
        Write(*,*) 'Warning!! : niter_max exceeded for 2nd x-point'
      Endif
    Enddo
    If ( .NOT. quiet ) Then
      Write(*,'(a,e12.3,2f12.3)') ' 2nd X-point. [Bp,R,Z] = ',bpx2,rx2,zx2
    Endif
  Endif
  Deallocate(Ztmp,Bout,Bout_tmp,phi_tmp)  
Endif

End Subroutine find_xpt_jdl



End Module util_routines
