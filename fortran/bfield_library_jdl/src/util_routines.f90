!-----------------------------------------------------------------------------
!
!   Routines/modules with utility functions that are (ideally) not specific
!    to one field method
!
!   Contains:
!     Subroutine calc_sep
!   
!-----------------------------------------------------------------------------
Module util_routines
Implicit None
Contains

!-----------------------------------------------------------------------------
!+ Calculates and outputs approximate curve(s) for the sepatrix
!-----------------------------------------------------------------------------
Subroutine calc_sep
!
! Description:
!  Calculates the separatrix curve(s) and outputs to file(s)
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
Use fieldline_follow_mod, Only: bfield_method
Implicit None
!Input/output
! Local Variables
Real(rknd) :: rx,zx,rx2,zx2
Integer(iknd) :: idir
! Local Parameters
Real(rknd), Parameter :: dx = 1.d-3
!- End of header -------------------------------------------------------------


Call find_xpt_jdl(.true.,.true.,1.d-8,.false.,rx,zx,rx2,zx2)
Do idir=-1,1,2
  Write(*,*) 'Working on direction ',idir
Enddo

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
Use m3dc1_routines_mod, Only: bfield_m3dc1
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
Real(rknd),parameter :: pi = 3.14159265358979323846_rknd
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
      Elseif ( bfield_method ==  3) Then ! g + m3dc1 
        phi_tmp(:) = my_phi_eval
        Call bfield_geq_bicub(rt,ztmp,n1,Bout,ierr)     
        Call bfield_m3dc1(rt,phi_tmp,ztmp,n1,Bout_tmp,ierr)
        Bout = Bout + Bout_tmp
      Elseif (bfield_method == 4) Then  ! m3dc1 total field
        phi_tmp(:) = my_phi_eval
        Call bfield_m3dc1(rt,phi_tmp,ztmp,n1,Bout,ierr)
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
        Elseif ( bfield_method ==  3) Then ! g + m3dc1 
          phi_tmp(:) = my_phi_eval
          Call bfield_geq_bicub(rt,ztmp,n1,Bout,ierr)     
          Call bfield_m3dc1(rt,phi_tmp,ztmp,n1,Bout_tmp,ierr)
          Bout = Bout + Bout_tmp
        Elseif (bfield_method == 4) Then  ! m3dc1 total field
          phi_tmp(:) = my_phi_eval
          Call bfield_m3dc1(rt,phi_tmp,ztmp,n1,Bout,ierr)
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
