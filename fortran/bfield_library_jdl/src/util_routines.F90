!-----------------------------------------------------------------------------
!
!   Routines/modules with utility functions that are (ideally) not specific
!    to one field method
!
!   Contains:
!     Module util_routines
!       Function num_lines_file
!       Subroutine calc_sep
!       Subroutine find_xpt_jdl
!       Subroutine get_psi_2d
!       Subroutine get_psiN_2d
!
!-----------------------------------------------------------------------------
Module util_routines
  Use bfield, Only : bfield_type
  Implicit None
  Private

  Public :: num_lines_file
  Public :: calc_sep
  Public :: find_xpt_jdl
  Public :: get_psi_2d
  Public :: get_psiN_2d
  
  Interface get_psi_2d
    Module Procedure get_psi_2d_scalar
    Module Procedure get_psi_2d_array
  End Interface get_psi_2d

  Interface get_psiN_2d
    Module Procedure get_psiN_2d_scalar
    Module Procedure get_psiN_2d_array
  End Interface Get_psiN_2d  

Contains

  

  !-----------------------------------------------------------------------------
  !+ Scalar procedure for get_psi_2d
  !-----------------------------------------------------------------------------  
  Function get_psi_2d_scalar(bfield,r,z,ierr) &
       Result(psi)
    Use kind_mod, Only : int32, real64
    Implicit None
    Type(bfield_type), Intent(In) :: bfield
    Real(real64), Intent(In) :: r,z
    Integer(int32), Intent(Out) :: ierr
    Real(real64) :: psi, psi_tmp(1)
    psi_tmp = get_psi_2d_array(bfield,[r],[z],1,ierr)
    psi = psi_tmp(1)
  End Function get_psi_2d_scalar


  !-----------------------------------------------------------------------------
  !+ Array procedure for get_psi_2d
  !-----------------------------------------------------------------------------  
  Function get_psi_2d_array(bfield,r,z,Npts,ierr) &
       Result(psi)
    Use kind_mod, Only: int32, real64  
    Use g3d_module, Only: get_psi_bicub
#ifdef HAVE_M3DC1
    Use m3dc1_routines_mod, Only: calc_psi_m3dc1_2d
#endif
    Implicit None
    Type(bfield_type), Intent(In) :: bfield
    Real(real64), Dimension(Npts), Intent(In) :: r,z
    Integer(int32), Intent(In) :: Npts
    Integer(int32), Intent(Out) :: ierr
    Real(real64), Dimension(Npts) :: psi
    Integer(int32) :: ierr_tmp

    If ( (bfield%method == 0)  .OR. &    ! g only
         (bfield%method == 1)  .OR. &    ! g+rmp coils
         (bfield%method == 2)  .OR. &    ! g+screening
         (bfield%method == 3)  .OR. &    ! g + m3dc1
         (bfield%method == 7)  .OR. &    ! ipec
         (bfield%method == 8)  .OR. &
         (bfield%method == 9)  .OR. &
         (bfield%method == 10) .OR. &    ! xpand
         (bfield%method == 11) &         
      ) Then    
      Call get_psi_bicub(bfield%g,r,z,Npts,psi,ierr_tmp)
#ifdef HAVE_M3DC1
    Elseif ( (bfield%method == 4) .OR. & ! m3dc1 total field
         (bfield%method == 5) ) Then ! m3dc1 total field (AS only)  
      Call calc_psi_m3dc1_2d(r,z,Npts,psi,ierr_tmp)
#endif
    Else
      Write(*,*) 'Bad value for bfield%method in get_psi_2d_array'
      Write(*,*) 'bfield%method is',bfield%method
      Stop
    Endif
    ierr = 0
    If (ierr_tmp .ne. 0) ierr = 1
  End Function get_psi_2d_array

  !-----------------------------------------------------------------------------
  !+ Scalar procedure for get_psiN_2d_scalar
  !-----------------------------------------------------------------------------  
  Function get_psiN_2d_scalar(bfield,r,z,ierr) &
       Result(psiN)
    Use kind_mod, Only : int32, real64
    Implicit None
    Type(bfield_type), Intent(In) :: bfield
    Real(real64), Intent(In) :: r,z
    Integer(int32), Intent(Out) :: ierr
    Real(real64) :: psiN, psiN_tmp(1)
    psiN_tmp = get_psiN_2d_array(bfield,[r],[z],1,ierr)
    psiN = psiN_tmp(1)
  End Function get_psiN_2d_scalar  
  
  !-----------------------------------------------------------------------------
  !+ Array procedure for get_psiN_2d
  !-----------------------------------------------------------------------------    
  Function get_psiN_2d_array(bfield,r,z,Npts,ierr) &
       Result(psiN)
    Use kind_mod, Only: int32, real64  
    Use g3d_module, Only: get_psiN_bicub
#ifdef HAVE_M3DC1
    Use m3dc1_routines_mod, Only: calc_psiN_m3dc1_2d
#endif
    Implicit None
    Type(bfield_type), Intent(In) :: bfield
    Real(real64), Dimension(Npts), Intent(In) :: r,z
    Integer(int32), Intent(In) :: Npts
    Integer(int32), Intent(Out) :: ierr
    Real(real64), Dimension(Npts) :: psiN
    Integer(int32) :: ierr_tmp

    If ( (bfield%method == 0)  .OR. &    ! g only
         (bfield%method == 1)  .OR. &    ! g+rmp coils
         (bfield%method == 2)  .OR. &    ! g+screening
         (bfield%method == 3)  .OR. &    ! g + m3dc1
         (bfield%method == 7)  .OR. &    ! ipec
         (bfield%method == 8)  .OR. &
         (bfield%method == 9)  .OR. &
         (bfield%method == 10) .OR. &   ! xpand      
         (bfield%method == 11) &         
      ) Then          
      Call get_psiN_bicub(bfield%g,r,z,Npts,psiN,ierr_tmp)
#ifdef HAVE_M3DC1
    Elseif ( &
         (bfield%method == 4) .OR. & ! m3dc1 total field
         (bfield%method == 5) &      ! m3dc1 total field (AS only)  
         ) Then 
      Call calc_psiN_m3dc1_2d(r,z,Npts,psiN,ierr_tmp)
#endif          
    Else
      Write(*,*) 'Bad value for bfield%method in get_psiN_2d_array'
      Write(*,*) 'bfield%method is',bfield%method
      Stop
    Endif
    ierr = 0
    If (ierr_tmp .ne. 0) ierr = 1
  End Function get_psiN_2d_array


  !-----------------------------------------------------------------------------
  !+ Returns the number of lines in a file
  !-----------------------------------------------------------------------------
  Function num_lines_file(fname) &
       Result(numl)
    Use kind_mod, Only :int32
    Implicit None
    Character(Len=*), Intent(In) :: fname
    Integer(int32) :: numl, iocheck
    Integer(int32), Parameter :: max_numl = 1000000
    !- End of header -------------------------------------------------------------
    Open(99,file = fname,IOSTAT=iocheck,STATUS="OLD")
    If (iocheck /= 0) Then
      Write(*,*) 'Error opening file: ', fname
      Stop 'Exiting: I/O error in function num_lines_file'
    Endif
    numl = 0
    Do
      Read(99,*,IOSTAT=iocheck)
      If (iocheck /= 0) Exit
      numl = numl + 1
      If ((numl .GE. Huge(numl) - 1) .OR. (numl .GE. max_numl)) Then
        Stop 'Error: Exceeded maximum numl in num_lines_file'
      Endif
    Enddo
    Close(99)
  End Function num_lines_file

  !-----------------------------------------------------------------------------
  !+ Calculates approximate curve(s) for the sepatrix
  !-----------------------------------------------------------------------------
  Subroutine calc_sep(bfield,second_sep,fname_out)
    !
    ! Description:
    !  Calculates the separatrix curve(s) 
    !
    !  bfield%method must be set first, and a gfile must also have been read
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
    Use kind_mod, Only: int32, real64
    Use fieldline_follow_mod, Only: follow_fieldlines_rzphi
    Use math_geo_module, Only : inside_poly
    Use phys_const, Only: pi
    Implicit None
    !Input/output
    Type(bfield_type), Intent(In) :: bfield
    Logical, Intent(In) :: second_sep
    Character(Len=*),Intent(In) :: fname_out
    ! Local Variables
    Real(real64) :: rx,zx,rx2,zx2
    Integer(int32) :: idir,maxtries,ii,i,nlim,nbox,in1,in2,nsep1,nsep0,nseps,isep
    Real(real64), Allocatable, Dimension(:) :: rsep0,zsep0,rsep0_tmp,zsep0_tmp,lim_r,lim_z
    Real(real64), Allocatable, Dimension(:) :: rsep1,zsep1
    Real(real64), Allocatable, Dimension(:) ::  box_r,box_z
    ! Fl folllowing vars
    Real(real64), Dimension(1) :: rstart,zstart,phistart
    Integer(int32), Dimension(1) :: fl_ierr,ilg
    Real(real64), Allocatable, Dimension(:,:) :: fl_r,fl_z,fl_p
    Real(real64) :: dphi_fl
    Integer(int32) :: nsteps_fl
    ! Local Parameters
    Real(real64), Parameter :: dx = 1.d-3
    Integer(real64), Parameter :: nsep_div = 1000  ! total number of sep points from fl follow divided by this to output
    !- End of header -------------------------------------------------------------

    Open(99,file=Trim(Adjustl(fname_out)))
    If (second_sep) Then
      nseps = 3
    Else
      nseps = 1
    Endif
    Write(99,*) nseps

    ! Set up boundaries to search
    ! 1) gfile lim
    Allocate(lim_r(bfield%g%limitr))
    Allocate(lim_z(bfield%g%limitr))
    ! Throw away zero points 
    nlim = 0
    Do i = 1,bfield%g%limitr
      If (bfield%g%lim(1,i) > 1.d-4) Then
        nlim = nlim + 1
        lim_r(nlim) = bfield%g%lim(1,i)
        lim_z(nlim) = bfield%g%lim(2,i)
      Endif
    Enddo
    If (nlim == 0) Then
      Write(*,*) 'Error: could not create limiter bdry in calc_sep'
      Stop
    Endif
    ! 2) box based on gfile domain size
    nbox = 4
    Allocate(box_r(nbox),box_z(nbox))
    box_r = (/bfield%g%r(1),bfield%g%r(bfield%g%mw),bfield%g%r(bfield%g%mw),bfield%g%r(1)/)
    box_z = (/bfield%g%z(1),bfield%g%z(1),bfield%g%z(bfield%g%mh),bfield%g%z(bfield%g%mh)/)


    Call find_xpt_jdl(bfield,.true.,.true.,1.d-8,.false.,rx,zx,rx2,zx2)

    dphi_fl = bfield%g%ip_sign*0.1d0*pi/180.d0
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
          Call follow_fieldlines_rzphi(bfield,rstart,zstart,phistart,1,dphi_fl*idir,nsteps_fl,fl_r,fl_z,fl_p,fl_ierr,ilg)
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
      Do i = 1,nsep0-1,Int((nsep0+nsep1)/nsep_div,int32)
        Write(99,*) rsep0(i),zsep0(i)
      Enddo
      Do i = 1,nsep1,Int((nsep0+nsep1)/nsep_div,int32)
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
  Subroutine find_xpt_jdl(bfield,second,refine,tol,quiet,rx,zx,rx2,zx2,dx)
    ! tol is magnitude of bp at xpt
    ! phi_eval_deg only used for m3dc1 fields
    ! Both cases require a gfile field for initial guess!!!
    Use kind_mod, Only: real64, int32
    Use math_geo_module, Only: rlinspace
#ifdef HAVE_M3DC1
    Use m3dc1_routines_mod, Only: bfield_m3dc1, bfield_m3dc1_2d
#endif
    Use g3d_module, Only: bfield_geq_bicub
    Implicit None
    Type(bfield_type), Intent(In) :: bfield
    Logical, Intent(in) :: second, refine, quiet
    Real(real64), Intent(in) :: tol
    Real(real64), Intent(in), Optional :: dx
    Real(real64), Intent(out) :: rx, zx, rx2, zx2

    Integer(int32), Parameter :: niter_max = 15
    Integer(int32), Parameter :: n1 = 100  ! grid dimension (square)

    Real(real64), Allocatable :: rtmp(:), ztmp(:),Bout(:,:)
    Real(real64) :: bp(n1,n1), rg(n1,n1), zg(n1,n1), rt(n1), zt(n1)
    Real(real64) :: bpx, err, bpx2, dx1_grid, dx2_grid
    Real(real64) :: Rmin_eval, Rmax_eval, Zmin_eval, Zmax_eval, dx1_save, dx2_save
    Integer(int32) :: icount, i, npts_bdry, ierr, ix,ixjx(2), niter
    ! Local parameters               
    !- End of header -------------------------------------------------------------

    ! Always start by guessing xpt from Bpmin in gfile bdry
    npts_bdry = Size(bfield%g%bdry,2)
    Allocate(rtmp(npts_bdry))
    Allocate(ztmp(npts_bdry))

    ! Throw away zero points (should not be any here anyway)
    icount = 0
    Do i = 1,npts_bdry
      If (bfield%g%bdry(1,i) > 1.d-4) Then
        icount = icount + 1
        rtmp(icount) = bfield%g%bdry(1,i)
        ztmp(icount) = bfield%g%bdry(2,i)
      Endif
    Enddo
    If (icount == 0) Then
      Write(*,*) 'Error 1 in find_xpt_jdl'
      Stop
    Endif

    ! set search area based on grid size
    dx1_grid = (bfield%g%r(bfield%g%mw) - bfield%g%r(1))*.05  ! These seem to work ok, really should make sure 
    dx2_grid = (bfield%g%z(bfield%g%mh) - bfield%g%z(1))*.05  ! boundary not exceeded by these guesses
    If (present(dx)) Then
      dx1_grid=dx
      dx2_grid=dx
    Endif

    dx1_save = dx1_grid
    dx2_save = dx2_grid

    Rmin_eval = bfield%g%r(3) + 1.d-3
    Zmin_eval = bfield%g%z(3) + 1.d-3
    Rmax_eval = bfield%g%r(bfield%g%mw-2) - 1.d-3
    Zmax_eval = bfield%g%z(bfield%g%mh-2) - 1.d-3

    Allocate(Bout(icount,3))

    ! Get minimum poloidal field on boundary (initial guess, still just gfile here)
    Call bfield_geq_bicub(bfield%g,rtmp(1:icount),ztmp(1:icount),icount,Bout,ierr)
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
    If (refine) Then
      err = bpx
      bp = 0.d0
      rg = 0.d0
      zg = 0.d0
      niter = 1
      Do While ((err .gt. tol) .AND. (niter .lt. niter_max))
        rt = rlinspace(Max(Rmin_eval,rx-dx1_grid),Min(Rmax_eval,rx+dx1_grid),n1)
        zt = rlinspace(Max(Zmin_eval,zx-dx2_grid),Min(Zmax_eval,zx+dx2_grid),n1)
     
        Do i = 1,n1 
          ztmp(:) = zt(i)
          If ( &
               bfield%method == 0  .OR. & ! g only
               bfield%method == 1  .OR. & ! g+rmp coils
               bfield%method == 2  .OR. & ! g+screening
               bfield%method == 3  .OR. & ! g+m3dc1
               bfield%method == 7  .OR. & ! ipec eq only
               bfield%method == 8  .OR. & ! ipec vac
               bfield%method == 9  .OR. & ! ipec pert
               bfield%method == 10 .OR. & ! xpand pert
               bfield%method == 11      & ! xpand vac
               ) Then      
            Call bfield_geq_bicub(bfield%g,rt,ztmp,n1,Bout,ierr)
#ifdef HAVE_M3DC1
          Elseif ( &
               bfield%method == 4 .OR. & ! m3dc1 total field
               bfield%method == 5      & ! m3dc1 total field (AS only)
               ) Then 
            Call bfield_m3dc1_2d(rt,ztmp,n1,Bout,ierr)
#endif          
          Else
            Write(*,*) 'Bad value for bfield%method in find_xpt_jdl'
            Write(*,*) 'bfield%method is',bfield%method
            Stop
          Endif

          rg(i,:) = rt
          zg(i,:) = ztmp
          bp(i,:) = Sqrt(Bout(:,1)**2 + Bout(:,2)**2)
        Enddo

        ixjx = Minloc(bp)
        bpx = bp(ixjx(1),ixjx(2))
!        de = Sqrt( (rg(ixjx(1),ixjx(2))-rx)**2 + (zg(ixjx(1),ixjx(2))-zx)**2 )
        dx1_grid = Abs(rg(ixjx(1),ixjx(2))-rx)
        dx2_grid = Abs(zg(ixjx(1),ixjx(2))-zx)
        err = bpx

        rx = rg(ixjx(1),ixjx(2))
        zx = zg(ixjx(1),ixjx(2))

        niter = niter + 1
        If ( niter >= niter_max ) Then
          Write(*,*) 'Warning!! : niter_max exceeded for 1st x-point'
        Endif
      Enddo
      If ( .NOT. quiet ) Then
        Write(*,'(a,e12.3,2f12.5)') ' 1st X-point. [Bp,R,Z] = ',bpx,rx,zx
      Endif
    Endif
    Deallocate(Ztmp,Bout)

    ! Find second xpoint
    rx2 = 0.d0
    zx2 = 0.d0
    dx1_grid = dx1_save
    dx2_grid = dx2_save
    If (second) Then      
      ! Guess that config is approximately up-down symmetric
      rx2 = rx
      zx2 = -zx
      Allocate(Bout(1,3))
      Call bfield_geq_bicub(bfield%g,(/rx2/),(/zx2/),1,Bout,ierr)
      bpx2 = Sqrt( Bout(1,1)**2 + Bout(1,2)**2 )
      err = bpx2
      Deallocate(Bout)
      Allocate(ztmp(n1),Bout(n1,3))
      If (refine) Then
!        de = dx2_grid
        bp = 0.d0
        rg = 0.d0
        zg = 0.d0
        niter = 1

        Do While ((err .gt. tol) .AND. (niter .lt. niter_max))
!          rt = rlinspace(rx2-0.5d0*de,rx2+0.5d0*de,n1)
!          zt = rlinspace(zx2-0.5d0*de,zx2+0.5d0*de,n1)
          rt = rlinspace(Max(Rmin_eval,rx2-dx1_grid),Min(Rmax_eval,rx2+dx1_grid),n1)
          zt = rlinspace(Max(Zmin_eval,zx2-dx2_grid),Min(Zmax_eval,zx2+dx2_grid),n1)

          
          Do i = 1,n1 
            ztmp(:) = zt(i)

            If ( &
                 bfield%method == 0  .OR. & ! g only
                 bfield%method == 1  .OR. & ! g+rmp coils
                 bfield%method == 2  .OR. & ! g+screening
                 bfield%method == 3  .OR. & ! g+m3dc1
                 bfield%method == 7  .OR. & ! ipec eq only
                 bfield%method == 8  .OR. & ! ipec vac
                 bfield%method == 9  .OR. & ! ipec pert
                 bfield%method == 10 .OR. & ! xpand pert
                 bfield%method == 11      & ! xpand vac
                 ) Then      
              Call bfield_geq_bicub(bfield%g,rt,ztmp,n1,Bout,ierr)
#ifdef HAVE_M3DC1
            Elseif ( &
                 bfield%method == 4 .OR. & ! m3dc1 total field
                 bfield%method == 5      & ! m3dc1 total field (AS only)
                 ) Then 
              Call bfield_m3dc1_2d(rt,ztmp,n1,Bout,ierr)
#endif          
            Else
              Write(*,*) 'Bad value for bfield%method in find_xpt_jdl'
              Write(*,*) 'bfield%method is',bfield%method
              Stop
            Endif

            rg(i,:) = rt
            zg(i,:) = ztmp
            bp(i,:) = Sqrt(Bout(:,1)**2 + Bout(:,2)**2)
          Enddo

          ixjx = Minloc(bp)
          bpx2 = bp(ixjx(1),ixjx(2))
!          de = Sqrt( (rg(ixjx(1),ixjx(2))-rx2)**2 + (zg(ixjx(1),ixjx(2))-zx2)**2 )
          dx1_grid = Abs(rg(ixjx(1),ixjx(2))-rx2)
          dx2_grid = Abs(zg(ixjx(1),ixjx(2))-zx2)

          err = bpx2

          rx2 = rg(ixjx(1),ixjx(2))
          zx2 = zg(ixjx(1),ixjx(2))

          niter = niter + 1
          If ( niter .ge. niter_max ) Then
            Write(*,*) 'Warning!! : niter_max exceeded for 2nd x-point'
          Endif
        Enddo
        If ( .NOT. quiet ) Then
          Write(*,'(a,e12.3,2f12.5)') ' 2nd X-point. [Bp,R,Z] = ',bpx2,rx2,zx2
        Endif
      Endif
      Deallocate(Ztmp,Bout)  
    Endif

  End Subroutine find_xpt_jdl



End Module util_routines
