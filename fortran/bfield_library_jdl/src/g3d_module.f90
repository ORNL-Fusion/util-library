!-----------------------------------------------------------------------------
!
!   Routines/modules related to g3d/efit
!    --> Mostly ported from Canik's idl routines
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!+ Module for variables read from gfile
!-----------------------------------------------------------------------------
Module gfile_var_pass
!
! Description:
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     04/14/2011   Original Code.  JDL
! 
! Author(s): J.D. Lore 4/14/2011
!
! Modules used:
Use kind_mod                ! Import rknd, iknd specifications

Implicit none
Integer, Parameter :: g_bspl_ord = 3
Integer(iknd):: &
  g_mw,g_mh,g_nbdry,g_limitr
Character(len=6) :: g_ecase
Real(rknd) :: g_xdim,g_zdim,g_rzero,g_rgrid1,g_zmid, &
  g_rmaxis,g_zmaxis,g_ssimag,g_ssibry,g_bcentr,g_cpasma, &
  g_dr,g_dz,g_ip_sign
Real(rknd), Allocatable :: g_fpol(:), &
  g_pres(:), g_ffprim(:), g_pprime(:), g_psirz(:,:), &
  g_qpsi(:), g_bdry(:,:),g_lim(:,:),g_r(:),g_z(:),g_pn(:), &
  g_bicub_coeffs(:,:,:),g_pnknot(:),g_fpol_bscoef(:)    

End module gfile_var_pass



!-----------------------------------------------------------------------------
!+ Module containing routines for gfile psi, B interpolation
!-----------------------------------------------------------------------------
Module g3d_module
Implicit None
Private :: psi_bi, dsdr_bi, dsdz_bi

Contains

!-----------------------------------------------------------------------------
!+ returns Bcyl from gfile bfield at R,Z
!-----------------------------------------------------------------------------
Subroutine bfield_geq_bicub(R1,Z1,Npts,Bout,ierr,verbose)
! Description: 
!  
! Output:
!   Bout = (:,[Br,Bz,Bt])
!
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/14/2011  Ported from Canik's idl routines.  JDL
! 
! Author(s): J.D. Lore 04/14/2011 
!
! Modules used:
Use kind_mod, Only : iknd, rknd
Use gfile_var_pass
Use bspline, Only : dbsval
Implicit None

! Input/output                      !See above for descriptions
Integer(iknd),Intent(in) :: Npts
Real(rknd),Dimension(Npts),Intent(in)  :: R1,Z1
Real(rknd),Dimension(Npts,3),Intent(out) :: Bout(Npts,3)
Integer(iknd),Intent(out) :: ierr
Logical, Optional :: verbose
! Local Scalars
Integer(iknd) :: ir,iz,index,ii
Real(rknd) :: dir,diz
Real(rknd) :: psi1,dsdr1,dsdz1,br1,bz1,psiN
! Local parameters               
Real(rknd), Parameter :: TWO = 2._rknd, THREE = 3._rknd  

Real(rknd) :: fpol,bt1
Logical :: myverbose 
!- End of header -------------------------------------------------------------
myverbose = .true.
If (present(verbose)) Then
  myverbose = verbose
Endif

!If (.not. Allocated(g_r) ) Then
!  Write(*,*) 'G VARIABLES NOT ALLOCATED, EXITING FROM bfield_geq_bicub!'
!  Stop
!Endif

ierr = 0
Do ii = 1,Npts 

  ! update for vec.
  ir = Floor( (R1(ii)-g_r(1))/g_dr ) + 1
  iz = Floor( (Z1(ii)-g_z(1))/g_dz ) + 1

  ! Check for points off grid
  If ( (ir .le. 2) .or. (ir .ge. g_mw - 1) ) Then
    If (myverbose) &
         Write(*,'(3(a,f12.3),a)') 'bfield_geq: Point off grid in R: R = ',R1(ii), &
         '. [Rmin,Rmax] = [',g_r(1),',',g_r(g_mw),']'
    ierr = 1
    Bout = 0.d0
    return
  Endif
  If ( (iz .le. 1) .or. (iz .ge. g_mh - 1) ) Then
    If (myverbose) &
         Write(*,'(3(a,f12.3),a)') 'bfield_geq: Point off grid in Z: Z = ',Z1(ii), &
         '. [Zmin,Zmax] = [',g_z(1),',',g_z(g_mh),']'
    ierr = 1
    Bout = 0.d0
    return
  Endif

  dir = (R1(ii) - g_r(ir))/g_dr
  diz = (Z1(ii) - g_z(iz))/g_dz

  index = iz + g_mh*(ir-1)
  psi1 = psi_bi(index,dir,diz)
  dsdr1 = dsdr_bi(index,dir,diz)
  dsdz1 = dsdz_bi(index,dir,diz)
  
  br1 = -dsdz1/R1(ii)
  bz1 =  dsdr1/R1(ii)
  
  psiN = (psi1*g_ip_sign - g_ssimag)/(g_ssibry-g_ssimag)
    
  ! Toroidal field
  If (psiN .ge. 0._rknd .AND. psiN .le. 1._rknd) Then 
    fpol = dbsval(psiN,g_bspl_ord,g_pnknot,g_mw,g_fpol_bscoef)
    bt1 = fpol/R1(ii)
  Else
    bt1 = g_bcentr*g_rzero/R1(ii)
  Endif
    
  Bout(ii,1) = br1
  Bout(ii,2) = bz1
  Bout(ii,3) = bt1
Enddo

End Subroutine bfield_geq_bicub

!-----------------------------------------------------------------------------
!+ reads a g file
!-----------------------------------------------------------------------------
Subroutine readg_g3d(filename)
!
! Description: 
!  Abbreviated function to read gfiles.  Assumes formatted file
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/12/2011  Ported from Canik's idl routines.  JDL
!
! Author(s): J.D. Lore - 04/12/2011
!
! Modules used:
Use kind_mod                ! Import rknd, iknd specifications
Use gfile_var_pass
Use bspline
Implicit None

! Input/output                      !See above for descriptions
Character(len=100) :: filename

! Local scalars
Integer(iknd) :: iocheck,idum,i,j
Real(rknd) :: xdum
Logical,Parameter :: DIAGNO = .false.  ! screen output of read variables

! Local arrays (1D)
Character(len=100) :: sjunk

Logical :: debug = .false.
! Allocatable arrays

!- End of header -------------------------------------------------------------

Write(*,*) 'Reading gfile: ',filename

Open(UNIT=99,FILE=filename,STATUS="old",IOSTAT=iocheck)
If ( iocheck /= 0 ) Then
    Write(*,*) 'Error opening gfile: ', filename
    Stop 'Exiting: I/O Error in subroutine readg_g3d'
Endif
If (debug) Write(*,*) 'Debugging readg_g3d'
Read(99,'(a8,a42,i3,2i4)') g_ecase,sjunk,idum,g_mw,g_mh
If (debug) Then
  Write(*,*) 'g_ecase',g_ecase
  Write(*,*) 'g_mw,g_mh',g_mw,g_mh
Endif
Allocate(g_fpol(g_mw),g_pres(g_mw),g_ffprim(g_mw))
Allocate(g_pprime(g_mw),g_psirz(g_mw,g_mh),g_qpsi(g_mw))

read(99,'(5e16.9)') g_xdim,g_zdim,g_rzero,g_rgrid1,g_zmid
If (debug) Then
  Write(*,*) 'g_xdim',g_xdim
  Write(*,*) 'g_zdim',g_zdim
  Write(*,*) 'g_rzero',g_rzero
  Write(*,*) 'g_rgrid1',g_rgrid1
  Write(*,*) 'g_zmid',g_zmid
Endif
read(99,'(5e16.9)') g_rmaxis,g_zmaxis,g_ssimag,g_ssibry,g_bcentr
If (debug) Then
  Write(*,*) 'g_rmaxis',g_rmaxis
  Write(*,*) 'g_zmaxis',g_zmaxis
  Write(*,*) 'g_ssimag',g_ssimag
  Write(*,*) 'g_ssibry',g_ssibry
  Write(*,*) 'g_bcentr',g_bcentr
Endif
read(99,'(1e16.9)') g_cpasma
If (debug) Write(*,*) 'g_cpasma',g_cpasma
read(99,'(5e16.9)') xdum

read(99,'(5e16.9)') (g_fpol(i),i=1,g_mw)
If (debug) Write(*,*) 'g_fpol([1,g_mw])',g_fpol(1),g_fpol(g_mw)
read(99,'(5e16.9)') (g_pres(i),i=1,g_mw)
read(99,'(5e16.9)') (g_ffprim(i),i=1,g_mw)
read(99,'(5e16.9)') (g_pprime(i),i=1,g_mw)

read(99,'(5e16.9)') ((g_psirz(i,j),i=1,g_mw),j=1,g_mh)
read(99,'(5e16.9)') (g_qpsi(i),i=1,g_mw)
read(99,'(2i5)')    g_nbdry,g_limitr

Allocate(g_bdry(2,g_nbdry))
Allocate(g_lim(2,g_limitr))
read(99,'(5e16.9)') ((g_bdry(i,j),i=1,2),j=1,g_nbdry)
read(99,'(5e16.9)') ((g_lim(i,j),i=1,2),j=1,g_limitr)

Close(99)

If (DIAGNO .eqv. .true.) Then
  Write(*,*) '------------------------------------------------------'
  Write(*,*) '   DIAGNO is on in readg_g3d:'
  Write(*,*) 'g_ecase:',g_ecase
  Write(*,*) 'g_mw:',g_mw
  Write(*,*) 'g_mh:',g_mh
  Write(*,*) 'g_xdim',g_xdim
  Write(*,*) 'g_zdim',g_zdim
  Write(*,*) 'g_rzero',g_rzero
  Write(*,*) 'g_rgrid1',g_rgrid1
  Write(*,*) 'g_zmid',g_zmid
  Write(*,*) 'g_rmaxis',g_rmaxis
  Write(*,*) 'g_zmaxis',g_zmaxis
  Write(*,*) 'g_ssimag',g_ssimag
  Write(*,*) 'g_ssibry',g_ssibry
  Write(*,*) 'g_bcentr',g_bcentr
  Write(*,*) 'g_cpasma',g_cpasma
  Write(*,*) 'g_fpol:'
  Do i = 1,g_mw
    Write(*,*) g_fpol(i)
  Enddo
  Write(*,*) 'g_pres:'
  Do i = 1,g_mw
    Write(*,*) g_pres(i)
  Enddo
  Write(*,*) 'g_ffprim:'
  Do i = 1,g_mw
    Write(*,*) g_ffprim(i)
  Enddo
  Write(*,*) 'g_pprime:'
  Do i = 1,g_mw
    Write(*,*) g_pprime(i)
  Enddo
  Write(*,*) 'g_psirz:'
  Do i = 1,g_mw
    Write(*,*) g_psirz(i,1:g_mh)
  Enddo
  Write(*,*) 'g_qpsi:'
  Do i = 1,g_mw
    Write(*,*) g_qpsi(i)
  Enddo
  Write(*,*) 'g_nbdry',g_nbdry
  Write(*,*) 'g_limitr',g_limitr

  Write(*,*) '------------------------------------------------------'
Endif
 

!
! Postprocessing
!

g_dR = g_xdim/(g_mw-1)
g_dZ = g_zdim/(g_mh-1)

Allocate(g_r(g_mw),g_z(g_mh),g_pn(g_mw))

Do i=0,g_mw-1
  g_r(i+1) = g_rgrid1 + g_dR*i
  g_pn(i+1) = Real(i,rknd)/(g_mw-1)
Enddo

Do i=0,g_mh-1
  g_z(i+1) = g_zmid - 0.5_rknd*g_zdim + g_dZ*i
Enddo

g_ip_sign = -g_cpasma/dabs(g_cpasma)

Allocate(g_bicub_coeffs(g_mw*g_mh,4,4))
g_bicub_coeffs = get_psi_bicub_coeffs()

! B-Spline fit poloidal current function
Allocate(g_pnknot(g_bspl_ord+g_mw))
Allocate(g_fpol_bscoef(g_mw))
Call dbsnak(g_mw,g_pn,g_bspl_ord,g_pnknot)
Call dbsint(g_mw,g_pn,g_fpol,g_bspl_ord,g_pnknot,g_fpol_bscoef)

End Subroutine readg_g3d


!-----------------------------------------------------------------------------
!+ Returns array of coefficients for bicubic interpolation
!-----------------------------------------------------------------------------
Function get_psi_bicub_coeffs()  &
Result(psi_bicub_coeffs)
!
! Description: 
!  
! Function arguments:
! Output:
!
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/15/2011  Ported from Canik's idl routines.  JDL
! 
! Author(s): J.D. Lore 04/15/2011 
!
! Modules used:
Use kind_mod                  ! Import rknd, iknd specifications
Use gfile_var_pass, Only : &
 g_dr,g_dz,g_mw,g_mh,g_ip_sign,g_psirz
Use g3df_math_routines_mod, Only : &
! Imported routines
inversion_lu  
  
Implicit None

! Input/output                      !See above for descriptions
Real(rknd),dimension(g_mw*g_mh,4,4)   :: psi_bicub_coeffs

! Local Scalars
Integer(iknd) :: nr,nz,ir,iz,index,inv_err
! Local arrays 
Real(rknd),dimension(g_mw,g_mh) :: psi2d,dsdr,dsdz,d2sdrdz
Real(rknd),dimension(16,16) :: bicub_mat,bicub_mat_inv
Real(rknd),dimension(16) :: b,coeff
! Local parameters                 


!- End of header -------------------------------------------------------------

nr = g_mw
nz = g_mh

psi2d = g_ip_sign * g_psirz

dsdr = (cshift(psi2d,shift=1,dim=1) - cshift(psi2d,shift=-1,dim=1))/(2._rknd*g_dr)
dsdz = (cshift(psi2d,shift=1,dim=2) - cshift(psi2d,shift=-1,dim=2))/(2._rknd*g_dz)
d2sdrdz = (cshift(cshift(psi2d,shift=1,dim=1),shift=1,dim=2) & 
  - cshift(cshift(psi2d,shift=-1,dim=1),shift=1,dim=2) &
  - cshift(cshift(psi2d,shift=1,dim=1),shift=-1,dim=2) &
  + cshift(cshift(psi2d,shift=-1,dim=1),shift=-1,dim=2))/(4._rknd*g_dr*g_dz)

bicub_mat = get_bicub_mat()

Call Inversion_LU(bicub_mat,bicub_mat_inv,Size(bicub_mat,1),inv_err)
bicub_mat_inv = Transpose(bicub_mat_inv)

psi_bicub_coeffs = 0._rknd
Do ir = 1,nr-1
  Do iz = 1,nz-1
    index = iz + nz*(ir-1)
    b = (/psi2d(ir,iz),            psi2d(ir+1,iz),            psi2d(ir,iz+1),            psi2d(ir+1,iz+1),     &
          dsdr(ir,iz)*g_dr,        dsdr(ir+1,iz)*g_dr,        dsdr(ir,iz+1)*g_dr,        dsdr(ir+1,iz+1)*g_dr, &
          dsdz(ir,iz)*g_dz,        dsdz(ir+1,iz)*g_dz,        dsdz(ir,iz+1)*g_dz,        dsdz(ir+1,iz+1)*g_dz, &
          d2sdrdz(ir,iz)*g_dr*g_dz,d2sdrdz(ir+1,iz)*g_dr*g_dz,d2sdrdz(ir,iz+1)*g_dr*g_dz,d2sdrdz(ir+1,iz+1)*g_dr*g_dz/)
    coeff = Matmul(b,bicub_mat_inv)
    psi_bicub_coeffs(index,:,1) = coeff(1:4)
    psi_bicub_coeffs(index,:,2) = coeff(5:8)
    psi_bicub_coeffs(index,:,3) = coeff(9:12)
    psi_bicub_coeffs(index,:,4) = coeff(13:16)
  Enddo
Enddo

EndFunction get_psi_bicub_coeffs

!-----------------------------------------------------------------------------
!+ Returns bicubic interpolation matrix
!-----------------------------------------------------------------------------
Function get_bicub_mat()  &
Result(bicub_mat)
!
! Description: 
!  
! Function arguments:
! Output:
!
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/19/2011  Ported from Canik's idl routines.  JDL
! 
! Author(s): J.D. Lore 04/19/2011 
!
! Modules used:
Use kind_mod                  ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Real(rknd),dimension(16,16)        :: bicub_mat

! Local Parameters
Real(rknd), Parameter :: One   = 1._rknd, &
                         Two   = 2._rknd, &
                         Three = 3._rknd, &
                         Four  = 4._rknd, &
                         Six   = 6._rknd, &
                         Nine  = 9._rknd

!- End of header -------------------------------------------------------------

bicub_mat(:,:) = 0._rknd

! Function values at corners
bicub_mat(1,1) = One
bicub_mat(2,(/1,2,3,4/)) = (/One,One,One,One/)
bicub_mat(3,(/1,5,9,13/)) = (/One,One,One,One/)
bicub_mat(4,:) = One


! 1st derivatives at corners: x direction
bicub_mat(5,2) = One
bicub_mat(6,(/2,3,4/)) = (/One,Two,Three/)
bicub_mat(7,(/2,6,10,14/)) = (/One,One,One,One/)
bicub_mat(8,(/2,3,4/)) = (/One,Two,Three/)
bicub_mat(8,(/6,7,8/)) = (/One,Two,Three/)
bicub_mat(8,(/10,11,12/)) = (/One,Two,Three/)
bicub_mat(8,(/14,15,16/)) = (/One,Two,Three/)

! 1st derivatives at corners: y direction
bicub_mat(9,5) = One
bicub_mat(10,(/5,6,7,8/)) = (/One,One,One,One/)
bicub_mat(11,(/5,9,13/)) = (/One,Two,Three/)
bicub_mat(12,(/5,9,13/)) = (/One,Two,Three/)
bicub_mat(12,(/6,10,14/)) = (/One,Two,Three/)
bicub_mat(12,(/7,11,15/)) = (/One,Two,Three/)
bicub_mat(12,(/8,12,16/)) = (/One,Two,Three/)

! cross derivatives at corners
bicub_mat(13,6) = One
bicub_mat(14,(/6,7,8/)) = (/One,Two,Three/)
bicub_mat(15,(/6,10,14/)) = (/One,Two,Three/)
bicub_mat(16,(/6,10,14/)) = (/One,Two,Three/)
bicub_mat(16,(/7,11,15/)) = (/Two,Four,Six/)
bicub_mat(16,(/8,12,16/)) =(/Three,Six,Nine/)

EndFunction get_bicub_mat


!-----------------------------------------------------------------------------
!+ Find x-point(s) from gfile field
!-----------------------------------------------------------------------------
Subroutine find_xpt_jdl(second,refine,tol,quiet,rx,zx,rx2,zx2,phi_eval_deg,dx)
  ! tol is magnitude of bp at xpt
  ! phi_eval_deg only used for m3dc1 fields
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

!-----------------------------------------------------------------------------
!+ returns psi from gfile bfield at R,Z
!-----------------------------------------------------------------------------
Subroutine get_psi_bicub(R1,Z1,Npts,psiout,psiNout,ierr)
! Description: 
!  
! Output:
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/14/2011  Ported from Canik's idl routines.  JDL
! 
! Author(s): J.D. Lore 04/14/2011 
!
! Modules used:
Use kind_mod                  ! Import rknd, iknd specifications
Use gfile_var_pass, Only : g_ip_sign, g_ssimag, g_ssibry, g_ssimag, g_r, g_z, g_dr, g_dz, g_mw, g_mh
Use bspline
Implicit None

! Input/output                      !See above for descriptions
Integer(iknd),Intent(in) :: Npts
Real(rknd),Dimension(Npts),Intent(in)  :: R1,Z1
Real(rknd),Dimension(Npts),Intent(out) :: psiout, psiNout
Integer(iknd),Intent(out) :: ierr

! Local Scalars
Integer(iknd) :: ir,iz,ii
Real(rknd) :: dir,diz

!- End of header -------------------------------------------------------------

If (.not. Allocated(g_r) ) Then
  Write(*,*) 'G VARIABLES NOT ALLOCATED, EXITING FROM get_psi_bicub!'
  Stop
Endif

ierr = 0
Do ii = 1,Npts 

  ! update for vec.
  ir = Floor( (R1(ii)-g_r(1))/g_dr ) + 1
  iz = Floor( (Z1(ii)-g_z(1))/g_dz ) + 1

  ! Check for points off grid
  If ( (ir .le. 2) .or. (ir .ge. g_mw - 1) ) Then
    Write(*,'(3(a,f12.3),a)') 'psi_geq: Point off grid in R: R = ',R1(ii),&
         '. [Rmin,Rmax] = [',g_r(1),',',g_r(g_mw),']'
    ierr = 1
    psiout(ii:Npts) = 0.d0
    psiNout(ii:Npts) = 0.d0
    return
  Endif
  If ( (iz .le. 1) .or. (iz .ge. g_mh - 1) ) Then
    Write(*,'(3(a,f12.3),a)') 'psi_geq: Point off grid in Z: Z = ',Z1(ii),&
         '. [Zmin,Zmax] = [',g_z(1),',',g_z(g_mh),']'
    ierr = 1
    psiout(ii:Npts) = 0.d0
    psiNout(ii:Npts) = 0.d0
    return
  Endif

  dir = (R1(ii) - g_r(ir))/g_dr
  diz = (Z1(ii) - g_z(iz))/g_dz
  psiout(ii) = psi_bi(iz + g_mh*(ir-1),dir,diz)
  psiNout(ii) = (psiout(ii)*g_ip_sign - g_ssimag)/(g_ssibry-g_ssimag)

Enddo

End Subroutine get_psi_bicub

!-----------------------------------------------------------------------------
!+ returns psi and derivatives from gfile bfield at R,Z
!-----------------------------------------------------------------------------
Subroutine get_psi_derivs_bicub(R1,Z1,Npts,psiout,dpsidr,dpsidz,ierr)
! Description: 
!  
! Output:
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/14/2011  Ported from Canik's idl routines.  JDL
! 
! Author(s): J.D. Lore 04/14/2011 
!
! Modules used:
Use kind_mod                  ! Import rknd, iknd specifications
Use gfile_var_pass
Implicit None

! Input/output                      !See above for descriptions
Integer(iknd),Intent(in) :: Npts
Real(rknd),Dimension(Npts),Intent(in)  :: R1,Z1
Real(rknd),Dimension(Npts),Intent(out) :: psiout, dpsidr, dpsidz
Integer(iknd),Intent(out) :: ierr

! Local Scalars
Integer(iknd) :: ir,iz,index,ii
Real(rknd) :: dir,diz
! Local parameters               
Real(rknd), Parameter :: TWO = 2._rknd, THREE = 3._rknd  

!- End of header -------------------------------------------------------------

If (.not. Allocated(g_r) ) Then
  Write(*,*) 'G VARIABLES NOT ALLOCATED, EXITING FROM get_psi_derivs_bicub!'
  Stop
Endif

ierr = 0
Do ii = 1,Npts 

  ! update for vec.
  ir = Floor( (R1(ii)-g_r(1))/g_dr ) + 1
  iz = Floor( (Z1(ii)-g_z(1))/g_dz ) + 1

  ! Check for points off grid
  If ( (ir .le. 2) .or. (ir .ge. g_mw - 1) ) Then
    Write(*,'(3(a,f12.3),a)') 'psi deriv: Point off grid in R: R = ',R1(ii),&
         '. [Rmin,Rmax] = [',g_r(1),',',g_r(g_mw),']'
    ierr = 1
    psiout = 0.d0
    return
  Endif
  If ( (iz .le. 1) .or. (iz .ge. g_mh - 1) ) Then
    Write(*,'(3(a,f12.3),a)') 'psi deriv: Point off grid in Z: Z = ',Z1(ii),&
         '. [Zmin,Zmax] = [',g_z(1),',',g_z(g_mh),']'
    ierr = 1
    psiout = 0.d0
    return
  Endif

  dir = (R1(ii) - g_r(ir))/g_dr
  diz = (Z1(ii) - g_z(iz))/g_dz
  index = iz + g_mh*(ir-1)
  psiout(ii) = psi_bi(index,dir,diz)
  dpsidr(ii) = dsdr_bi(index,dir,diz)
  dpsidz(ii) = dsdz_bi(index,dir,diz)
Enddo

End Subroutine get_psi_derivs_bicub




Function psi_bi(index,dir,diz)
Use kind_mod
Use gfile_var_pass, Only : g_bicub_coeffs
Implicit None
Integer(iknd), Intent(In) :: index
Real(rknd), Intent(In) :: dir, diz
Real(rknd) :: psi_bi
psi_bi = g_bicub_coeffs(index,1,1)   + g_bicub_coeffs(index,2,1)*dir       + &
     g_bicub_coeffs(index,3,1)*dir*dir       + g_bicub_coeffs(index,4,1)*dir*dir*dir      + &
     g_bicub_coeffs(index,1,2)*diz   + g_bicub_coeffs(index,2,2)*dir*diz   + &
     g_bicub_coeffs(index,3,2)*dir*dir*diz   + g_bicub_coeffs(index,4,2)*dir*dir*dir*diz  + &
     g_bicub_coeffs(index,1,3)*diz*diz  + g_bicub_coeffs(index,2,3)*dir*diz*diz  + &
     g_bicub_coeffs(index,3,3)*dir*dir*diz*diz  + g_bicub_coeffs(index,4,3)*dir*dir*dir*diz*diz + &
     g_bicub_coeffs(index,1,4)*diz*diz*diz  + g_bicub_coeffs(index,2,4)*dir*diz*diz*diz  + &
     g_bicub_coeffs(index,3,4)*dir*dir*diz*diz*diz  + g_bicub_coeffs(index,4,4)*dir*dir*dir*diz*diz*diz
End Function psi_bi


Function dsdr_bi(index,dir,diz)
Use kind_mod
Use gfile_var_pass, Only : g_bicub_coeffs, g_dr
Implicit None
Integer(iknd), Intent(In) :: index
Real(rknd), Intent(In) :: dir, diz
Real(rknd) :: dsdr_bi
Real(rknd), Parameter :: TWO = 2._rknd, THREE = 3._rknd  
dsdr_bi = (g_bicub_coeffs(index,2,1)       + TWO*g_bicub_coeffs(index,3,1)*dir +    &
     THREE*g_bicub_coeffs(index,4,1)*dir*dir      + &
     g_bicub_coeffs(index,2,2)*diz   + TWO*g_bicub_coeffs(index,3,2)*dir*diz   + &
     THREE*g_bicub_coeffs(index,4,2)*dir*dir*diz  + &
     g_bicub_coeffs(index,2,3)*diz*diz  + TWO*g_bicub_coeffs(index,3,3)*dir*diz*diz  + &
     THREE*g_bicub_coeffs(index,4,3)*dir*dir*diz*diz + &
     g_bicub_coeffs(index,2,4)*diz*diz*diz  + TWO*g_bicub_coeffs(index,3,4)*dir*diz*diz*diz + &
     THREE*g_bicub_coeffs(index,4,4)*dir*dir*diz*diz*diz)/g_dr
End Function dsdr_bi


Function dsdz_bi(index,dir,diz)
Use kind_mod, Only : iknd, rknd
Use gfile_var_pass, Only : g_bicub_coeffs, g_dz
Implicit None
Integer(iknd), Intent(In) :: index
Real(rknd), Intent(In) :: dir, diz
Real(rknd) :: dsdz_bi
Real(rknd), Parameter :: TWO = 2._rknd, THREE = 3._rknd  
dsdz_bi = (g_bicub_coeffs(index,1,2)                 + g_bicub_coeffs(index,2,2)*dir       + &
                g_bicub_coeffs(index,3,2)*dir*dir            + g_bicub_coeffs(index,4,2)*dir*dir*dir      + &
            TWO*g_bicub_coeffs(index,1,3)*diz         + TWO*g_bicub_coeffs(index,2,3)*dir*diz   + &
            TWO*g_bicub_coeffs(index,3,3)*dir*dir*diz    + TWO*g_bicub_coeffs(index,4,3)*dir*dir*dir*diz  + &
          THREE*g_bicub_coeffs(index,1,4)*diz*diz      + THREE*g_bicub_coeffs(index,2,4)*dir*diz*diz  + &
        THREE*g_bicub_coeffs(index,3,4)*dir*dir*diz*diz + THREE*g_bicub_coeffs(index,4,4)*dir*dir*dir*diz*diz)/g_dz
End Function dsdz_bi


End Module g3d_module



