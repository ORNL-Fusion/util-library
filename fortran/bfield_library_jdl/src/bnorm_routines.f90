Module bnorm_routines
  Implicit None
  Private
  Public :: get_pest_coords
Contains
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  Subroutine get_pest_coords(g,pnwant,ntheta,closeit,rpest,zpest,jpest)
    Use kind_mod, Only : int32, real64
    Use phys_const, Only : pi
    Use bfield, Only : bfield_type, g_type
    Use fieldline_follow_mod, Only : follow_fieldlines_rzphi
    Use math_geo_module, Only : linear_interp
    Use g3d_module, Only: calc_RZ_at_psiN_theta1d
    Implicit None
    Type(g_type), Intent(In) :: g
    Real(real64), Intent(In) :: pnwant(:)
    Integer(int32), Intent(In) :: ntheta
    Logical, Intent(In) :: closeit
    Real(real64), Intent(Out) :: rpest(:,:)
    Real(real64), Intent(Out) :: zpest(:,:)
    Real(real64), Intent(Out) :: jpest(:,:)

    Type(bfield_type) :: bfield
    Real(real64), Allocatable :: phistart_arr(:), fl_r(:,:), fl_z(:,:), fl_p(:,:)
    Integer(int32), Allocatable :: ilg(:), fl_ierr(:)

    Real(real64) :: dphi
    Integer(int32) :: npn,i,ierr,nsteps
    Real(real64), Allocatable :: rstart(:), zstart(:), qpsi(:)


    npn = Size(pnwant)
    Allocate(rstart(npn),zstart(npn),qpsi(npn))
    
    Call calc_RZ_at_psiN_theta1d(g,pnwant,0.d0,rstart,zstart)

    Do i = 1,npn
      Call linear_interp(g%pn,g%qpsi,g%mw,pnwant(i),qpsi(i),ierr)
    Enddo
    
    Allocate(ilg(npn),fl_ierr(npn),phistart_arr(npn))
    Allocate(fl_r(npn,nsteps+1),fl_z(npn,nsteps+1),fl_p(npn,nsteps+1))

    phistart_arr = 0.d0
    dphi = 0.5d0*pi/180.d0
    nsteps = Nint(1.1d0*2.d0*pi*Maxval(qpsi)/dphi)
    Write(*,*) 'JDL nsteps',nsteps

    bfield%g = g
    bfield%method = 0
    Call follow_fieldlines_rzphi(bfield,rstart,zstart,phistart_arr,npn,&
         dphi,nsteps,fl_r,fl_z,fl_p,fl_ierr,ilg)

    !f.theta = atan2(f.z - g.zmaxis,f.r - g.rmaxis);
    !f.theta(f.theta<0) = f.theta(f.theta<0)+2*pi;
    
    Deallocate(ilg,fl_ierr,phistart_arr,fl_r,fl_z,fl_p)
    Deallocate(rstart,zstart,qpsi)

    rpest = 0.
    zpest = 0.
    jpest = 0.
    
  End Subroutine get_pest_coords


End Module bnorm_routines
