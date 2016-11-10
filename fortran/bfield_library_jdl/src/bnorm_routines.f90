Module bnorm_routines
  Implicit None
  Private
  Public :: get_pest_coords
Contains
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  Subroutine get_pest_coords(g,pnwant,ntheta,rpest,zpest,jpest)
    Use kind_mod, Only : int32, real64
    Use phys_const, Only : pi
    Use bfield, Only : bfield_type, g_type
    Use fieldline_follow_mod, Only : follow_fieldlines_rzphi
    Use math_geo_module, Only : linear_interp
    Use g3d_module, Only: calc_RZ_at_psiN_theta1d,bfield_geq_bicub
    Implicit None
    Type(g_type), Intent(In) :: g
    Real(real64), Intent(In) :: pnwant(:)
    Integer(int32), Intent(In) :: ntheta

    Real(real64), Intent(Out) :: rpest(:,:)
    Real(real64), Intent(Out) :: zpest(:,:)
    Real(real64), Intent(Out) :: jpest(:,:)

    Type(bfield_type) :: bfield
    Real(real64), Allocatable :: phistart_arr(:), &
         fl_r(:,:), fl_z(:,:), fl_p(:,:), fl_theta(:,:)
    Integer(int32), Allocatable :: ilg(:), fl_ierr(:)

    Real(real64) :: dphi,phi_th0,phi_pest,bpol(ntheta),Bval(ntheta,3)
    Integer(int32) :: npn,i,ierr,nsteps,icross,j
    Real(real64), Allocatable :: rstart(:), zstart(:), qpsi(:)

    npn = Size(pnwant)
    Allocate(rstart(npn),zstart(npn),qpsi(npn))
    
    Call calc_RZ_at_psiN_theta1d(g,pnwant,0.d0,rstart,zstart)
    Do i = 1,npn
      Call linear_interp(g%pn,g%qpsi,g%mw,pnwant(i),qpsi(i),ierr)
    Enddo   
    Allocate(ilg(npn),fl_ierr(npn),phistart_arr(npn))
    phistart_arr = 0.d0
    dphi = 0.5d0*pi/180.d0
    nsteps = Nint(1.1d0*2.d0*pi*Maxval(qpsi)/dphi)
    
    Allocate(fl_r(npn,nsteps+1),fl_z(npn,nsteps+1),fl_p(npn,nsteps+1))
    Allocate(fl_theta(npn,nsteps+1))    
    
    bfield%g = g
    bfield%method = 0
    Call follow_fieldlines_rzphi(bfield,rstart,zstart,phistart_arr,npn,&
         dphi,nsteps,fl_r,fl_z,fl_p,fl_ierr,ilg)    

    fl_theta = Atan2(fl_z - g%zmaxis,fl_r - g%rmaxis)
    fl_theta = Mod(fl_theta + 2.d0*pi,2.d0*pi)
    
    Do i = 1,npn

      ! Get one period in theta: 0:2*pi
      If (fl_theta(i,2) < fl_theta(i,1)) Then
        Write(*,*) "Need to take reverse theta into account"
        Stop "Quitting in get_pest_coords"
      Endif
      icross = 0
      Do j = 2,nsteps+1
        If (fl_theta(i,j) - fl_theta(i,j-1) .LT. 0.d0) Then
          icross = j
          Exit
        Endif
      Enddo
      If (icross .eq. 0) Then
        Write(*,*) 'Did not find theta crossing, maybe increase nsteps?'
        Stop "Quitting in get_pest_coords"
      Endif

      ! Find phi where theta crosses 2*pi
      fl_theta(i,icross) = fl_theta(i,icross) + 2.d0*pi
      Call linear_interp(fl_theta(i,1:icross),fl_p(i,1:icross),icross,&
           2.d0*pi,phi_th0,ierr)

      Do j = 0,ntheta-1
        phi_pest = phi_th0*Real(j,real64)/Real(ntheta,real64)
        Call linear_interp(fl_p(i,1:icross),fl_r(i,1:icross),icross, &
             phi_pest,rpest(i,j+1),ierr)
        if (ierr .eq. 1) Then
          Write(*,*) 'linear interp err!!!!!!!!!!!!!!!!'
        Endif
        Call linear_interp(fl_p(i,1:icross),fl_z(i,1:icross),icross, &
             phi_pest,zpest(i,j+1),ierr)
        if (ierr .eq. 1) Then
          Write(*,*) 'linear interp err!!!!!!!!!!!!!!!!'
        Endif
      Enddo

      Bval = 0.d0
      Call bfield_geq_bicub(g,rpest(i,:),zpest(i,:),ntheta,Bval,ierr)
      bpol = Sqrt(Bval(:,1)**2 + Bval(:,2)**2) 
      jpest(i,:) = Abs(qpsi(i)/(g%rzero*g%bcentr))*bpol(:)*rpest(i,:)**3

    Enddo
    
    Deallocate(ilg,fl_ierr)
    Deallocate(phistart_arr,fl_r)
    Deallocate(fl_z,fl_p,fl_theta)
    Deallocate(rstart,zstart,qpsi)
  End Subroutine get_pest_coords


End Module bnorm_routines
