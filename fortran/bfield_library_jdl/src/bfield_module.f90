Module bfield_typedef
  Use kind_mod, Only : int32, real64
  Use g_typedef, Only : g_type
  Use coil_typedef, Only : coil_type
  Implicit None
  Private
  Type, Public :: bfield_type
    Integer(int32) :: method = -1
    Type(g_type) :: g
    Type(coil_type) :: coil
    Integer(int32) :: method_2d = -1    ! Method corresponding to AS fields
    Integer(int32) :: method_pert = -1  ! Method corresponding to pert only
    Integer(int32) :: method_save = -1  ! Used to save standard method
    Logical :: method_switched = .false.
  End Type bfield_type
End Module bfield_typedef

!   Set bfield_method to control the fieldline deriviative calls.
!   Appropriate loading must be done before calls to any fieldline following routine 
!   (e.g., gfile loading, rmp coil generation)
!
!    bfield%method == 
!                     0 -- gfile field only
!                     1 -- gfile + rmp coil field
!                     2 -- gfile + Pavel's screened fields with bspline interpolation (not fully implemented!)
!                     3 -- gfile + M3DC1 perturbed field
!                     4 -- M3DC1 only
!                     5 -- M3DC1 2D field only
!                     6 -- Just coils
!                     7 -- ipec eq only
!                     8 -- ipec eq+vacuum
!                     9 -- ipec eq+response
!                    10 -- xpand perturbed
!                    11 -- xpand vacuum
!                    12 -- ipec vacuum -- pert only
!                    13 -- ipec response -- pert only
!                    14 -- VMEC coils file with extcur
!                    15 -- xdr file
!                    15 -- bgrid file (FLARE format)
! 


Module bfield
  Use bfield_typedef, Only : bfield_type
  Use coil_typedef, Only : coil_type
  Use g3d_module, Only : g_type
  Implicit None
  Private
  Public :: bfield_type, coil_type, g_type
  Public :: set_bfield_2d, reset_bfield, set_bfield_pert_only
  Public :: calc_B_rzphi_general

Contains

  Subroutine set_bfield_2d(bfield)
    Implicit None
    Type(bfield_type), Intent(InOut) :: bfield
    Call reset_bfield(bfield)
    bfield%method_save = bfield%method
    bfield%method      = bfield%method_2d
    bfield%method_switched = .true.
    Write(*,*) 'Setting bfield method to 2d:',bfield%method
  End Subroutine set_bfield_2d

  Subroutine set_bfield_pert_only(bfield)
    Implicit None
    Type(bfield_type), Intent(InOut) :: bfield
    Call reset_bfield(bfield)
    bfield%method_save = bfield%method
    bfield%method      = bfield%method_pert
    bfield%method_switched = .true.
    Write(*,*) 'Setting bfield method to pert only:',bfield%method
  End Subroutine set_bfield_pert_only

  Subroutine reset_bfield(bfield)
    Implicit None
    Type(bfield_type), Intent(InOut) :: bfield
    If (bfield%method_switched) Then
      bfield%method = bfield%method_save
      bfield%method_switched = .false.
      Write(*,*) 'Resetting bfield method to:',bfield%method
    Endif
  End Subroutine reset_bfield

  Subroutine calc_B_rzphi_general(bfield,r,z,phi,n,br,bz,bphi,ierr_out)
    Use kind_mod, Only: real64, int32
    Use g3d_module, Only : bfield_geq_bicub
#ifdef HAVE_M3DC1
    Use M3DC1_routines_mod, Only : bfield_m3dc1, bfield_m3dc1_2d
#endif
    Use ipec_module, Only : bfield_ipec
    Use xpand_module, Only: bfield_xpand
    Use biotsavart_module, Only : bfield_bs_cyl
    Use VMEC_routines_mod, Only : bfield_vmec_coils
#ifdef HAVE_FXDR    
    Use xdr_routines_mod, Only : bint_xdr_n
#endif    
    Use bgrid_module, Only : bfield_bgrid
    Implicit None
    Type(bfield_type), Intent(In) :: bfield
    Integer(int32), Intent(In) :: n
    Real(real64), Intent(In) :: r(n),z(n),phi(n)
    Real(real64), Intent(Out) :: br(n),bz(n),bphi(n)
    Integer(int32), Intent(Out), Optional :: ierr_out
    Real(real64) :: btmp(n,3)
    Integer(int32) :: ierr

    ierr = 0
    If (Present(ierr_out)) ierr_out = 0

    br   = 0._real64
    bz   = 0._real64
    bphi = 0._real64
    btmp = 0._real64
    
    Select Case (bfield%method)
    Case (0)
      Call bfield_geq_bicub(bfield%g,r,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
    Case (1)
      Call bfield_bs_cyl(r,phi,z,n,bfield%coil,br,bphi,bz)
      Call bfield_geq_bicub(bfield%g,r,z,n,btmp,ierr)
      br   = btmp(:,1) + br
      bz   = btmp(:,2) + bz 
      bphi = btmp(:,3) + bphi
    Case (2)
      Write(*,*) 'method == 2, pavel screening not implemented'
      Stop "Quitting from bfield general"
    Case (3) ! g+m3dc1
#ifdef HAVE_M3DC1  
      Call bfield_geq_bicub(bfield%g,r,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
      Call bfield_m3dc1(r,phi,z,n,btmp,ierr)
      br   = btmp(:,1) + br
      bz   = btmp(:,2) + bz 
      bphi = btmp(:,3) + bphi
#else
      Stop "Compiled without m3dc1 support"
#endif
   Case (4)  ! m3dc1 only
#ifdef HAVE_M3DC1      
      Call bfield_m3dc1(r,phi,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
#else
      Stop "Compiled without m3dc1 support"
#endif      
    Case (5)  ! m3dc1 2d
#ifdef HAVE_M3DC1  
      Call bfield_m3dc1_2d(r,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
#else
      Stop "Compiled without m3dc1 support"
#endif      
    Case (6) ! just coils
      Call bfield_bs_cyl(r,phi,z,n,bfield%coil,br,bphi,bz)
    Case (7) ! ipec eq only
      Call bfield_ipec(r,phi,z,n,btmp,ierr,0)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
    Case (8) ! ipec vac
      Call bfield_ipec(r,phi,z,n,btmp,ierr,1)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)      
    Case (9) ! ipec pert
      Call bfield_ipec(r,phi,z,n,btmp,ierr,2)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)      
    Case (10) ! xpand eq+response
      Call bfield_xpand(r,phi,z,n,btmp,ierr,0)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)      
    Case (11) ! xpand eq+vac
      Call bfield_xpand(r,phi,z,n,btmp,ierr,1)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
    Case (12) ! ipec vac (pert only)
      Call bfield_ipec(r,phi,z,n,btmp,ierr,3)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)      
    Case (13) ! ipec response (pert only)
      Call bfield_ipec(r,phi,z,n,btmp,ierr,4)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
    Case (14) ! VMEC coils
      Call bfield_vmec_coils(r,phi,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
    Case (15) ! Xdr
#ifdef HAVE_FXDR         
      Call bint_xdr_n(r,phi,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,3)  ! Note order!
      bphi = btmp(:,2)
#else
      Stop "Compiled without fxdr support"
#endif      
    Case (16) ! bgrid
      Call bfield_bgrid(r,phi,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
    Case Default
      Write(*,*) 'Unknown bfield%method:',bfield%method
      Stop "Exiting from bfield general"
    End Select
    If (Present(ierr_out)) ierr_out = ierr
  End Subroutine calc_B_rzphi_general
End Module bfield

