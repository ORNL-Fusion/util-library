!-----------------------------------------------------------------------------
!
!   Routines/modules related to M3DC1
!
!
!   Contains:
!     Subroutine prepare_m3dc1_fields
!     Subroutine bfield_m3dc1
!     Subroutine bfield_m3dc1_2d
!     Subroutine bfield_m3dc1_pert
!     Subroutine close_m3dc1_fields
!     Subroutine calc_psi_m3dc1_2d
!     Subroutine calc_psiN_m3dc1_2d
!
!    m3dc1_itime  : time slice to read (0 = vacuum only, 1 = with response)
!    m3dc1_factors : factor by which to multiply perturbed (linear) part of solution
!    m3dc1_toroidal_on_err : On error flag from m3dc1 routine, return B=Bt=1
!    m3dc1_field_type : return total field (0) or perturbed part only (1)
!-----------------------------------------------------------------------------
Module M3DC1_routines_mod
Use kind_mod, Only: real64, int32
Use fusion_io
Implicit None
Integer(int32), Parameter  :: num_sets_max = 10
Integer(int32), Save :: m3dc1_itime = -1
Integer(int32), Save :: m3dc1_field_type = -1
Real(real64), Save   :: m3dc1_factors(num_sets_max) = 0.d0
Real(real64), Save   :: m3dc1_phases_deg(num_sets_max)  = 0.d0
Logical :: m3dc1_toroidal_on_err = .false.

Integer, Private, Save :: num_sets
Integer(int32), Private, Save :: isrc(num_sets_max), imag(num_sets_max), imag_2d, ivec_2d
Real(real64), Private, Save   :: psi_axis, psi_lcfs
type(fio_search_hint), Private, Save :: hint(num_sets_max)

Contains

!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine prepare_m3dc1_fields(filenames)
!
! Description:
! Input:
!
! Output:
!
! Calls:
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!   
! Author(s): J.D. Lore 1/27/2014
Use kind_mod, Only: int32
Use fusion_io
Implicit None
Character(len=*), Intent(In) :: filenames(:)
Integer(int32) :: ierr, ipsi_axis, ipsi_lcfs, iset
!- End of header -------------------------------------------------------------

! Open M3D-C1 source
! ~~~~~~~~~~~~~~~~~~
num_sets = size(filenames)
Write(*,*) 'In M3DC1 fields, reading ',num_sets,' files'
Do iset = 1,num_sets
  Write(*,'(a,a)') ' Reading ', Trim(filenames(iset))
  Call fio_open_source_f(FIO_M3DC1_SOURCE, Trim(filenames(iset)), isrc(iset), ierr)
  If (ierr .ne. 0) Then
    Write(*,*) 'Error opening m3dc1 source file, exiting from prepare_m3dc1_fields'
    Stop
  Endif
  
  ! Get and set options
  ! ~~~~~~~~~~~
  Call fio_get_options_f(isrc(iset), ierr)
  If ( m3dc1_itime == -1 ) Then
    Write(*,*) 'm3dc1_itime has not been set!',m3dc1_itime
    Stop
  Endif
  Call fio_set_int_option_f(FIO_TIMESLICE, m3dc1_itime, ierr)

  if (iset == 1) Then
    ! read 2D equilibrium field data
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Call fio_set_int_option_f(FIO_PART, FIO_EQUILIBRIUM_ONLY, ierr)
    Call fio_get_field_f(isrc(iset), FIO_MAGNETIC_FIELD,   imag_2d, ierr)  
    Call fio_get_field_f(isrc(iset), FIO_VECTOR_POTENTIAL, ivec_2d, ierr)
    Call fio_get_series_f(isrc(iset), FIO_MAGAXIS_PSI, ipsi_axis, ierr)
    Call fio_get_series_f(isrc(iset), FIO_LCFS_PSI, ipsi_lcfs, ierr)

    Call fio_eval_series_f(ipsi_axis, 0._real64, psi_axis, ierr)
    Call fio_eval_series_f(ipsi_lcfs, 0._real64, psi_lcfs, ierr)
!    Write(*,*) 'Psi at magnetic axis: ', psi_axis
!    Write(*,*) 'Psi at lcfs: ', psi_lcfs

    Call fio_close_series_f(ipsi_axis, ierr)
    Call fio_close_series_f(ipsi_lcfs, ierr)
    
  Endif


  ! read 3D field data
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  If (m3dc1_field_type .eq. 0) Then
    Write(*,*) 'M3DC1 returning total field'
    Call fio_set_int_option_f(FIO_PART, FIO_TOTAL, ierr) 
  Elseif (m3dc1_field_type .eq. 1) Then
    Write(*,*) 'M3DC1 returning perturbed field only'
    Call fio_set_int_option_f(FIO_PART, FIO_PERTURBED_ONLY, ierr)
  Else
    Write(*,*) 'Bad value for m3dc1_field_type:',m3dc1_field_type
    Stop
  Endif
  Write(*,*) 'M3DC1 Amplitude scale:  ', m3dc1_factors(iset)
  Write(*,*) 'M3DC1 phase shift (deg):', m3dc1_phases_deg(iset)
  Call fio_set_real_option_f(FIO_LINEAR_SCALE, m3dc1_factors(iset), ierr)
  Call fio_get_field_f(isrc(iset), FIO_MAGNETIC_FIELD, imag(iset), ierr)


!  ! Set up pert only field
!  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Call fio_set_int_option_f(FIO_PART, FIO_PERTURBED_ONLY, ierr)
!  Call fio_set_real_option_f(FIO_LINEAR_SCALE, m3dc1_factors(iset), ierr)
!  Call fio_get_field_f(isrc(iset), FIO_MAGNETIC_FIELD, imag_pert(iset), ierr)

  ! Initialize hint
  call fio_allocate_search_hint_f(isrc(iset), hint(iset), ierr)

  
Enddo



End Subroutine prepare_m3dc1_fields


!-----------------------------------------------------------------------------
!+ Evaluate B(r,phi,z) using m3dc1 perturbed (or total) fields, phi in radians.
!-----------------------------------------------------------------------------
Subroutine bfield_m3dc1(r,phi,z,Npts,Bout,ierr)
! Output:
!   Bout = (:,[Br,Bz,Bt])
Use fusion_io
Use kind_mod, Only: int32, real64
Use phys_const, Only: pi
Implicit None
Real(Real64), Intent(In), Dimension(Npts) :: r, z, phi
Integer(int32), Intent(In) :: Npts
Real(real64), Intent(Out), Dimension(Npts,3) :: Bout
Integer(int32), Intent(Out) :: ierr
! Local variables
Real(real64) :: x(3), b_tmp(3)
Integer(int32) :: ierr_b, i, iset
ierr = 0
Bout = 0._real64
Do i=1,Npts
  x(1) = r(i)
  x(2) = phi(i)
  x(3) = z(i)  
  
  Do iset = 1,num_sets
    x(2) = x(2) + m3dc1_phases_deg(iset)*pi/180._real64
    b_tmp = 0._Real64
    Call fio_eval_field_f(imag(iset), x, b_tmp, ierr_b, hint=hint(iset))  ! b_tmp(R,phi,Z)
    
    If (ierr_b .ne. 0) Then
      If (m3dc1_toroidal_on_err) Then
        b_tmp = 0._real64
        b_tmp(2) = 1._real64 
      Else
        ierr = 1
      Endif
    Endif
    Bout(i,1) = Bout(i,1) + b_tmp(1) !Br
    Bout(i,2) = Bout(i,2) + b_tmp(3) !Bz
    Bout(i,3) = Bout(i,3) + b_tmp(2) !Bphi
  Enddo
Enddo

End Subroutine bfield_m3dc1

!-----------------------------------------------------------------------------
!+ Evaluate B(r,phi,z) using m3dc1 perturbed ONLY fields, phi in radians.
!-----------------------------------------------------------------------------
!Subroutine bfield_m3dc1_pert(r,phi,z,Npts,Bout,ierr)
! Will need to define imag_pert above
!End Subroutine bfield_m3dc1_pert


!-----------------------------------------------------------------------------
!+ Evaluate B(r,z) using Equilibrium only M3DC1 fields
!-----------------------------------------------------------------------------
Subroutine bfield_m3dc1_2d(r,z,Npts,Bout,ierr)
! Output:
!   Bout = (:,[Br,Bz,Bt])
Use fusion_io
Use kind_mod, Only: int32, real64
Implicit None
Real(Real64), Intent(In), Dimension(Npts) :: r, z
Integer(int32), Intent(In) :: Npts
Real(real64), Intent(Out), Dimension(Npts,3) :: Bout
Integer(int32), Intent(Out) :: ierr
! Local variables
Real(real64) :: x(3), b_tmp(3)
Integer(int32) :: ierr_b, i
ierr = 0
Bout = 0._real64
Do i=1,Npts
  x(1) = r(i)
  x(2) = 0._real64
  x(3) = z(i)
  b_tmp = 0._real64

  Call fio_eval_field_f(imag_2d, x, b_tmp, ierr_b, hint=hint(1))  ! b_tmp(R,phi,Z)

  If (ierr_b .ne. 0) Then
    If (m3dc1_toroidal_on_err) Then
      b_tmp = 0._real64
      b_tmp(2) = 1._real64 
    Else
      ierr = 1
    Endif
  Endif
  Bout(i,1) = b_tmp(1)
  Bout(i,2) = b_tmp(3)
  Bout(i,3) = b_tmp(2)  
Enddo

End Subroutine bfield_m3dc1_2d



!-----------------------------------------------------------------------------
!+ Evaluate psiN(r,z) using Equilibrium only M3DC1 fields
!-----------------------------------------------------------------------------
Subroutine calc_psiN_m3dc1_2d(r,z,Npts,PsiN,ierr)
Use fusion_io
Use kind_mod, Only: int32, real64
Implicit None
Real(Real64), Intent(In), Dimension(Npts) :: r, z
Integer(int32), Intent(In) :: Npts
Real(real64), Intent(Out), Dimension(Npts) :: PsiN
Integer(int32), Intent(Out) :: ierr
! Local variables
Real(real64) :: psi(Npts)
Integer(int32) :: ierr_psi
ierr = 0
Call calc_psi_m3dc1_2d(r,z,Npts,psi,ierr_psi)
If (ierr_psi .ne. 0) ierr = 1
psiN = (psi - psi_axis)/(psi_lcfs - psi_axis)
End Subroutine calc_psiN_m3dc1_2d

!-----------------------------------------------------------------------------
!+ Evaluate psi(r,z) using Equilibrium only M3DC1 fields
!-----------------------------------------------------------------------------
Subroutine calc_psi_m3dc1_2d(r,z,Npts,Psi,ierr)
Use fusion_io
Use kind_mod, Only: int32, real64
Implicit None
Real(Real64), Intent(In), Dimension(Npts) :: r, z
Integer(int32), Intent(In) :: Npts
Real(real64), Intent(Out), Dimension(Npts) :: Psi
Integer(int32), Intent(Out) :: ierr
! Local variables
Real(real64) :: x(3), a_tmp(3)
Integer(int32) :: ierr_a, i
ierr = 0
Psi = 0._real64
Do i=1,Npts
  x(1) = r(i)
  x(2) = 0._real64
  x(3) = z(i)
  a_tmp = 0._real64
  
  Call fio_eval_field_f(ivec_2d, x, a_tmp, ierr_a, hint=hint(1))  ! A(R,phi,Z)
  If (ierr_a .ne. 0) Then
    ierr = 1
  Endif
  psi(i) = a_tmp(2)*r(i)
Enddo

End Subroutine calc_psi_m3dc1_2d


!-----------------------------------------------------------------------------
!+ Close files and deallocate variables associated with M3DC1
!-----------------------------------------------------------------------------
Subroutine close_m3dc1_fields
!
! Description:
! Input:
!
! Output:
!
! Calls:
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!   
! Author(s): J.D. Lore 1/27/2014
Use kind_mod, Only: int32
Use fusion_io
Implicit None
Integer(int32) :: ierr, iset
!- End of header -------------------------------------------------------------

! Close fields and source
! !~~~~~~~~~~~~~~~~~~~~~~
!Call fio_close_field_f(ine, ierr)
!Call fio_close_field_f(ini, ierr)
!Call fio_close_field_f(ipres, ierr)
!Call fio_close_field_f(ipe, ierr)
Do iset = 1,num_sets
  Call fio_close_field_f(imag(iset), ierr)
  If (iset == 1) Then
    Call fio_close_field_f(imag_2d, ierr)
    Call fio_close_field_f(ivec_2d, ierr)
  Endif
  Call fio_close_source_f(isrc(iset), ierr)
  call fio_deallocate_search_hint_f(isrc(iset), hint(iset), ierr)
Enddo


End Subroutine close_m3dc1_fields

End Module M3DC1_routines_mod
