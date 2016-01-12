!-----------------------------------------------------------------------------
!
!   Routines/modules related to M3DC1
!
!
!   Contains:
!     Subroutine prepare_m3dc1_fields
!     Subroutine bfield_m3dc1
!     Subroutine close_m3dc1_fields
!
!    m3dc1_itime  : time slice to read (0 = vacuum only, 1 = with response)
!    m3dc1_factor : factor by which to multiply perturbed (linear) part of solution
!    m3dc1_toroidal_on_err : On error flag from m3dc1 routine, return B=Bt=1
!    m3dc1_field_type : return total field (0) or perturbed part only (1)
!-----------------------------------------------------------------------------
Module M3DC1_routines_mod
Use kind_mod, Only: real64, int32
Implicit None
Integer(int32), Save :: m3dc1_itime = -1
Integer(int32), Save :: m3dc1_field_type = -1
Real(real64), Save   :: m3dc1_factor = 0.d0
Logical :: m3dc1_toroidal_on_err = .false.

Integer(int32), Private, Save :: isrc, imag, imag_2d, ivec_2d
Real(real64), Private, Save   :: psi_axis, psi_lcfs

Contains

!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine prepare_m3dc1_fields(filename)
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
Character(len=120), Intent(In) :: filename
Integer(int32) :: ierr, ipsi_axis, ipsi_lcfs
!- End of header -------------------------------------------------------------

! Open M3D-C1 source
! ~~~~~~~~~~~~~~~~~~
Write(*,'(a,a)') 'Reading ', Trim(filename)
Call fio_open_source_f(FIO_M3DC1_SOURCE, Trim(filename), isrc, ierr)
If (ierr .ne. 0) Then
  Write(*,*) 'Error opening m3dc1 source file, exiting from prepare_m3dc1_fields'
  Stop
Endif

! Set options
! ~~~~~~~~~~~
Call fio_get_options_f(isrc, ierr)
If ( m3dc1_itime == -1 ) Then
  Write(*,*) 'm3dc1_itime has not been set!',m3dc1_itime
  Stop
Endif
Call fio_set_int_option_f(FIO_TIMESLICE, m3dc1_itime, ierr)
If (m3dc1_field_type .eq. 0) Then
  Write(*,*) 'm3dc1 returning total field'
  Call fio_set_int_option_f(FIO_PART, FIO_TOTAL, ierr) 
Elseif (m3dc1_field_type .eq. 1) Then
  Write(*,*) 'm3dc1 returning perturbed field'
  Call fio_set_int_option_f(FIO_PART, FIO_PERTURBED_ONLY, ierr)
Else
  Write(*,*) 'Bad value for m3dc1_field_type:',m3dc1_field_type
  Stop
Endif

Write(*,*) 'M3DC1 factor:', m3dc1_factor
Call fio_set_real_option_f(FIO_LINEAR_SCALE, m3dc1_factor, ierr)

! read fields
! ~~~~~~~~~~~
! magnetic field
Call fio_get_field_f(isrc, FIO_MAGNETIC_FIELD, imag, ierr)


! read 2D equilibrium field data
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Call fio_set_int_option_f(FIO_PART, FIO_EQUILIBRIUM_ONLY, ierr)
Call fio_get_field_f(isrc, FIO_MAGNETIC_FIELD,   imag_2d, ierr)  
Call fio_get_field_f(isrc, FIO_VECTOR_POTENTIAL, ivec_2d, ierr)
Call fio_get_series_f(isrc, FIO_MAGAXIS_PSI, ipsi_axis, ierr)
Call fio_get_series_f(isrc, FIO_LCFS_PSI, ipsi_lcfs, ierr)

Call fio_eval_series_f(ipsi_axis, 0._real64, psi_axis, ierr)
Call fio_eval_series_f(ipsi_lcfs, 0._real64, psi_lcfs, ierr)
Write(*,*) 'Psi at magnetic axis: ', psi_axis
Write(*,*) 'Psi at lcfs: ', psi_lcfs

Call fio_close_series_f(ipsi_axis, ierr)
Call fio_close_series_f(ipsi_lcfs, ierr)

End Subroutine prepare_m3dc1_fields


!-----------------------------------------------------------------------------
!+ Evaluate B(r,phi,z) using m3dc1 perturbed fields, phi in radians.
!-----------------------------------------------------------------------------
Subroutine bfield_m3dc1(r,phi,z,Npts,Bout,ierr)
! Output:
!   Bout = (:,[Br,Bz,Bt])
Use fusion_io
Use kind_mod, Only: int32, real64
Implicit None
Real(Real64), Intent(In), Dimension(Npts) :: r, z, phi
Integer(int32), Intent(In) :: Npts
Real(real64), Intent(Out), Dimension(Npts,3) :: Bout
Integer(int32), Intent(Out) :: ierr
! Local variables
Real(real64) :: x(3), b_tmp(3)
Integer(int32) :: ierr_b, i
ierr = 0
Do i=1,Npts
  x(1) = r(i)
  x(2) = phi(i)
  x(3) = z(i) 
  
  Call fio_eval_field_f(imag, x, b_tmp, ierr_b)  ! b_tmp(R,phi,Z)

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

End Subroutine bfield_m3dc1

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
Do i=1,Npts
  x(1) = r(i)
  x(2) = 0._real64
  x(3) = z(i) 
  
  Call fio_eval_field_f(imag_2d, x, b_tmp, ierr_b)  ! b_tmp(R,phi,Z)

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
!+ Evaluate B(r,z) using Equilibrium only M3DC1 fields
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
Real(real64) :: x(3), a_tmp(3)
Integer(int32) :: ierr_a, i
ierr = 0
Do i=1,Npts
  x(1) = r(i)
  x(2) = 0._real64
  x(3) = z(i) 
  
  Call fio_eval_field_f(ivec_2d, x, a_tmp, ierr_a)  ! A(R,phi,Z)
  If (ierr_a .ne. 0) Then
    ierr = 1
  Endif
  psiN(i) = (a_tmp(2)*r(i) - psi_axis)/(psi_lcfs - psi_axis)
Enddo

End Subroutine calc_psiN_m3dc1_2d


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
Integer(int32) :: ierr
!- End of header -------------------------------------------------------------

! Close fields and source
! !~~~~~~~~~~~~~~~~~~~~~~
!Call fio_close_field_f(ine, ierr)
!Call fio_close_field_f(ini, ierr)
!Call fio_close_field_f(ipres, ierr)
!Call fio_close_field_f(ipe, ierr)
Call fio_close_field_f(imag, ierr)
Call fio_close_field_f(imag_2d, ierr)
Call fio_close_field_f(ivec_2d, ierr)
Call fio_close_source_f(isrc, ierr)
End Subroutine close_m3dc1_fields

End Module M3DC1_routines_mod
