!-----------------------------------------------------------------------------
!
! Tests read of afile, spline interpolation to B
! -- JDL
!
Program test_Pavel_beval

! Description: 
!
! History:
! Version   Date        Comment
! -------   ----        -------
! 1.0     05/24/2013   Original code. JL 
! 
! Author(s): J. Lore - current

Use kind_mod                     ! Import rknd, iknd specifications
Use screening_module
Implicit None

Real  :: tarray(2),tres
Real  :: tarray0(2),tres0
Real(rknd) :: Bout(3)
Integer(iknd) :: ierr

Character(Len=200) :: Afile_path


!- End of header -------------------------------------------------------------


Afile_path = '/home/jjl/Pavel/screened-8-19/'
Call setup_screening_vars(Afile_path)
!call etime(tarray0,tres0)
Call read_Afiles(Afile_path)
Call prepare_Afile_splinefits !(Afile_path)
call etime(tarray,tres)
!print *,'Time: ',tres - tres0

!call etime(tarray0,tres0)
!Call read_spline_data(Afile_path)
!call etime(tarray,tres)
!print *,'Time: ',tres - tres0


Call bfield_bspline((/1.2d0/),(/0.1d0/),(/-0.2d0/),1,Bout,ierr)
Write(*,*) Bout,'<--------'

End program test_Pavel_beval

