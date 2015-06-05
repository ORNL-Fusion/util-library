!-----------------------------------------------------------------------------
!
! Tests xpoint finding routine
! -- JDL
!
Program test_read_gfile
  ! 6/4/15
  ! Author(s): J. Lore - current

Use kind_mod                     ! Import rknd, iknd specifications
Use g3d_module
Use util_routines, Only: calc_sep
Use fieldline_follow_mod, Only : bfield_method
Use m3dc1_routines_mod
Implicit None

Character(Len=120) :: gfilename, fname, fname_out
!- End of header -------------------------------------------------------------


!gfilename = '/home/jjl/gfiles/DIII-D/g160884.03014_251'
gfilename = '/home/jjl/gfiles/DIII-D/g160884.05009_537'

Call readg_g3d(gfilename)
fname_out = "sep_data_gfile.out"
Call calc_sep(.false.,fname_out)


m3dc1_itime = 1
m3dc1_factor = 0.d0
!fname = '/home/jjl/m3dc1/160884/3000/n=3_2f_odd/C1.h5'
fname = '/home/jjl/m3dc1/160884/5000/n=3_2f_odd/C1.h5'
m3dc1_field_type = 0
bfield_method = 4
Call prepare_m3dc1_fields(fname)
fname_out = 'sep_data_m3dc1_tot.out'
Call calc_sep(.false.,fname_out)
Call close_m3dc1_fields


Write(*,'(/a)') '-----------------------'
m3dc1_field_type = 1
bfield_method = 3
Call prepare_m3dc1_fields(fname)
fname_out = 'sep_data_m3dc1_super.out'
Call calc_sep(.false.,fname_out)
Call close_m3dc1_fields



End program test_read_gfile

