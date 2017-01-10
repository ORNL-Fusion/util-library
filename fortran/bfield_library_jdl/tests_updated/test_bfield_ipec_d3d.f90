!-----------------------------------------------------------------------------
!
! -- JDL
!
Program test_bfield_ipec_d3d
Use kind_mod, Only: int32, real64
Use ipec_module
Implicit None

Character(Len=100) :: run_path
Integer(int32) :: ierr,ifield_type, ipec_type

Real(real64) :: Bout(1,3)

!- End of header -------------------------------------------------------------

run_path = '/home/jjl/IPEC/164723/low/gpec'
!fname = '/home/jjl/IPEC/164723/low/gpec/ipec_eqbrzphi_n3.out'
!Write(*,*) 'Reading ipec file:', fname
!Call read_ipec_field_file(fname,0)
ipec_type = 1
Call open_ipec_fields(run_path,ipec_type)
ifield_type = 1

Call bfield_ipec((/2.1d0/),(/0.1d0/),(/0.05d0/),1,Bout,ierr,ifield_type)
WRite(*,*) Bout
Call close_ipec_fields



End program test_bfield_ipec_d3d

