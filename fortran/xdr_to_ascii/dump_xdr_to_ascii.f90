!
!  routine to read the bfield from Geiger/w7's xdr output file
!  and dump the data to an ascii file.
!
! JDL 7/13/2011 (first part modified from Geiger's(?) readbgrid)
!

program dump_xdr_to_ascii
use kind_mod
Use xdr_routines_mod, Only : readbgrid_xdr, write_xdr_to_ascii
Use setup_bfield_module !, Only : setup_bfield_xdr, xdr_check, xdr_fname, xdr_verbose
Use bfield
implicit none

Real(real64), dimension(1) :: r,z,phi,br,bz,bphi

integer(int32) :: ierr_out
character(len=300) :: filename_in
character(len=300) :: filename_out
!!!logical :: xdr_check 


Namelist / settings_nml / xdr_check, filename_in, filename_out

Open(99,file="settings.nml",status="old")
Read(99,nml=settings_nml)
Close(99)

!!Call readbgrid_xdr(filename_in,xdr_check,.true.)


xdr_fname = filename_in
xdr_verbose = .true.
Call setup_bfield_xdr



!r = 4.3d0
!z = 0.1d0
!phi = 0.2d0
!Call calc_B_rzphi_general(bfield,r,z,phi,1,br,bz,bphi,ierr_out)
!writE(*,*) 'r,z,phi = ',r,z,phi
!writE(*,*) 'Br,Bz,Bphi = ',br,bz,bphi
!stop

Call write_xdr_to_ascii(filename_out)

write(*,*) '>>>>>>done writing file<<<<<'





end program dump_xdr_to_ascii
