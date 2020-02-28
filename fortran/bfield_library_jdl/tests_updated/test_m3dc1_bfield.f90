program test_m3dc1_bfield
Use kind_mod
Use m3dc1_routines_mod
Implicit None
Real(real64), Allocatable :: bg2(:,:), rg(:), zg(:), pg(:)
Integer(int32) :: ierr_b, num_pts, i, m3dc1_nsets
Character(250) :: m3dc1_filenames(4)
Real(Kind=4)  :: tarray(2),tres,tres0

m3dc1_nsets = 4
m3dc1_filenames(1) = '/home/jjl/M3DC1_runs/orlov_2300/n=1/eb1_1f_probeg/C1.h5'
m3dc1_filenames(2) = '/home/jjl/M3DC1_runs/orlov_2300/n=2/eb1_1f_probeg/C1.h5'
m3dc1_filenames(3) = '/home/jjl/M3DC1_runs/orlov_2300/n=3/eb1_1f_probeg/C1.h5'
m3dc1_filenames(4) = '/home/jjl/M3DC1_runs/orlov_2300/n=4/eb1_1f_probeg/C1.h5'
m3dc1_factors(1:4) = 2.
m3dc1_phases_deg(1:4) = 0.
m3dc1_itime = 1  ! 0 vacuum, 1 response
m3dc1_field_type = 1 ! total (0) or perturbed only (1)
m3dc1_toroidal_on_err = .false.

Call Etime(tarray,tres0)
Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)

Call Etime(tarray,tres)
Write(*,*) ' Initialization took ',tres-tres0,' seconds'
Call Etime(tarray,tres0)



num_pts = 1
Allocate(bg2(num_pts,3),rg(num_pts),zg(num_pts),pg(num_pts))
rg = 2.2d0
zg = 0.d0
pg = 30.d0*3.14159d0/180.d0


Write(*,*) 'Evaluating B at [R,Z,phi] = ',rg,zg,pg
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)     
write(*,*) 'Bout [Br,Bz,Bphi] = ',bg2
Call Etime(tarray,tres)
Write(*,*) ' First call took ',tres-tres0,' seconds'
Call Etime(tarray,tres0)


rg = 2.4d0
zg = 0.d0
pg = 30.d0*3.14159d0/180.d0


Write(*,*) 'Evaluating B at [R,Z,phi] = ',rg,zg,pg
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)     
write(*,*) 'Bout [Br,Bz,Bphi] = ',bg2
Call Etime(tarray,tres)
Write(*,*) ' Second call took ',tres-tres0,' seconds'
Call Etime(tarray,tres0)


Deallocate(bg2,rg,zg,pg)

End Program Test_m3dc1_bfield
