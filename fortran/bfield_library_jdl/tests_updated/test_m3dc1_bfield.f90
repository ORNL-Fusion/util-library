Program test_m3dc1_bfield
Use kind_mod
Use m3dc1_routines_mod
Implicit None
Real(real64), Allocatable :: bg2(:,:), rg(:), zg(:), pg(:)
Integer(int32) :: ierr_b, num_pts, i, m3dc1_nsets
Character(250) :: m3dc1_filenames(2)



m3dc1_nsets = 2  ! MUST EQUAL LENGTH OF FILENAMES AND SCALE_FACTORS!!
m3dc1_filenames(1) = '/home/jjl/m3dc1/164723/03059/n=3/eb1_2f_iu/C1.h5'
m3dc1_filenames(2) = '/home/jjl/m3dc1/164723/03059/n=3/eb1_2f_il/C1.h5'

m3dc1_itime = 1
m3dc1_factors(1) = 3.69
m3dc1_factors(2) = 3.69
m3dc1_toroidal_on_err = .true.
m3dc1_field_type = 0
Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)




num_pts = 1
Allocate(bg2(num_pts,3),rg(num_pts),zg(num_pts),pg(num_pts))
rg = 2.2
zg = 0.
pg = 0.

Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)     

write(*,*) bg2
Deallocate(bg2,rg,zg,pg)

End Program Test_m3dc1_bfield
