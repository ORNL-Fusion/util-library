Program test_m3dc1_bfield

Use kind_mod
Use m3dc1_routines_mod
Implicit None
Real(rknd), Allocatable :: bg2(:,:), rg(:), zg(:), pg(:)
Integer(iknd) :: ierr_b, num_pts, i
Character(120) :: fname


num_pts = 1
Allocate(bg2(num_pts,3),rg(num_pts),zg(num_pts),pg(num_pts))

rg = 1.101516412700766d0         
zg = -0.242571188190119d0
pg = 0.314159265358979d0

fname = '/home/jjl/m3dc1/mesh21a_kap6_amu6_d0/C1.h5'

m3dc1_itime = 0
m3dc1_factor = 4.d0

Call prepare_m3dc1_fields(fname)

Do i = 1,1000
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)     
Write(*,'(3f21.15)') bg2
Enddo 

Deallocate(bg2,rg,zg,pg)

End Program Test_m3dc1_bfield
