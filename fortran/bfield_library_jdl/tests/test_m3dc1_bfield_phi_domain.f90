Program test_m3dc1_bfield_phi_domain

Use kind_mod
Use m3dc1_routines_mod
Implicit None

Real(rknd), Allocatable :: bg2(:,:), rg(:), zg(:), pg(:)
Integer(iknd) :: ierr_b, num_pts, i
Character(120) :: fname
real(rknd) :: pi = 3.14159

num_pts = 1
Allocate(bg2(num_pts,3),rg(num_pts),zg(num_pts),pg(num_pts))


fname = '/home/jjl/m3dc1/160884/3000/n=3_2f_odd/C1.h5'

m3dc1_itime = 1
m3dc1_factor = 5.09296d0

Call prepare_m3dc1_fields(fname)


write(*,*) '    R (m)       Phi (deg)   Z (m)       Br (T)      Bphi (T)    Bz (T)'
rg = 2.2d0
zg = 0.d0
pg = 0.d0*pi/180.d0
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2

rg = 2.2d0
zg = 0.d0
pg = 120.d0*pi/180.d0
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2

rg = 2.2d0
zg = 0.d0
pg = 240.d0*pi/180.d0
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2

rg = 2.2d0
zg = 0.d0
pg = 360.d0*pi/180.d0
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2

rg = 2.2d0
zg = 0.d0
pg = -120.d0*pi/180.d0
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2


rg = 2.2d0
zg = 0.d0
pg = 1200.d0*pi/180.d0
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2

rg = 2.2d0
zg = 0.d0
pg = -1200.d0*pi/180.d0
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2

rg = 2.2d0
zg = 0.d0
pg = 22.d0*pi/180.d0
Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2


!write(*,*) 'adsfsadfsda'
!pg = 0.d0
!do i = 1,1000
!  Call bfield_m3dc1(rg,pg,zg,num_pts,bg2,ierr_b)
!  write(*,'(6f12.5)') rg,pg*180./pi,zg,bg2
!  pg = pg + 1.
!enddo


Deallocate(bg2,rg,zg,pg)

End Program test_m3dc1_bfield_phi_domain
