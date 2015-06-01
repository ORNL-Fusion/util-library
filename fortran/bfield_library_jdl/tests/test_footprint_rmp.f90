!-----------------------------------------------------------------------------
!
! Tests footprints
! -- JDL
!
Program test_footprint_rmp

! Description: 
!
! History:
! Version   Date        Comment
! -------   ----        -------
! 1.0     05/24/2013   Original code. JL 
! 
! Author(s): J. Lore - current

Use kind_mod                     ! Import rknd, iknd specifications
Use gfile_var_pass
Use g3d_module
Use rmp_module
Use fieldline_follow_mod
Use NSTX_routines_mod
Use math_geo_module

Implicit None

! Local variables (scalar)
Real(rknd),Allocatable :: Rstarts(:),Zstarts(:),phistarts(:),Lc(:,:)
Real(rknd) :: Rstart(1),Zstart(1),phistart(1)
Integer(iknd) :: nphi, Nstarts, i, j, nsteps,ierr_b(1), ilg(1), inside,ithit,hit_jpt,jpt

Character(Len=100) :: gfilename

Integer :: ntest,ntorpts
Real(rknd) :: dphi, taper(6)
Real(rknd),Allocatable,Dimension(:,:) :: rout,zout,phiout

!- End of header -------------------------------------------------------------


!gfilename = '../bfield_files/g135183.00433'
gfilename = '../bfield_files/g140508.00403'
Call readg_g3d(gfilename)

Write(*,*) 'building coils'
taper = (/1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0/)
taper = taper*1000.d0
ntorpts = 5
Allocate(rmp_coil(6*(2*ntorpts+1),3))
Allocate(rmp_current(6*(2*ntorpts+1)))
Call build_nstx_rwmcoils_jl(taper,ntorpts,rmp_coil,rmp_current,rmp_ncoil_pts)

! footprint settings
Nstarts = 400
nphi = 500
Allocate(Rstarts(Nstarts),Zstarts(Nstarts))
Allocate(phistarts(nphi))
Allocate(Lc(Nstarts,nphi))
Rstarts = rlinspace(0.35d0,0.5d0,Nstarts)
Zstarts = -1.603d0
phistarts = rlinspace(0.d0,2.0944d0,nphi)

! fieldline following settings
rmp_method = 1
dphi = 0.5d0*3.141592653589793d0/180.d0
nsteps = floor(100.d0*2.d0*3.141592653589793d0/dphi)
Allocate(rout(1,nsteps+1))
Allocate(zout(1,nsteps+1))
Allocate(phiout(1,nsteps+1))

open(99,file='Lc_data.out')
Do j = 1,nphi
  Write(*,*) 'Working on j = ',j,' of ',nphi
  phistart(1) = phistarts(j)
  Do i = 1,Nstarts
!    Write(*,*) 'Working on i = ',i,' of ',Nstarts
    Rstart(1) = Rstarts(i)
    Zstart(1) = Zstarts(i)

!    write(*,*) 'starting line from',Rstart,Zstart,phistart
    Call ffr(Rstart,Zstart,phistart,1,dphi,nsteps,rout,zout,phiout,ierr_b,ilg)
    
    ithit = 0
    hit_jpt = ilg(1)
    Do jpt = 1,ilg(1)
      inside = inside_poly(rout(1,jpt),zout(1,jpt),g_lim(1,:),g_lim(2,:),g_limitr)
      if (inside .ne. 1) Then 
        ithit = 1
        hit_jpt = jpt
        Exit
      Endif
    Enddo
!    write(*,*) 'hit here',hit_jpt

    Lc(i,j) = Sum(sqrt( &
           (rout(1,1:hit_jpt-1)*cos(phiout(1,1:hit_jpt-1)) - rout(1,2:hit_jpt)*cos(phiout(1,2:hit_jpt)) )**2 &
         + (rout(1,1:hit_jpt-1)*sin(phiout(1,1:hit_jpt-1)) - rout(1,2:hit_jpt)*sin(phiout(1,2:hit_jpt)) )**2 &
         + (zout(1,1:hit_jpt-1) - zout(1,2:hit_jpt) )**2 &
           ))
    write(99,*) Lc(i,j)

  Enddo
Enddo
close(99)



End program test_footprint_rmp

