!-----------------------------------------------------------------------------
!
! Tests footprints
! -- JDL
!
Program test_footprint_screened

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
Use screening_module
Use math_geo_module

Implicit None

! Local variables (scalar)
Real(rknd),Allocatable :: Rstarts(:),Zstarts(:),phistarts(:),Lc(:,:)
Real(rknd) :: Rstart(1),Zstart(1),phistart(1)
Integer(iknd) :: nphi_fp, Nstarts, i, j, nsteps,ierr_b(1), ilg(1), inside,ithit,hit_jpt,jpt

Character(Len=100) :: gfilename

Integer :: ntest,ntorpts
Real(rknd) :: dphi, taper(6)
Real(rknd),Allocatable,Dimension(:,:) :: rout,zout,phiout
Character(Len=200) :: Afile_path
!- End of header -------------------------------------------------------------


!gfilename = '../bfield_files/g135183.00433'
gfilename = '../bfield_files/g140508.00403'
Call readg_g3d(gfilename)

Afile_path = '/home/jjl/Pavel/screened-8-19/'
Call setup_screening_vars(Afile_path)
Call read_Afiles(Afile_path)
Call prepare_Afile_splinefits !(Afile_path)

! footprint settings
Nstarts = 400
nphi_fp = 500
Allocate(Rstarts(Nstarts),Zstarts(Nstarts))
Allocate(phistarts(nphi_fp))
Allocate(Lc(Nstarts,nphi_fp))
Rstarts = rlinspace(0.35d0,0.5d0,Nstarts)
Zstarts = -1.603d0
phistarts = rlinspace(0.d0,2.0944d0,nphi_fp)

! fieldline following settings
rmp_method = 2
dphi = 0.5d0*3.141592653589793d0/180.d0
nsteps = floor(100.d0*2.d0*3.141592653589793d0/dphi)
Allocate(rout(1,nsteps+1))
Allocate(zout(1,nsteps+1))
Allocate(phiout(1,nsteps+1))

open(99,file='Lc_data.out')
Do j = 1,nphi_fp
  Write(*,*) 'Working on j = ',j,' of ',nphi_fp
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



End program test_footprint_screened

