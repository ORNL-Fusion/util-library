!
!  routine to read the bfield from Geiger/w7's xdr output file
!  and dump the data to an ascii file.
!
! JDL 7/13/2011 (first part modified from Geiger's(?) readbgrid)
!

program dump_xdr_to_ascii

use fxdr_module
implicit none

character(len=9) :: filename_in
character(len=10) :: filename_out
logical :: l_val=.false., l_check=.false., xdr_check 
integer :: ierr, ixdrs
integer :: i1, i2
integer :: iz,ir,i
integer :: nperio,ialfa,knull,iganz
integer :: k2, iald21
real(kind(1.0d+0)) :: bmean
real(kind(1.d+0)) :: rnull,ronull,eta,fnull,ampby0,bz0,bfak
real(kind(1.d+0)), dimension(:,:,:), allocatable :: brg
real(kind(1.d+0)), dimension(:,:,:), allocatable :: bfg
real(kind(1.d+0)), dimension(:,:,:), allocatable :: bzg



xdr_check = .false.
filename_in = 'input.xdr'
filename_out = 'output.dat'


ixdrs = initxdr(filename_in,'r',l_val)
ierr = ixdrint(ixdrs,i1)
ierr = ixdrint(ixdrs,i2)
ierr = ixdrdouble(ixdrs,rnull)
ierr = ixdrdouble(ixdrs,ronull)
ierr = ixdrdouble(ixdrs,eta)
ierr = ixdrdouble(ixdrs,fnull)
ierr = ixdrint(ixdrs,nperio)
ierr = ixdrint(ixdrs,ialfa)
ierr = ixdrint(ixdrs,knull)
ierr = ixdrint(ixdrs,iganz)
ierr = ixdrdouble(ixdrs,ampby0)
ierr = ixdrdouble(ixdrs,bz0)
ierr = ixdrdouble(ixdrs,bfak)
k2=2*knull+1
iald21=ialfa/2+1
write(6,*) "Parameters read:"
write(6,*) "i1 : ", i1
write(6,*) "i2 : ", i2
write(6,*) "rnull : ", rnull
write(6,*) "ronull: ", ronull
write(6,*) "eta   : ", eta
write(6,*) "fnull : ", fnull
write(6,*) "nperio: ", nperio
write(6,*) "ialfa : ", ialfa
write(6,*) "knull : ", knull
write(6,*) "iganz : ", iganz
write(6,*) "ampby0 : ", ampby0
write(6,*) "bz0 : ", bz0
write(6,*) "bfak : ", bfak
write(6,*) "Number of R-grid points: ",k2
write(6,*) "Number of z-grid points: ",k2

allocate(brg(iald21,k2,k2),bfg(iald21,k2,k2),bzg(iald21,k2,k2))
brg = 0.d+0 ; bfg = 0.d+0 ; bzg = 0.d+0
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      ierr = ixdrdouble(ixdrs,bfg(i,ir,iz))
    enddo
  enddo
enddo
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      ierr = ixdrdouble(ixdrs,brg(i,ir,iz))
    enddo
  enddo
enddo
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      ierr = ixdrdouble(ixdrs,bzg(i,ir,iz))
    enddo
  enddo
enddo
ierr = ixdrclose(ixdrs)
write(6,*) "End of reading!"

if (xdr_check .eqv. .true.) then
do i=1,iald21
  do iz=2,k2-1
    do ir=2,k2-1
      if( bfg(i,ir-1,iz) /= 0.d+0 .and. &
          bfg(i,ir+1,iz) /= 0.d+0 .and. &
          bfg(i,ir,iz-1) /= 0.d+0 .and. &
          bfg(i,ir,iz+1) /= 0.d+0 ) then
        bmean=0.25d+0*(bfg(i,ir-1,iz)+bfg(i,ir+1,iz)+ &
                       bfg(i,ir,iz-1)+bfg(i,ir,iz+1))
        if( abs(bfg(i,ir,iz)) < 0.5d0*abs(bmean) .or. &
            abs(bfg(i,ir,iz)) > 1.5d0*abs(bmean) ) then
          write(6,*)"Point detected:"
          write(6,*)"bfg (",i,",",ir,",",iz,") = ",bfg(i,ir,iz)
          write(6,*)"bmean(",i,",",ir,",",iz,") = ",bmean
          brg(i,ir,iz)=0.25d+0*(brg(i,ir-1,iz)+brg(i,ir+1,iz)+ &
                       brg(i,ir,iz-1)+brg(i,ir,iz+1))
          bfg(i,ir,iz)=0.25d+0*(bfg(i,ir-1,iz)+bfg(i,ir+1,iz)+ &
                       bfg(i,ir,iz-1)+bfg(i,ir,iz+1))
          bzg(i,ir,iz)=0.25d+0*(bzg(i,ir-1,iz)+bzg(i,ir+1,iz)+ &
                       bzg(i,ir,iz-1)+bzg(i,ir,iz+1))
          l_check=.true.
        endif
      endif
    enddo
  enddo
enddo
if(l_check) then
  write(*,*) '----------------------------------'
  write(*,*) 'Bad points have been "corrected"'
  write(*,*) 'Instead of rewriting the xdr file'
  write(*,*) 'I am just dumping the corrected data'
  write(*,*) '7-13-2011 JDL'
  write(*,*) '----------------------------------'
endif
else
write(*,*) 'NOT CHECKING FOR BAD POINTS!!!!!!!!'
endif

Open(12, file=filename_out,Position="asis",status="unknown")

write(*,*) 'Writing file ',filename_out
write(12,*) i1
write(12,*) i2
write(12,*) rnull
write(12,*) ronull
write(12,*) eta
write(12,*) fnull
write(12,*) nperio
write(12,*) ialfa
write(12,*) knull
write(12,*) iganz
write(12,*) ampby0
write(12,*) bz0
write(12,*) bfak
write(12,*) k2
write(12,*) iald21
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      write(12,*) brg(i,ir,iz)
    enddo
  enddo
enddo
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      write(12,*) bfg(i,ir,iz)
    enddo
  enddo
enddo
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      write(12,*) bzg(i,ir,iz)
    enddo
  enddo
enddo
!write(12,*) brg   ! br
!write(12,*) bfg   ! bphi
!write(12,*) bzg   ! bz
Close(12)

write(*,*) '>>>>>>done writing file<<<<<'





end program dump_xdr_to_ascii
