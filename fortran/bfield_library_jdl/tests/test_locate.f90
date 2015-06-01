Program test_locate
Use kind_mod
Use math_geo_module
Implicit None

Integer(iknd), Parameter :: narr = 5
Real(rknd) :: xarr(narr), x
Integer(iknd) :: j, ierr

xarr = rlinspace(0.d0,-1.d0,narr)
x = -1.0d0

Call locate_bisect(xarr,narr,x,j,ierr)


Write(*,*) 'xarr',xarr
If (ierr .eq. 0) Then
Write(*,*) 'x,j,xarr(j:j+1)',x,j,xarr(j),xarr(j+1)
Endif
Write(*,*) 'Error? (ierr,j)',ierr, j




End program test_locate
