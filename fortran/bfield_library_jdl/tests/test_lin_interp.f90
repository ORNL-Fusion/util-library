Program test_lin_interp
Use kind_mod
Use math_geo_module
Implicit None

Integer(iknd), Parameter :: narr = 5
Real(rknd) :: xarr(narr),yarr(narr), x, yout
Integer(iknd) :: j, ierr

!xarr = rlinspace(0.d0,-1.d0,narr)
xarr = rlinspace(0.d0,-1.d0,narr)
yarr = rlinspace(0.d0,-10.d0,narr)
yarr(narr) = 1000.d0
x = -.99999d0

Call linear_interp(xarr,yarr,narr,x,yout,ierr)


Write(*,*) 'xarr',xarr
Write(*,*) 'yarr',yarr
Write(*,*) 'Interpolate at x = ',x
Write(*,*) 'Answer ',yout
Write(*,*) 'Error? ',ierr




End program test_lin_interp
