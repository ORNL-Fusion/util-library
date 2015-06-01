Program test_quicksort
Use kind_mod
Use math_geo_module
Implicit None

Integer(iknd), Parameter :: narr = 8
Real(rknd) :: xarr(narr), xarr_org(narr)
Integer(iknd) :: j, order(narr)

xarr = [-.2d0,5.d0,-5.d0,0.d0,-5.d0,0.d0,0.d0,-0.d0]
Do j = 1,narr
  order(j) = j
Enddo
xarr_org = xarr
Write(*,*) 'Original array:',xarr_org
Call quicksort(xarr,narr,order)
Write(*,*) 'Sorted array  :',xarr
Write(*,*) 'Order',order
Write(*,*) 'Ordered original',xarr_org(order)





End program test_quicksort
