!-----------------------------------------------------------------------------
!+ Math/geometric routines
! 
!  Subroutine quicksort
!  Subroutine int_curve_curve
!  Subroutine int_line_curve
!  Subroutine int_two_lines
!  Subroutine move_L_on_C
!  Subroutine linear_interp
!  Subroutine locate_bisect
!  Function rlinspace
!  Function inside_poly
!
!-----------------------------------------------------------------------------
Module math_geo_module
Implicit None
Private
Public :: quicksort
Public :: int_curve_curve
Public :: int_line_curve
Public :: int_two_lines
Public :: move_L_on_C
Public :: linear_interp
Public :: rlinspace
Public :: inside_poly
Contains

!------------------------------------------------------------------------------
!+ Quicksort
!------------------------------------------------------------------------------
Recursive Subroutine quicksort(array,narray, order)
  ! order must be an integer "index" array. order=[1:narray]
  ! When finished, order is an index array that will order the
  ! original array. I.e., array_orig(order) = array_sorted
  ! Also, array2(order) = array_sorted will equal array_orig
Use kind_mod, Only: real64, int32
Integer(int32), Intent(in) :: narray
Real(real64), Intent(inout)  :: array(narray)
Integer(int32), Intent(inout)  :: order(narray)

Integer(int32) :: left, right, itemp, marker
Real(real64) :: temp, pivot

If (narray .le. 1) Return
!! QQ check if order(1) = 1, (2)=2, else define it here
left = 0
right = narray + 1
pivot = array((left+right)/2)  ! could replace by random number

! Move pointers until they cross
Do While (left < right) 
  ! Find elements on either side of pivot
  right = right - 1
  Do While (array(right) > pivot)
    right = right - 1
  Enddo
  left = left + 1
  Do While (array(left) < pivot)
    left = left + 1
  Enddo
  ! Swap out of place elements
  If (left < right) Then
     temp = array(left); array(left) = array(right); array(right) =  temp
    itemp = order(left); order(left) = order(right); order(right) = itemp
  Endif
Enddo

! If you ended on the pivot then it is in place, repeat with sub arrays
If (left .eq. right) Then
  marker = left + 1
Else
  marker = left
Endif

Call quicksort(array(1:marker-1),marker-1,order(1:marker-1))
Call quicksort(array(marker:narray),narray-marker+1,order(marker:narray))
End Subroutine quicksort

!------------------------------------------------------------------------------
!+ Linear interpolation of two curves C
!------------------------------------------------------------------------------
Subroutine int_curve_curve(rline1,zline1,nline1,rline2,zline2,nline2,first,pint1,ierr,found_ind1,found_ind2,int_count)
! First curve is stepped and both curves linearly interpolated
! Curve is defined by array of points
! JL 2/2011
Use kind_mod, Only: real64, int32
Implicit None
Logical, Intent(in) :: first
Integer(int32), Intent(in) :: nline1, nline2
Real(real64), Intent(in) :: rline1(nline1), zline1(nline1), rline2(nline2), zline2(nline2)
Real(real64), Intent(out) :: pint1(2)
Integer(int32), Intent(out), Optional :: ierr, found_ind1, found_ind2, int_count

Integer(int32) :: ii1, ierr_tmp, int_count_tmp
Integer(int32) :: ierr_0, found_ind1_0, found_ind2_0, int_count_0
Real(real64) :: p1(2), p2(2)
Real(real64) :: pint1_tmp(2)

int_count_0 = 0
found_ind1_0 = 0
found_ind2_0 = 0

Do ii1 = 1,nline1 - 1
  p1 = [rline1(ii1),zline1(ii1)]
  p2 = [rline1(ii1+1),zline1(ii1+1)]
  Call int_line_curve(p1,p2,rline2,zline2,nline2,first,pint1_tmp,ierr_tmp,found_ind2_0,int_count_tmp)

  If (ierr_tmp == 0) Then
    ierr_0 = 0
    pint1 = pint1_tmp
    found_ind1_0 = ii1
    int_count_0 = int_count_0 + int_count_tmp

    If (first) Then
      If (Present(ierr)) ierr = ierr_0
      If (Present(found_ind1)) found_ind1 = found_ind1_0
      If (Present(found_ind2)) found_ind2 = found_ind2_0
      If (Present(int_count)) int_count = int_count_0
      Return
    Endif
  Endif
Enddo

If (int_count_0 .eq. 0) Then
    pint1 = [0.d0,0.d0]
    ierr_0 = 1
    If (Present(ierr)) ierr = ierr_0
    If (Present(found_ind1)) found_ind1 = 0
    If (Present(found_ind2)) found_ind2 = 0
    If (Present(int_count)) int_count = int_count_0    
Endif
If (int_count_0 > 1) Then
    Write(*,*) 'Warning: More than one intersection found in int_curve_curve. Returning last point'
Endif

End Subroutine int_curve_curve


!------------------------------------------------------------------------------
!+ Linear interpolation of curve C to find intersection with line L
!------------------------------------------------------------------------------
Subroutine int_line_curve(p1,p2,rline,zline,nline,first,pint1,ierr,found_ind,int_count_out)
! Curve is defined by array of points, L by p1,p2
! JL 2/2011
Use kind_mod, Only: real64, int32
Implicit None
Logical, Intent(in) :: first
Integer(int32), Intent(in) :: nline
Real(real64), Intent(in) :: p1(2), p2(2), rline(nline), zline(nline)
Real(real64), Intent(out) :: pint1(2)
Integer(int32), Intent(out), Optional :: ierr, found_ind, int_count_out
Logical :: test
Integer(int32) :: ii, ierr2, int_count
Real(real64) :: p3(2), p4(2), u1, u2

int_count = 0
If (Present(found_ind)) found_ind = 0
Do ii = 1,nline-1
  p3 = [rline(ii),zline(ii)]
  p4 = [rline(ii+1),zline(ii+1)]
  Call int_two_lines(p1,p2,p3,p4,u1,u2,ierr2)

  If (ierr2 .ne. 0) Cycle

  If (ii == nline -1) Then
    test = ((u1 .ge. 0._real64) .AND. (u1 .le. 1._real64) .AND. &
         (u2 .ge. 0._real64) .AND. (u2 .le. 1._real64))
  Else
    test = ((u1 .ge. 0._real64) .AND. (u1 .le. 1._real64) .AND. &
         (u2 .ge. 0._real64) .AND. (u2 .lt. 1._real64))
  Endif

  If (test) Then
    int_count = int_count + 1
    pint1 = p1 + u1*(p2-p1)
    If(Present(ierr)) ierr = 0
    If (Present(found_ind)) found_ind = ii
    If (first) Then
      If (Present(int_count_out)) int_count_out = int_count
      Return
    Endif
  Endif  
Enddo

If (int_count .eq. 0) Then
    pint1 = [0.d0,0.d0]
    If(Present(ierr)) ierr = 1
Endif
If (int_count > 1) Then
    Write(*,*) 'Warning: More than one intersection found in int_line_curve. Returning last point'
Endif
If (Present(int_count_out)) int_count_out = int_count

End Subroutine int_line_curve


!-----------------------------------------------------------------------------
!+ Checks for the intersection of two line segments
!-----------------------------------------------------------------------------
Subroutine int_two_lines(p1,p2,p3,p4,u1,u2,ierr,pint)
! Calculates the intersection of lines p1 to p2, and p3 to p4.  Each
! point is an x-y pair.  Returns two values, which are the distance
! between p1 and p2 of the intersection, and between p3 and p4.  This
! distance is normalized to the length of the lines.

Use kind_mod, Only: real64, Int32
implicit none

Real(real64), Intent(in),Dimension(2) :: p1,p2,p3,p4
Real(real64), Intent(out) :: u1,u2
Real(real64), Intent(out), Optional :: pint(2)
Real(real64) :: denom
Real(real64), Parameter :: tol = 2.e-16_real64
Integer(int32), Intent(out) :: ierr

denom = (p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2))

If ( Abs(denom) .lt. tol ) then  ! parallel lines
  u1=0._real64
  u2=0._real64
  If (Present(pint)) pint = 0._real64
  ierr=1
Else
  u1 = ((p4(1)-p3(1))*(p1(2)-p3(2)) - (p4(2)-p3(2))*(p1(1)-p3(1)))/denom
  u2 = ((p2(1)-p1(1))*(p1(2)-p3(2)) - (p2(2)-p1(2))*(p1(1)-p3(1)))/denom
  ierr=0
  If (Present(pint)) pint = p1 + u1*(p2-p1)
Endif

Endsubroutine int_two_lines

!------------------------------------------------------------------------------
!+ Move a distance on a curve linear - piecewise
!------------------------------------------------------------------------------
Subroutine move_L_on_C(L,rline,zline,nline,ic_near_L,err_near_L,R_L,Z_L,ierr)
! Moves a distance L along the curve defined by
! arrays rline,zline.  dL is defined as the
! linear distance between curve points.
! JL 2/2011
Use kind_mod, Only: real64, int32
Implicit None
Integer(int32), Intent(in) :: nline
Real(real64), Intent(in) :: L, rline(nline), zline(nline)
Real(real64), Intent(out) :: R_L, Z_L
Real(real64), Intent(out), Optional :: err_near_L
Integer(int32), Intent(out), Optional :: ierr, ic_near_L

Real(real64) :: dL(nline-1), Ltot, SumL, f
Integer(int32) :: i


dL = Sqrt( (rline(1:nline-1)-rline(2:nline))**2 + (zline(1:nline-1)-zline(2:nline))**2 )

Ltot = Sum(dL)
If ( Ltot .lt. L ) Then
    Write(*,*) 'Requested length, total length ',L,Ltot
    Write(*,*) 'Error: Ltot < L'
    If (Present(ierr)) ierr = 1
    Return
Endif

SumL = 0.d0
Do i = 1,nline-1
  If (SumL + dL(i) .gt. L) Exit
  SumL = SumL + dL(i)
Enddo
f = (L - SumL)/dL(i)
R_L = f*(rline(i+1)-rline(i))+rline(i)
Z_L = f*(zline(i+1)-zline(i))+zline(i)
If (f .gt. 0.5_real64) Then
  If (Present(ic_near_L)) ic_near_L = i+1
  If (Present(err_near_L)) err_near_L = (1._real64 - f)*dL(i)
Else
  If (Present(ic_near_L)) ic_near_L = i
  If (Present(err_near_L)) err_near_L = f*dL(i)
Endif
If (Present(ierr)) ierr = 0
End Subroutine move_L_on_C

!------------------------------------------------------------------------------
!+ Linear interpolation
!------------------------------------------------------------------------------
Subroutine linear_interp(xarr,yarr,narr,xin,yout,ierr)
Use kind_mod, Only: real64, int32
Implicit None
Integer(int32), Intent(In) :: narr
Real(real64), Intent(In) :: xarr(narr), yarr(narr), xin
Integer(int32), Intent(Out) :: ierr
Real(real64), Intent(Out) :: yout

Integer(int32) :: j

If (xin .lt. Minval(xarr)) Then
  Write(*,*) 'Error from linear_interp, x out of bounds lower (x,Minval(xarr))=',xin,Minval(xarr)
  ierr = 1
  Return
!  Stop "quitting"
Endif
If (xin .gt. Maxval(xarr)) Then
  Write(*,*) 'Error from linear_interp, x out of bounds upper (x,Maxval(xarr))=',xin,Maxval(xarr)
  ierr = 1
  Return
!  Stop "quitting"
Endif
Call locate_bisect(xarr,narr,xin,j,ierr)

If (j .eq. 0) Then
  If (abs(xin - xarr(1)) < 1.e-10) Then
    yout = yarr(1)
    ierr = 0
    Return
  Endif
Endif
If (j .eq. narr) Then
  If (abs(xin - xarr(narr)) < 1.e-10) Then
    yout = yarr(narr)
    ierr = 0
    Return
  Endif
Endif

If (ierr .eq. 1) Then
  Write(*,*) "xarr([1,end])",xarr(1),xarr(narr)
  Write(*,*) "xin",xin
  Write(*,*) 'xin - xarr(1)',xin-xarr(1)
  Stop "error linear interp!, ierr set from locate_bisect"  
  yout = 0.d0
  Return
Endif

yout = yarr(j) + ((xin-xarr(j))/(xarr(j+1)-xarr(j)))*(yarr(j+1)-yarr(j))

End Subroutine linear_interp



!------------------------------------------------------------------------------
!+ Search an ordered list of xarr for x, returns j where xarr(j) <= x <= xarr(j+1)
!------------------------------------------------------------------------------
Subroutine locate_bisect(xarr,narr,x,j,ierr)
!
! If x is off the table (including exactly equal to end values!) then
!  j is set to 0 (off left end) or narr (off right end)
!
Use kind_mod, Only: real64, int32
Implicit None
Integer(int32), Intent(In) :: narr
Real(real64), Intent(In) :: xarr(narr), x
Integer(int32), Intent(Out) :: j, ierr
Logical :: test1
Integer(int32) :: jlo, jup, jmid
ierr = 0
jlo = 0
jup = narr + 1
test1 = ( xarr(narr) .gt. xarr(1) )

!If ( ( x .gt. Maxval(xarr) ) .OR. ( x .lt. Minval(xarr) ) ) Then
!  ierr = 1
!  Stop "x out of range in locate_bisect"
!Endif

Do While (jup - jlo .gt. 1) 
  jmid = (jup + jlo)/2  ! compute midpoint
  If ( test1 .eqv. (x .gt. xarr(jmid)) ) Then
    jlo = jmid
  Else
    jup = jmid
  Endif
Enddo
j = jlo

If (( j .eq. 0 ) .OR. (j .eq. narr)) ierr = 1

End Subroutine locate_bisect

!------------------------------------------------------------------------------
!+ Returns a linearly spaced real vector given endpoints and number of elements
!------------------------------------------------------------------------------
Function rlinspace(xstart,xend,numel)  & 
Result(rlinvec)
!
! Description: 
!   This function returns a real vector of length(numel) with linearly spaced
!   values from xstart to xend.  Similar to the Matlab function.
!
! Inputs:
!  xstart,xend: Values of the first and last points of the array [real]
!  numel: Number of elements in the array [integer]
! Outputs:
!  rlinvec: The linearly spaced array
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/22/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3. JL
!  1.2     07/18/2011  Ported to w7 routines. JL
!
! Author(s): J. Lore 07/2009 - 7/18/2011
!
! Modules used:
Use kind_mod, Only: real64, int32

Implicit None

! Input/output                      !See above for descriptions
Real(real64),    Intent(in) :: xstart  
Real(real64),    Intent(in) :: xend
Integer(int32), Intent(in) :: numel
Real(real64)                :: rlinvec(numel)

! Local scalars
Integer(int32)   ::  ii
!- End of header -------------------------------------------------------------

If (numel .eq. 1) Then
  rlinvec(1) = xstart
  Return
Endif

Do ii = 1,numel
  rlinvec(ii) = ( Real(ii,real64) - 1._real64 ) * ( xend - xstart ) &
       / ( numel - 1._real64 ) + xstart
Enddo

EndFunction rlinspace


Function inside_poly(x,y,px,py,npoly) &
Result(inside)
!
! Inputs:
!   npoly - size of vectors px and py (integer)
!   x - The x coordinate of the point.
!   y - The y coordinate of the point.
!   px - The x coordinates of the polygon.
!   py - The y coordinates of the polygon.
! 
!  The return value of the function is 1 if the point is inside the
!  polygon and 0 if it is outside the polygon.

!
! Polygon is closed by connecting first and last points. If polygon is already
! closed this will not affect the results (theta contribution is zero)  
!

Use kind_mod, Only: real64, int32
Use phys_const, Only : pi
Implicit none

Real(real64), Intent(in) :: x,y
Integer(int32), Intent(in) :: npoly
Real(real64), Dimension(npoly), Intent(in) :: px,py

Integer(int32) :: inside

! Local variables
Real(real64), Dimension(npoly+1) :: tmp_px, tmp_py
Real(real64), Dimension(npoly) :: theta,dp,cp

!- End of header -------------------------------------------------------------

! Close polygon
tmp_px(1:npoly) = px
tmp_px(npoly+1) = px(1)
tmp_py(1:npoly) = py
tmp_py(npoly+1) = py(1)

dp = (tmp_px(1:npoly)-x)*(tmp_px(2:npoly+1)-x) + (tmp_py(1:npoly)-y)*(tmp_py(2:npoly+1)-y) ! dot product
cp = (tmp_px(1:npoly)-x)*(tmp_py(2:npoly+1)-y) - (tmp_py(1:npoly)-y)*(tmp_px(2:npoly+1)-x) ! cross product

theta = atan2(cp,dp)

inside = 0
If ( Abs(Sum(theta)) .gt. pi ) Then
  inside = 1
Endif

Endfunction inside_poly

End Module math_geo_module
