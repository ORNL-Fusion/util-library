!-----------------------------------------------------------------------------
!+ Math/geometric routines
!-----------------------------------------------------------------------------
Module math_geo_module
Use kind_mod
Implicit None
Contains

!------------------------------------------------------------------------------
!+ Quicksort
!------------------------------------------------------------------------------
Recursive Subroutine quicksort(array,narray, order)
! Order should be an integer array 1:narray
Use kind_mod
Integer(iknd), Intent(in) :: narray
Real(rknd), Intent(inout)  :: array(narray)
Integer(iknd), Intent(inout)  :: order(narray)

Integer(iknd) :: left, right, itemp, marker
Real(rknd) :: temp, pivot

If (narray .le. 1) Return

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
Use kind_mod
Implicit None
Logical, Intent(in) :: first
Integer(iknd), Intent(in) :: nline1, nline2
Real(rknd), Intent(in) :: rline1(nline1), zline1(nline1), rline2(nline2), zline2(nline2)
Real(rknd), Intent(out) :: pint1(2)
Integer(iknd), Intent(out) :: ierr, found_ind1, found_ind2, int_count

Integer(iknd) :: ii1, ii2
Real(rknd) :: p1(2), p2(2), p3(2), p4(2), dL1, dL2, dR1, dZ1, theta1, u1, u2

int_count = 0
found_ind1 = 0
found_ind2 = 0

Do ii1 = 1,nline1 - 1
  p1 = [rline1(ii1),zline1(ii1)]
  p2 = [rline1(ii1+1),zline1(ii1+1)]
  Do ii2 = 1,nline2 - 1
    p3 = [rline2(ii2),zline2(ii2)]
    p4 = [rline2(ii2+1),zline2(ii2+1)]
    Call int_two_lines(p1,p2,p3,p4,u1,u2)
    dL1 = Sqrt(Sum((p2-p1)**2))
    dL2 = Sqrt(Sum((p4-p3)**2))        
    If ( (Sign(1.d0,dL1) .ne. Sign(1.d0,u1)) .OR. (Sign(1.d0,dL2) .ne. Sign(1.d0,u2)) &
         .OR. (Abs(u1) .gt. 1.d0) .OR. (Abs(u2) .gt. 1.d0) ) Then
    Else
      dR1 = p2(1)-p1(1)
      dZ1 = p2(2)-p1(2)
      theta1 = Atan2(dZ1,dR1)
      pint1 = [dL1*u1*Cos(theta1)+p1(1),dL1*u1*Sin(theta1)+p1(2)]
      ierr = 0
      found_ind1 = ii1
      found_ind2 = ii2
      int_count = int_count + 1
      If (first) Return 
    Endif
  Enddo
Enddo

If (int_count .eq. 0) Then
    pint1 = [0.d0,0.d0]
    ierr = 1
Endif
If (int_count > 1) Then
    Write(*,*) 'Warning: More than one intersection found in int_line_curve. Returning last point'
Endif

End Subroutine int_curve_curve


!------------------------------------------------------------------------------
!+ Linear interpolation of curve C to find intersection with line L
!------------------------------------------------------------------------------
Subroutine int_line_curve(p1,p2,rline,zline,nline,first,pint1,ierr,found_ind,int_count)
! Curve is defined by array of points, L by p1,p2
! JL 2/2011
Use kind_mod
Implicit None
Logical, Intent(in) :: first
Integer(iknd), Intent(in) :: nline
Real(rknd), Intent(in) :: p1(2), p2(2), rline(nline), zline(nline)
Real(rknd), Intent(out) :: pint1(2)
Integer(iknd), Intent(out) :: ierr, found_ind, int_count

Integer(iknd) :: ii
Real(rknd) :: p3(2), p4(2), dL1, dL2, dR1, dZ1, theta1, u1, u2

int_count = 0
found_ind = 0
Do ii = 1,nline-1
  p3 = [rline(ii),zline(ii)]
  p4 = [rline(ii+1),zline(ii+1)]
  Call int_two_lines(p1,p2,p3,p4,u1,u2)
  dL1 = Sqrt(Sum((p2-p1)**2))
  dL2 = Sqrt(Sum((p4-p3)**2))
        
  If ( (Sign(1.d0,dL1) .ne. Sign(1.d0,u1)) .OR. (Sign(1.d0,dL2) .ne. Sign(1.d0,u2)) &
       .OR. (Abs(u1) .gt. 1.d0) .OR. (Abs(u2) .gt. 1.d0) ) Then
  Else
    dR1 = p2(1)-p1(1)
    dZ1 = p2(2)-p1(2)
    theta1 = Atan2(dZ1,dR1)
    pint1 = [dL1*u1*Cos(theta1)+p1(1),dL1*u1*Sin(theta1)+p1(2)]
    ierr = 0
    found_ind = ii
    int_count = int_count + 1
    If (first) Return      
  Endif  
Enddo

If (int_count .eq. 0) Then
    pint1 = [0.d0,0.d0]
    ierr = 1
Endif
If (int_count > 1) Then
    Write(*,*) 'Warning: More than one intersection found in int_line_curve. Returning last point'
Endif

End Subroutine int_line_curve


!-----------------------------------------------------------------------------
!+ Checks for the intersection of two line segments
!-----------------------------------------------------------------------------
Subroutine int_two_lines(p1,p2,p3,p4,u1,u2)
! Calculates the intersection of lines p1 to p2, and p3 to p4.  Each
! point is an x-y pair.  Returns two values, which are the distance
! between p1 and p2 of the intersection, and between p3 and p4.  This
! distance is normalized to the length of the lines.

Use kind_mod
implicit none

Real(rknd), Intent(in),Dimension(2) :: p1,p2,p3,p4
Real(rknd), intent(out) :: u1,u2
Real(rknd) :: denom

denom = (p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2))

if ( denom .eq. 0._rknd ) then 
  u1 = 1.d30
  u2 = 1.d30
else
  u1 = ((p4(1)-p3(1))*(p1(2)-p3(2)) - (p4(2)-p3(2))*(p1(1)-p3(1)))/denom
  u2 = ((p2(1)-p1(1))*(p1(2)-p3(2)) - (p2(2)-p1(2))*(p1(1)-p3(1)))/denom
endif

Endsubroutine int_two_lines

!------------------------------------------------------------------------------
!+ Move a distance on a curve linear - piecewise
!------------------------------------------------------------------------------
Subroutine move_L_on_C(L,rline,zline,nline,ic_near_L,err_near_L,R_L,Z_L)
! Moves a distance L along the curve defined by
! arrays rline,zline.  dL is defined as the
! linear distance between curve points.
! JL 2/2011
Use kind_mod
Implicit None
Integer(iknd), Intent(in) :: nline
Real(rknd), Intent(in) :: L, rline(nline), zline(nline)
Real(rknd), Intent(out) :: err_near_L, R_L, Z_L
Integer(iknd), Intent(out) :: ic_near_L

Real(rknd) :: dL(nline), Ltot, SumL(nline), delt(nline), diff, theta, dL_step
Integer(iknd) :: ind_diff, ii, i

dL = 0.d0
dL(2:nline) = Sqrt( (rline(1:nline-1)-rline(2:nline))**2 + (zline(1:nline-1)-zline(2:nline))**2 )

Ltot = Sum(dL)
If ( Ltot .lt. L ) Then
    Write(*,*) 'Requested length, total length ',L,Ltot
    Write(*,*) 'Error: Ltot < L'
    Stop
Endif

SumL(1) = 0.d0
Do i = 2,nline
  SumL(i) = SumL(i-1) + dL(i)
Enddo
delt = L - SumL
ind_diff = Minloc(Abs(delt),1)
diff = delt(ind_diff)

If ( diff .gt. 0.d0 ) Then
    dL_step = dL(ind_diff+1)
    ii = ind_diff
    theta = Atan2(zline(ii+1)-zline(ii),rline(ii+1)-rline(ii))
    If ( dL_step .eq. 0.d0 ) Then
      Write(*,*) ' dL_step 1 cannot be zero in move_L_on_C',dL_step
      Stop
    Endif
    R_L = rline(ii) + diff*Cos(theta)
    Z_L = zline(ii) + diff*Sin(theta)
Elseif ( diff .lt. 0.d0 ) Then
    dL_step = dL(ind_diff)
    ii = ind_diff
    theta = Atan2(zline(ii-1)-zline(ii),rline(ii-1)-rline(ii))
    If ( dL_step == 0.d0 ) Then
      Write(*,*) ' dL_step 2 cannot be zero in move_L_on_C',dL_step,ind_diff
      Stop
    Endif
    R_L = rline(ii) - diff*Cos(theta)
    Z_L = zline(ii) - diff*Sin(theta)
Else
    R_L = rline(ind_diff)
    Z_L = zline(ind_diff)
Endif

ic_near_L = ind_diff
err_near_L = diff

End Subroutine move_L_on_C

!------------------------------------------------------------------------------
!+ Linear interpolation
!------------------------------------------------------------------------------
Subroutine linear_interp(xarr,yarr,narr,xin,yout,ierr)
Use kind_mod
Implicit None
Integer(iknd), Intent(In) :: narr
Real(rknd), Intent(In) :: xarr(narr), yarr(narr), xin
Integer(iknd), Intent(Out) :: ierr
Real(rknd), Intent(Out) :: yout

Integer(iknd) :: j

Call locate_bisect(xarr,narr,xin,j,ierr)

If (ierr .eq. 1) Then 
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
Use kind_mod
Implicit None
Integer(iknd), Intent(In) :: narr
Real(rknd), Intent(In) :: xarr(narr), x
Integer(iknd), Intent(Out) :: j, ierr
Logical :: test1
Integer(iknd) :: jlo, jup, jmid
ierr = 0
jlo = 0
jup = narr + 1
test1 = ( xarr(narr) .gt. xarr(1) )

!If ( ( x .gt. Maxval(xarr) ) .OR. ( x .lt. Minval(xarr) ) ) Then
!  ierr = 1
!  j = 0
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
Use kind_mod                ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Real(rknd),    Intent(in) :: xstart  
Real(rknd),    Intent(in) :: xend
Integer(iknd), Intent(in) :: numel
Real(rknd)                :: rlinvec(numel)

! Local scalars
Integer(iknd)   ::  ii
!- End of header -------------------------------------------------------------

Do ii = 1,numel
  rlinvec(ii) = ( Real(ii,rknd) - 1._rknd ) * ( xend - xstart ) &
       / ( numel - 1._rknd ) + xstart
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

Use kind_mod
Use phys_const, Only : pi
Implicit none

Real(rknd), Intent(in) :: x,y
Integer(iknd), Intent(in) :: npoly
Real(rknd), Dimension(npoly), Intent(in) :: px,py

Integer(iknd) :: inside

! Local variables
Real(rknd), Dimension(npoly+1) :: tmp_px, tmp_py
Real(rknd), Dimension(npoly) :: theta,dp,cp

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
