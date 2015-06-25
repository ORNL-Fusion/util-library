function [j,ierr] = locate_bisect(xarr,x)
% !------------------------------------------------------------------------------
% !+ Search an ordered list of xarr for x, returns j where xarr(j) <= x <= xarr(j+1)
% !------------------------------------------------------------------------------
% !
% ! If x is off the table (including exactly equal to end values!) then
% !  j is set to 0 (off left end) or narr (off right end)
% !
ierr = 0;
jlo = 0;
narr = length(xarr);
jup = narr + 1;
test1 = ( xarr(narr) > xarr(1) );

while (jup - jlo > 1)
  jmid = floor((jup + jlo)/2);  %! compute midpoint
  if test1 == (x > xarr(jmid))
    jlo = jmid;
  else
    jup = jmid;
  end
end
j = jlo;

if ( j == 0 ) || (j == narr)
    ierr = 1;
end

