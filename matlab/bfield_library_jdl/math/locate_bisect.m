function [j,ierr] = locate_bisect(xarr,x)
% ------------------------------------------------------------------------------
%  Search an ordered list of xarr for x, returns j where xarr(j) <= x = xarr(j+1)
%  Goal is to return the "left" index, which should only be violated when 
%  ierr = 1. Array must be monotonic.
% ------------------------------------------------------------------------------
% 
% If x is off the table (including exactly equal to end points) then
%  j is set to 0 (off left end) or narr (off right end) and ierr = 1
%
ierr = 0;
jlo = 0;
narr = length(xarr);
jup = narr + 1;
is_ascend = ( xarr(narr) > xarr(1) );  % Check if ascending

while (jup - jlo > 1)
  jmid = floor((jup + jlo)/2);  % compute midpoint
  if is_ascend == (x > xarr(jmid))
    jlo = jmid;
  else
    jup = jmid;
  end
end
j = jlo;

if ( j == 0 ) || (j == narr)
    ierr = 1;
end

