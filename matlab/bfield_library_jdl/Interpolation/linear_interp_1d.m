function [y,ierr] = linear_interp_1d(xarr,yarr,x)
% Linear interpolate y = f(x) given yarr = f(xarr)
% Add: error handling

if abs(x - xarr(1)) < eps
    y = yarr(1);
    ierr = 0;
    return;
end

if abs(x - xarr(end)) < eps
    y = yarr(end);
    ierr = 0;
    return;
end

if x < min(xarr) || x > max(xarr)
    error('x out of bounds')
end

[j,ierr] = locate_bisect(xarr,x);

if ierr ~= 0
    error('Error from locate_bisect')
end

t = (x - xarr(j))/(xarr(j+1)-xarr(j));
y = (1-t)*yarr(j) + t*yarr(j+1);

ierr = 0;