function [z,ierr] = linear_interp_2d(xarr,yarr,z2d,x,y)
% bilinearly interpolate z = f(x,y) given z2d = f(xarr,yarr)
% Add: error handling

if x < min(xarr) || x > max(xarr)
    error('x out of bounds')
end
if y < min(yarr) || y > max(yarr)
    error('y out of bounds')
end

[j,ierr] = locate_bisect(xarr,x);
if ierr ~= 0
    error('Error from locate_bisect in x')
end
[k,ierr] = locate_bisect(yarr,y);
if ierr ~= 0
    error('Error from locate_bisect in y')
end

t = (x - xarr(j))/(xarr(j+1)-xarr(j));
u = (y - yarr(k))/(yarr(k+1)-yarr(k));
z = (1-t)*(1-u)*z2d(j,k) + t*(1-u)*z2d(j+1,k) + t*u*z2d(j+1,k+1) + (1-t)*u*y(j,k+1);

ierr = 0;