clearvars;


xarr = linspace(10,1,10);

yarr = xarr.^2;

x = [5.5,6.5]


[y,ierr] = linear_interp_1d(xarr,yarr,x)