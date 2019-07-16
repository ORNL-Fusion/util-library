clearvars;

TESTNUM = 1;

% For original definition
nx = 10;
ny = 11;

% For interpolation
nxf = 50;
nyf = 51;

% define evaluation points
x1d = linspace(0,1,nx);
y1d = linspace(0,1,ny);
[x,y] = meshgrid(x1d,y1d);
xf1d = linspace(0,1,nxf);
yf1d = linspace(0,1,nyf);
[xf,yf] = meshgrid(xf1d,yf1d);

% Evaluate test function
z = eval_interp_test(x,y,TESTNUM);
zf_0 = eval_interp_test(x,y,TESTNUM);

% Interpolate
% interp2d_method = 'nearest';
interp2d_method = 'linear';
% interp2d_method = 'cubic';
% interp2d_method = 'makima';
% zf = interp2(x,y,z,xf,yf,interp2d_method);
[zf,ierr] = linear_interp_2d(x1d,y1d,z,xf1d,yf1d);

figure; set(gcf,'pos',[621 678 1117 420])
subplot(1,2,1); hold on; view(102,22); colorbar; title('input');
surf(x,y,z)
subplot(1,2,2); hold on; view(102,22); colorbar; title('interp2')
surf(xf,yf,zf)
% view(3)
