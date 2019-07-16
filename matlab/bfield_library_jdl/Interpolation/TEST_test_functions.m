clearvars;

nx = 10;
ny = 11;
[x,y] = meshgrid(linspace(0,1,nx),linspace(0,1,ny));



TESTNUM = 7;

f = eval_interp_test(x,y,TESTNUM);

figure; hold on; 
surf(x,y,f)
view(3)
colorbar;