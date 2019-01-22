clearvars;

infile = 'C:\Work\FLARE\w7x\lc.dat';

gridfile ='C:\Work\FLARE\w7x\grid.dat';

griddata = dlmread(gridfile,'',2,0);
lcdata = dlmread(infile,'',10,0);

n = size(griddata,1);
x = griddata(:,1);
y = griddata(:,2);
z = griddata(:,3);

lc = lcdata(:,1) + lcdata(:,2);

% figure;
% plot3(x,y,z,'o')

figure;
ssize = ones(size(x))*200;
scatter3(x,y,z,ssize,lc,'.')