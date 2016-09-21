clearvars;

% run_path = 'C:\Work\fortran\test_poincare\164723\AS_2d\'; gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
run_path = 'C:\Work\fortran\test_poincare\164723\VAC_2d\'; gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
g = readg_g3d(gfile_name);

fname = fullfile(run_path,'poincare_output_2d.out');

data = dlmread(fname);

phi = data(1,1);
n = data(1,2);
n2 = sqrt(n);

ind = data(2:end,1);
r = data(2:end,2);
z = data(2:end,3);
psimin = data(2:end,4);
Lc = data(2:end,5);

r2 = reshape(r,[n2,n2]);
z2 = reshape(z,[n2,n2]);
psimin2 = reshape(psimin,[n2,n2]);
Lc2d = reshape(Lc,[n2,n2]);

figure; hold on;
plot(r,z,'k.')
h=pcolor(r2,z2,psimin2);
% h=pcolor(r2,z2,Lc2d);
set(h,'linestyle','none');
plot(g.lim(1,:),g.lim(2,:),'k-')