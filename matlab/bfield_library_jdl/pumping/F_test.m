clearvars;

Hpump = 0.03;

Rent = 1.6;

Rprof = linspace(1.5,1.59,100);
dR = Rprof(2)-Rprof(1);

F = (1 - 1./sqrt(1 + (Hpump./(Rent - Rprof)).^2))/2;
% th = atan(Hpump/(Rprof-Rent));
% F2 = (1-cos(th))/2