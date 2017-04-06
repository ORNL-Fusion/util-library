    clearvars;
shots = 7400 + [0,3:6,8,10,12:13,16,17,18]; mytitle = 'I_A = 6368 A, no skimmer';  x0_guess = -.5512;y0_guess = -2.533; force_guess = 1;
% shots = 7400 + [3,8,12,18]; mytitle = 'I_A = 3300 A, no skimmer';
% shots = 7400 + [77,87,88,92:98,100,101,103]; mytitle = 'I_A = 3300 A, with skimmer';x0_guess = -.5877;y0_guess = -2.8914; force_guess = 1;
ns = length(shots);
% x0_guess = -.4;y0_guess = -3.2;



scale_H = 1;

plotitIR = 1;
for i = 1:ns
    shot = shots(i);
    [helicon_current(i),current_A,current_B,config,skimmer] = get_Proto_current(shot);
    [rr{i},dd{i},radius(i),angle(i),x0_final(i),y0_final(i)] = plot_IR_data_raw(shot,plotitIR,x0_guess,y0_guess,force_guess);
    [maxIR(i),mind] = max(dd{i});
    radius_maxIR(i) = rr{i}(mind);
end
asdgfaad
if scale_H == 1
    helicon_current = helicon_current*(3300/current_A);
    htitle = 'I_H(3300/I_A) [A]';
else
    htitle = 'I_H [A]';
end
figure; hold on; box on;
plot(shots,helicon_current,'bo')
xlabel('Shot','fontsize',14)
ylabel(htitle,'fontsize',14)
set(gca,'fontsize',14)
title(mytitle,'fontsize',14)


cmax = 600;
figure; hold on; box on;
cf = colorflipper(256,'parula');
cvals = linspace(0,cmax,size(cf,1));
for i = 1:ns
    cind = round(interp1(cvals,1:length(cvals),helicon_current(i)));
    plot(rr{i},dd{i},'color',cf(cind,:))    
end
colorbar;
set(gca,'clim',[0,cmax])
colormap(cf);
xlabel('r [cm]','fontsize',14)
ylabel('\DeltaT [C]','fontsize',14)
title(mytitle,'fontsize',14)

figure; hold on; box on;
plot(helicon_current,radius/100,'k.')
plot(helicon_current,radius_maxIR/100,'bo')
ylabel('radius [m]','fontsize',14)
xlabel(htitle,'fontsize',14)
set(gca,'fontsize',14)
title(mytitle,'fontsize',14)


figure; hold on; box on;
plot(helicon_current,angle*180/pi,'bo')
ylabel('angle [deg]','fontsize',14)
xlabel(htitle,'fontsize',14)
set(gca,'fontsize',14)
title(mytitle,'fontsize',14)

figure; hold on; box on;
plot(helicon_current,maxIR,'bo')
ylabel('max \DeltaT [deg C]','fontsize',14)
xlabel(htitle,'fontsize',14)
set(gca,'fontsize',14)
title(mytitle,'fontsize',14)

figure; hold on; box on;
plot(x0_final,y0_final,'bo')
set(gca,'fontsize',14)
title(mytitle,'fontsize',14)