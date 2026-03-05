clearvars;
pc = phys_const;

suffix = '89'; element='w'; mass_amu = 183.84;
% suffix = '96'; element='c';
% suffix = '50'; element='w';
% suffix = '96'; element='ne';

ADAS_dir = 'C:\Users\jjl\ORNL Dropbox\Jeremy Lore\ADAS\adf11_all';

fprintf('Reading scd\n');
fname = fullfile(ADAS_dir,strcat('scd',suffix),strcat('scd',suffix,'_',element,'.dat'));  % Effective ionization coefficients (cm^-3/s)
scd = read_adas_adf11_file(fname);



B = 4;

indReac = 1;
RateCoeff = scd;
te = linspace(3,20,100);
ne = 1e18;

vp = sqrt(2*300*pc.kB/(pc.amu0*mass_amu));
gyro = pc.amu0*mass_amu*vp/(pc.e*B)*1000
[mfp1] = calc_mfp(te,ne,RateCoeff,indReac,vp);

vp = sqrt(2*1*pc.e/(pc.amu0*mass_amu));
gyro = pc.amu0*mass_amu*vp/(pc.e*B)*1000
[mfp2] = calc_mfp(te,ne,RateCoeff,indReac,vp);

vp = sqrt(2*3*pc.e/(pc.amu0*mass_amu));
gyro = pc.amu0*mass_amu*vp/(pc.e*B)*1000
[mfp3] = calc_mfp(te,ne,RateCoeff,indReac,vp);

figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14)
plot(te,mfp1*100,'LineWidth',2)
plot(te,mfp2*100,'LineWidth',2)
plot(te,mfp3*100,'LineWidth',2)
set(gca,'yscale','log')
xlabel('T_e (eV)')
ylabel('\lambda_{mfp} (cm)')
legend('W_0 = 300 K','W_0 = 1 eV','W_0 = 3 eV')
set(gca,'YTick',[1,2,5,10,20,50,100])
