clearvars;

element = 'N'; % 'Ar', 'C', 'Kr', 'N', 'Ne'

Data = read_Lz(element);
Te_eval = logspace(log10(.1),log10(1000),1000);  % in eV
tau_eval = [0.1,1,10,1000];                  % in ms


for i = 1:length(tau_eval)
    Lz_eval(:,i) = exp(interp2(log(Data.Te_eV),log(Data.tau_ms),log(Data.Lz_Wm3),log(Te_eval),log(tau_eval(i))));
    myLeg{i} = strcat('n_e\tau',sprintf(' = %.1f [%.0e m^{-3}ms]',tau_eval(i),Data.neRef));
end

figure; hold on; box on; grid on; set(gcf,'color','w');set(gca,'fontsize',14);
plot(Te_eval,Lz_eval,'LineWidth',2)
set(gca,'xscale','log','yscale','log')
xlabel('T_e (eV)')
ylabel(strcat('L_z^{',element,'} (Wm^3)'))
legend(myLeg,'fontsize',10)
title(element)


TeNormVal_eV = 15;
TeNorm = Te_eval./TeNormVal_eV;
art_rad =1;

rqa = art_rad*2e-31./(TeNorm.^1.5+1./TeNorm.^3);  % Heat loss rate coefficient
rrd = art_rad*2e-31./(TeNorm.^1.5+1./TeNorm.^3);  % Line radiation rate coefficient
plot(Te_eval,rqa)
function Data = read_Lz(element)

%% Read file
switch element
    case 'Ar'
        fName = 'Ar_lz_tau.dat';
    case 'C'
        fName = 'C__lz_tau.dat';
    case 'Kr'
        fName = 'Kr_lz_tau.dat';
    case 'N'
        fName = 'N__lz_tau.dat';
    case 'Ne'
        fName = 'Ne_lz_tau.dat';
    otherwise
        error('No file for element %s',element)
end

fprintf('Reading file: %s\n',fName);

fid = fopen(fName,'r');
Data.header{1} = fgetl(fid);
fprintf('%s\n',Data.header{1});
Data.header{2} = fgetl(fid);
temp = fgetl(fid);
i1 = strfind(temp,'=');
i2 = strfind(temp,'m');
Data.neRef = sscanf(temp(i1+1:i2-1),'%e');
temp = fgetl(fid);
temp = sscanf(temp,'%d',2);
Data.nTe = temp(1);
Data.nTau = temp(2);
fprintf('nTe = %d, nTau = %d\n',Data.nTe,Data.nTau);
fgetl(fid);
Data.Te_eV = fscanf(fid,'%e\n',Data.nTe);
Data.Lz_Wm3 = nan(Data.nTau,Data.nTe);
Data.tau_ms = nan(Data.nTau,1);
for i = 1:Data.nTau
    temp = fgetl(fid);
    i1 = strfind(temp,'=');
    Data.tau_ms(i) = sscanf(temp(i1+1:end),'%f');
    Data.Lz_Wm3(i,:) = fscanf(fid,'%e\n',Data.nTe);
end

fclose(fid);

end