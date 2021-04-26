clearvars;

gfile_name = 'C:\Work\DIII-D\APS 2016\g156855.04526_694';
g = readg_g3d(gfile_name);

files{1} = 'C:\Work\DIII-D\APS 2016\profs_156855_4500_awlhm.mat'; mytitle{1} = '156855';
files{2} = 'C:\Work\DIII-D\APS 2016\profs_156859_4500_awlhm.mat'; mytitle{2} = '156859';
% files{2} = 'C:\Work\DIII-D\APS 2016\profs_156861_4500_awlhm.mat'; mytitle{2} = '156861';
files{3} = 'C:\Work\DIII-D\APS 2016\profs_156867_4500_awlhm.mat'; mytitle{3} = '156867';
files{4} = 'C:\Work\DIII-D\APS 2016\profs_156871_4500_awlhm.mat'; mytitle{4} = '156871';


for i = 1:length(files)
    myprofs{i} = load(files{i});
end

psiNfit = linspace(0.95,1.05,1000);

psiN_core_SOLPS = 0.9;

tewants = [63,63,85,83];
shifts = zeros(size(tewants));

for i = 1:length(files)
    tewant = tewants(i)/1000;
    profs = myprofs{i}.profs;
    tefit = evaluate_tanh_fit(profs.tetanh,psiNfit);    
    Te_sep_noshift = evaluate_tanh_fit(profs.tetanh,1);
    shifts(i) = 1-interp1(tefit,psiNfit,tewant);
    Te_check = evaluate_tanh_fit(profs.tetanh,1-shifts(i));
    [rs,zs] = calc_RZ_at_psiN_theta(g,[1,1-shifts(i)],0);
    real_shift = -diff(rs);
    
    % Find n at core SOLPS boundary
    nfit = evaluate_tanh_fit(profs.ntanh,psiN_core_SOLPS-shifts(i));
    fprintf('\n')
    fprintf('%s\n',mytitle{i})
    fprintf('Te_sep, no shift = %f [eV]\n',Te_sep_noshift*1000)
    fprintf('Want Te_sep      = %f [eV]\n',tewant*1000)
    fprintf('Got Te_sep       = %f [eV]\n',Te_check*1000)
    fprintf('psiN shift       = %f [eV]\n',shifts(i)*1000)
    fprintf('dR shift [mm]    = %f [eV]\n',real_shift*1000)
    fprintf('ne at core bdry  = %f [e20] \n',nfit)
%     fprintf('sqrt(pshift)[~mm]= %f [eV]\n',sqrt(shifts(i)*1000))
end