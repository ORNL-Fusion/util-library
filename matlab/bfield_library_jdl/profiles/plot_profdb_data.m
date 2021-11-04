function plot_profdb_data(fname,xType,g)
%  te, ti max in kev ,psin_range,te_max,ti_max,ne_max,newfig,BIN_TS
%  ne max in 10^20
%
if nargin < 1
    error('fname required')
end
if nargin < 2
    xType = 1;
end
if xType == 2
    if nargin < 3
        error('For xType 2 (R-Rsep), g must be provided')
    end
else
    xType = 1;
end


% if nargin < 2
    psin_range = [0.90,1.1];
% end
% if nargin < 3
    te_max = [];
% end
% if nargin < 4
    ti_max = [];
% end
% if nargin < 5
    ne_max = [];
% end
% if nargin < 6
    newfig = 1;
% end
% if nargin < 7
    BIN_TS = 0;
% end


psiNshift = 0.0;

col = 'k';

SCALE_TI = 0;
ti_exp_scale_min = 5; % eV
ti_exp_scale_width = 0.03; % delta psiN

% BIN_TS = 1;
% ;; Settings for binning TS data
binleft=0.85;    %; Range in psi_N
binright=1.1;
dbin = 0.005;

psiNfit = linspace(0.85,1.1,100);
xrange = psin_range;


if xType == 2
    setx = @(psiN) calc_rho_midplane_from_psiN(g,psiN);
else
    setx = @(psiN) psiN;
end
    

%%
% load([data_path,fname])
load(fname)

if BIN_TS
    nbin = round((binright-binleft)/dbin);
    psiNbin = binleft+dbin*[0:nbin-1]+0.5*dbin;
    for i=1:nbin
        ibin = find( profs.psi_te + psiNshift >= binleft+(i-1)*dbin & profs.psi_te + psiNshift < binleft+i*dbin);
        tebin(i) = mean(profs.tedat(ibin));
        teberr(i) = std(profs.tedat(ibin));
    end
    for i=1:nbin
        ibin = find( profs.psi_ne + psiNshift >= binleft+(i-1)*dbin & profs.psi_ne + psiNshift < binleft+i*dbin);
        nebin(i) = mean(profs.nedat(ibin));
        neberr(i) = std(profs.nedat(ibin));
    end        
end

tefit = evaluate_tanh_fit(profs.tetanh,psiNfit-psiNshift);
if profs.titanhok
    tifit = evaluate_tanh_fit(profs.titanh,psiNfit-psiNshift);
end
nefit = evaluate_tanh_fit(profs.ntanh,psiNfit-psiNshift);        

if SCALE_TI == 1
    iscale = psiNfit >= 1.;
    psiNfit_sep = psiNfit(iscale); psiNfit_sep = psiNfit_sep(1);
    ti_sep = evaluate_tanh_fit(profs.titanh,1-psiNshift);
    tifit(iscale) = (ti_sep-ti_exp_scale_min/1000)*exp(-(psiNfit(iscale)-psiNfit_sep)/ti_exp_scale_width) + ti_exp_scale_min/1000;
end

%% Plot
if newfig == 1
    figure; hold on; box on; set(gcf,'color','w');
end
n1_sub = 3; n2_sub = 1;

subplot(n1_sub,n2_sub,1); hold on; box on; grid on; set(gca,'fontsize',14);
set(gca,'xlim',setx(xrange));
if ~isempty(te_max) 
    set(gca,'ylim',[0,te_max]);
end
if BIN_TS
    errorbar(setx(psiNbin+psiNshift),tebin,teberr,col,'marker','none','linestyle','none')
else
    errorbar(setx(profs.psi_te+psiNshift),profs.tedat,profs.te_err,col,'marker','none','linestyle','none')
end
plot(setx(psiNfit+psiNshift),tefit,[col,'-'],'linewidth',3)
ylabel('T_e [keV]')
if xType == 2
    title('Profiles vs dR_{sep}^{OMP}')    
else
    title('Profiles vs \psi_N')
end
subplot(n1_sub,n2_sub,2); hold on; box on; grid on; set(gca,'fontsize',14);
set(gca,'xlim',setx(xrange));
if ~isempty(ne_max)
    set(gca,'ylim',[0,ne_max]);
end
if BIN_TS
    errorbar(setx(psiNbin+psiNshift),nebin,neberr,col,'marker','none','linestyle','none')
else
    errorbar(setx(profs.psi_ne+psiNshift),profs.nedat,profs.ne_err,col,'marker','none','linestyle','none')
end
plot(setx(psiNfit+psiNshift),nefit,[col,'-'],'linewidth',3)
ylabel('n_e [10^{20} m^{-3}]')

subplot(n1_sub,n2_sub,3); hold on; box on;  grid on; set(gca,'fontsize',14);
set(gca,'xlim',setx(xrange));
if ~isempty(ti_max)
    set(gca,'ylim',[0,ti_max]);
end
errorbar(setx(profs.psi_ti+psiNshift),profs.tidat,profs.ti_err,col,'marker','none','linestyle','none')
if profs.titanhok
    plot(setx(psiNfit+psiNshift),tifit,[col,'-'],'linewidth',3)
end
if xType == 2
    xlabel('dR_{sep}^{OMP}')
else
    xlabel('\psi_N')
end
ylabel('T_i [keV]')






