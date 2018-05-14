function plot_profdb_data(fname,psin_range,te_max,ti_max,ne_max,newfig,BIN_TS)
%  te, ti max in kev
%  ne max in 10^20
%
if nargin < 1
    error('fname required')
end
if nargin < 2
    psin_range = [0.90,1.1];
end
if nargin < 3
    te_max = [];
end
if nargin < 4
    ti_max = [];
end
if nargin < 5
    ne_max = [];
end
if nargin < 6
    newfig = 1;
end
if nargin < 7
    BIN_TS = 0;
end


shift = 0.0;

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


%--------------------------------------------------------------------------

% load([data_path,fname])
load(fname)

if BIN_TS
    nbin = round((binright-binleft)/dbin);
    psiNbin = binleft+dbin*[0:nbin-1]+0.5*dbin;
    for i=1:nbin
        ibin = find( profs.psi_te + shift >= binleft+(i-1)*dbin & profs.psi_te + shift < binleft+i*dbin);
        tebin(i) = mean(profs.tedat(ibin));
        teberr(i) = std(profs.tedat(ibin));
    end
    for i=1:nbin
        ibin = find( profs.psi_ne + shift >= binleft+(i-1)*dbin & profs.psi_ne + shift < binleft+i*dbin);
        nebin(i) = mean(profs.nedat(ibin));
        neberr(i) = std(profs.nedat(ibin));
    end        
end

tefit = evaluate_tanh_fit(profs.tetanh,psiNfit-shift);
if profs.titanhok
    tifit = evaluate_tanh_fit(profs.titanh,psiNfit-shift);
end
nefit = evaluate_tanh_fit(profs.ntanh,psiNfit-shift);        

if SCALE_TI == 1
    iscale = psiNfit >= 1.;
    psiNfit_sep = psiNfit(iscale); psiNfit_sep = psiNfit_sep(1);
    ti_sep = evaluate_tanh_fit(profs.titanh,1-shift);
    tifit(iscale) = (ti_sep-ti_exp_scale_min/1000)*exp(-(psiNfit(iscale)-psiNfit_sep)/ti_exp_scale_width) + ti_exp_scale_min/1000;
end

if newfig == 1
    figure; hold on; box on;
end
n1_sub = 3; n2_sub = 1;
subplot(n1_sub,n2_sub,1); hold on; box on; set(gca,'xlim',xrange); 
if ~isempty(te_max) 
    set(gca,'ylim',[0,te_max]);
end
if BIN_TS
    plot(psiNbin+shift,tebin)
    errorbar(psiNbin+shift,tebin,teberr,col,'marker','none','linestyle','none')
else
    plot(profs.psi_te+shift,profs.tedat)
    errorbar(profs.psi_te+shift,profs.tedat,profs.te_err,col,'marker','none','linestyle','none')
end
plot(psiNfit+shift,tefit,[col,'-'],'linewidth',3)
ylabel('T_e [keV]','fontsize',12)
title('Profiles vs \psi_N','fontsize',12)
set(gca,'fontsize',12)

subplot(n1_sub,n2_sub,2); hold on; box on; set(gca,'xlim',xrange); 
if ~isempty(ne_max)
    set(gca,'ylim',[0,ne_max]);
end
if BIN_TS
    errorbar(psiNbin+shift,nebin,neberr,col,'marker','none','linestyle','none')
else
    errorbar(profs.psi_ne+shift,profs.nedat,profs.ne_err,col,'marker','none','linestyle','none')
end
plot(psiNfit+shift,nefit,[col,'-'],'linewidth',3)
ylabel('n_e [10^{20} m^{-3}]','fontsize',12)
set(gca,'fontsize',12)

subplot(n1_sub,n2_sub,3); hold on; box on; set(gca,'xlim',xrange); 
if ~isempty(ti_max)
    set(gca,'ylim',[0,ti_max]);
end
errorbar(profs.psi_ti+shift,profs.tidat,profs.ti_err,col,'marker','none','linestyle','none')
if profs.titanhok
    plot(psiNfit+shift,tifit,[col,'-'],'linewidth',3)
end
xlabel('\psi_N','fontsize',12)
ylabel('T_i [keV]','fontsize',12)
set(gca,'fontsize',12)





