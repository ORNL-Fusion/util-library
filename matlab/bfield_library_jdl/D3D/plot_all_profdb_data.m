clearvars;

% files{1} = 'C:\Work\DIII-D\APS 2016\profdb_data\profs_156855_4500_awlhm.mat'; mytitle{1} = '156855';
% files{2} = 'C:\Work\DIII-D\APS 2016\profdb_data\profs_156859_4500_awlhm.mat'; mytitle{2} = '156859';
% % files{2} = 'C:\Work\DIII-D\APS 2016\profdb_data\profs_156861_4500_awlhm.mat'; mytitle{2} = '156861';
% files{3} = 'C:\Work\DIII-D\APS 2016\profdb_data\profs_156867_4500_awlhm.mat'; mytitle{3} = '156867';
% files{4} = 'C:\Work\DIII-D\APS 2016\profdb_data\profs_156871_4500_awlhm.mat'; mytitle{4} = '156871';


files{1} = 'C:\Work\DIII-D\APS 2016\profdb_data\profs_156855_4500_awlhm.mat'; mytitle{1} = '156855';
% files{2} = 'C:\Work\DIII-D\APS 2016\profdb_data\profs_156855_3500_awlhm.mat'; mytitle{2} = '156855';


for i = 1:length(files)
    myprofs{i} = load(files{i});
end


BIN_TS = 1;
myshifts = [0,0,0,0];
psiNfit = linspace(0.85,1.1,100);
psin_range = [0.85,1.1];
xrange = psin_range;
newfig = 1;
te_max = 0.7;
ne_max = 1;
ti_max = 0.8;
col_data = 'k';
col_spl = 'r';
col_tanh = 'b';

binleft=0.85;    %; Range in psi_N
binright=1.1;
dbin = 0.005;

for i = 1:length(files)
    clear myleg
    profs = myprofs{i}.profs;           
    shift = myshifts(i);
    if BIN_TS
        nbin = round((binright-binleft)/dbin);
        psiNbin = binleft+dbin*[0:nbin-1]+0.5*dbin;
        for j=1:nbin
            ibin = find( profs.psi_te + shift >= binleft+(j-1)*dbin & profs.psi_te + shift < binleft+j*dbin);
            tebin(j) = mean(profs.tedat(ibin));
            teberr(j) = std(profs.tedat(ibin));
        end
        for j=1:nbin
            ibin = find( profs.psi_ne + shift >= binleft+(j-1)*dbin & profs.psi_ne + shift < binleft+j*dbin);
            nebin(j) = mean(profs.nedat(ibin));
            neberr(j) = std(profs.nedat(ibin));
        end
    end    
    tefit = evaluate_tanh_fit(profs.tetanh,psiNfit-shift);
    if profs.titanhok
        tifit = evaluate_tanh_fit(profs.titanh,psiNfit-shift);
    end
    tifit_spl = interp1(profs.psi_tispl,profs.tispl,psiNfit-shift);
    nefit = evaluate_tanh_fit(profs.ntanh,psiNfit-shift);
    
    %
    % TE
    %  
    if newfig == 1
        figure; hold on; box on;
    end    
    n1_sub = 3; n2_sub = 1;
    subplot(n1_sub,n2_sub,1); hold on; box on; set(gca,'xlim',xrange);
    title(mytitle{i});
    if ~isempty(te_max)
        set(gca,'ylim',[0,te_max]);
    end
    if BIN_TS
        errorbar(psiNbin+shift,tebin,teberr,col_data,'marker','o','linestyle','none')
    else
        errorbar(profs.psi_te+shift,profs.tedat,profs.te_err,col_data,'marker','none','linestyle','none')
    end
    plot(psiNfit+shift,tefit,[col_tanh,'-'],'linewidth',3)        
    ylabel('T_e [keV]','fontsize',12)
    set(gca,'fontsize',12)
    
    
    %
    % NE
    %
    subplot(n1_sub,n2_sub,2); hold on; box on; set(gca,'xlim',xrange);
    if ~isempty(ne_max)
        set(gca,'ylim',[0,ne_max]);
    end
    if BIN_TS
        errorbar(psiNbin+shift,nebin,neberr,col_data,'marker','o','linestyle','none')
    else
        errorbar(profs.psi_ne+shift,profs.nedat,profs.ne_err,col_data,'marker','none','linestyle','none')
    end
    plot(psiNfit+shift,nefit,[col_tanh,'-'],'linewidth',3)
    ylabel('n_e [10^{20} m^{-3}]','fontsize',12)
    set(gca,'fontsize',12)

    %
    % TI
    %
    subplot(n1_sub,n2_sub,3); hold on; box on; set(gca,'xlim',xrange);
    if ~isempty(ti_max)
        set(gca,'ylim',[0,ti_max]);
    end
    myleg{1} = 'data';
    errorbar(profs.psi_ti+shift,profs.tidat,profs.ti_err,col_data,'marker','none','linestyle','none')
    if profs.titanhok
        plot(psiNfit+shift,tifit,[col_tanh,'-'],'linewidth',3)
        myleg = [myleg,{'tanh'}];
    end
    plot(psiNfit+shift,tifit_spl,[col_spl,'-'],'linewidth',3)
    myleg = [myleg,{'spline'}];
    legend(myleg)
    
    xlabel('\Psi_N','fontsize',12)
    ylabel('T_i [keV]','fontsize',12)
    set(gca,'fontsize',12)
end

