clearvars; pc = phys_const;
EXP_DATA = 1;

PLOT_UPSTREAM_PSIN = 1;
PLOT_DOWNSTREAM_R = 1;  % Plot osp and isp

%%
temp = load('/home/jjl/180533-LMode/hf_180533_processed.mat');
expt.downstream.R = temp.Data.R;
expt.downstream.Z = temp.Data.Z;
expt.downstream.qperp_MWm2 = temp.Data.irDataLinearOffset;
expt.downstream.qperp_err_MWm2 = temp.Data.irErr;
clear temp;

temp =  load('/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/EXP_DATA/2300_upstream_fit_data_imp.mat');
expt.upstream = temp;
clear temp;
 
temp = load('/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/EXP_DATA/lp_180533_2250_3000.mat');
LP = temp.LP;
clear temp;

%%
myleg = '';

i = 0;
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2400/P5MW_n1p96e19_D+C_nosput';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/Proposed_shape_LSN_v2/P5MW_n2p5e19_xport_copy';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/Proposed_shape_LSN_v2/P5MW_n3e19_xport_copy';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/Proposed_shape_LSN_v2/P5MW_profile_fit';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_Csput';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_Csput_Ychem1pc';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_shift';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_shift_out';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_shift_Csput';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_shift_Csput_Y2pc';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_shift_out_Csput_Y2pc';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_shift_out_Csput_Y4pc';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_shift_out_Csput_Y4pc_allwall';
% % i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_Csput_Ychem2pc';
% % i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_Csput_Ychem3pc';
% % i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_Csput_Ychem6pc';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_Csput_Ychem10pc';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_Csput_Ychem5pc_coredensityBC';
% varAux = [0,0.5,1,2,3,6,10,11];
% myleg = {'D only','D+C phys only','D+C Y_{chem} 1%','D+C Y_{chem} 2%','D+C Y_{chem} 3%','D+C Y_{chem} 6%','D+C Y_{chem} 10%','EXPT'};

% myleg = {'case 1','case 2','case 3'};

% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_n2.5e19';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_n3.0e19';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_n3.5e19';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_n4e19';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_n5e19';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_n4e19_Ychem6pc';
% i = i+1; run_path{i} = '/data1/jjl/solps-iter/runs/d3d/NegD/180533/2300/test_Csput_Ychem6pc';

myleg = {'case 1','case 2','case 3','case 4','case 5','case 6'};
%% Plot settings
colors = lines;
styles = styflipper(length(run_path));
syms = symflipper(length(run_path));

fontSize = 14;
lineWidth = 2;

%%

% g = readg_g3d(find_solps_file(run_path{1},'gfile'));


%     osp.psiN = calc_psiN(g,osp.R,osp.Z);

%%
for i = 1:length(run_path)
    fprintf('Working on %s\n',run_path{i})
    
    g = readg_g3d(find_solps_file(run_path{i},'gfile'));
    
    Geo{i} = read_b2fgmtry(find_solps_file(run_path{i},'b2fgmtry'),0,1);
    Inds{i} = get_solps_jxa_jxi(run_path{i},Geo{i});
    DsData{i} = get_solps_ds_data(Geo{i},run_path{i});
    
    profs_jxa{i}.r1d_left = Geo{i}.r2d_left(Inds{i}.jxaMatlab,:).';
    profs_jxa{i}.z1d_left = Geo{i}.z2d_left(Inds{i}.jxaMatlab,:).';
    profs_jxa{i}.dsa = DsData{i}.dsa;
        
    profs_odiv{i}.r1d_left = Geo{i}.r2d_left(Inds{i}.jOuterFluxMatlab,:).';
    profs_odiv{i}.z1d_left = Geo{i}.z2d_left(Inds{i}.jOuterFluxMatlab,:).';
    profs_odiv{i}.L = [0;cumsum(sqrt(diff(profs_odiv{i}.r1d_left).^2 + diff(profs_odiv{i}.z1d_left).^2))];
    profs_odiv{i}.dsr = DsData{i}.dsr;
    
    profs_idiv{i}.r1d_left = Geo{i}.r2d_left(Inds{i}.jInnerFluxMatlab,:).';
    profs_idiv{i}.z1d_left = Geo{i}.z2d_left(Inds{i}.jInnerFluxMatlab,:).';
    profs_idiv{i}.dsl = DsData{i}.dsl;
       
    profs_jxa{i}.psiN  = calc_psiN(g,profs_jxa{i}.r1d_left,profs_jxa{i}.z1d_left);
    profs_odiv{i}.psiN = calc_psiN(g,profs_odiv{i}.r1d_left,profs_odiv{i}.z1d_left);
    profs_odiv{i}.Rsep = interp1(profs_odiv{i}.psiN,profs_odiv{i}.r1d_left,1);
    profs_odiv{i}.Zsep = interp1(profs_odiv{i}.psiN,profs_odiv{i}.z1d_left,1);
    profs_odiv{i}.Lsep = interp1(profs_odiv{i}.psiN,profs_odiv{i}.L,1);
    profs_odiv{i}.psiN = calc_psiN(g,profs_odiv{i}.r1d_left,profs_odiv{i}.z1d_left);
    
    State{i} = read_b2fstate(fullfile(run_path{i},'b2fstate'));
    
    State{i}.fht = calc_fht(State{i},Geo{i});
    Ion{i} = identify_solps_ion(State{i});
    
    profs_jxa{i}.te_keV = State{i}.te(Inds{i}.jxaMatlab,:)/pc.eV/1e3;
    profs_jxa{i}.ti_keV = State{i}.ti(Inds{i}.jxaMatlab,:)/pc.eV/1e3;
    profs_jxa{i}.ne_m3  = State{i}.ne(Inds{i}.jxaMatlab,:);
    profs_jxa{i}.ni_m3  = State{i}.na(Inds{i}.jxaMatlab,:,Ion{i}.isMainMatlab);
    
    profs_jxa{i}.zfzc6 = 6*State{i}.na(Inds{i}.jxaMatlab,:,end)./profs_jxa{i}.ne_m3;
    
    
    profs_odiv{i}.te_eV = State{i}.te(Inds{i}.jOuterVolumeMatlab,:)/pc.eV;
    profs_odiv{i}.ti_eV = State{i}.ti(Inds{i}.jOuterVolumeMatlab,:)/pc.eV;
    profs_odiv{i}.ne_m3 = State{i}.ne(Inds{i}.jOuterVolumeMatlab,:);
    profs_odiv{i}.ni_m3 = State{i}.na(Inds{i}.jOuterVolumeMatlab,:,Ion{i}.isMainMatlab);
    profs_odiv{i}.qperp_MWm2 = State{i}.fht(Inds{i}.jOuterFluxMatlab,:,1)./Geo{i}.sx(Inds{i}.jOuterFluxMatlab,:)./1e6;
    
    
    profs_idiv{i}.qperp_MWm2 = -State{i}.fht(Inds{i}.jInnerFluxMatlab,:,1)./Geo{i}.sx(Inds{i}.jInnerFluxMatlab,:)./1e6;
    profs_idiv{i}.te_eV = State{i}.te(Inds{i}.jInnerVolumeMatlab,:)/pc.eV;
    profs_idiv{i}.ne_m3 = State{i}.ne(Inds{i}.jInnerVolumeMatlab,:);
    

end


xlim = [0.9,1.05];


%%
if PLOT_UPSTREAM_PSIN
        
    ylim = [0,0.7];
    figure; hold on; box on; grid on;set(gcf,'color','w'); set(gca,'fontsize',fontSize);
    xlabel('\psi_N','fontsize',fontSize)
    ylabel('T_e (keV)','fontsize',fontSize)        
    errorbar(expt.upstream.teData.data.psiN,expt.upstream.teData.data.te./1e3,expt.upstream.teData.data.teErr./1e3,'k.')
    plot(expt.upstream.teData.fit.psiN,expt.upstream.teData.fit.te./1e3,'k','linewidth',2,'HandleVisibility','off')
    for i = 1:length(run_path)
        plot(profs_jxa{i}.psiN(2:end),profs_jxa{i}.te_keV(2:end),'linewidth',lineWidth,'color',colors(i,:))
    end    
    legend(['expt',myleg])
    axis([xlim(1),xlim(2),ylim(1),ylim(2)])
    
    ylim = [0,3];
    figure; hold on; box on; grid on;set(gcf,'color','w'); set(gca,'fontsize',fontSize);
    xlabel('\psi_N','fontsize',fontSize)
    ylabel('n_e (10^{19} m^{-3})','fontsize',fontSize)        
    errorbar(expt.upstream.neData.data.psiN,expt.upstream.neData.data.ne./1e19,expt.upstream.neData.data.neErr./1e19,'k.')
    plot(expt.upstream.neData.fit.psiN,expt.upstream.neData.fit.ne./1e19,'k','linewidth',2,'HandleVisibility','off')
    for i = 1:length(run_path)
        plot(profs_jxa{i}.psiN(2:end),profs_jxa{i}.ne_m3(2:end)/1e19,'linewidth',lineWidth,'color',colors(i,:))
    end    
    legend(['expt',myleg])
    axis([xlim(1),xlim(2),ylim(1),ylim(2)])       
    
    
    ylim = [0,0.15];
    figure; hold on; box on; grid on;set(gcf,'color','w'); set(gca,'fontsize',fontSize);
    xlabel('\psi_N','fontsize',fontSize)
    ylabel('Zf_z','fontsize',fontSize)        
    errorbar(expt.upstream.zfzData.data.psiN,expt.upstream.zfzData.data.zfz,expt.upstream.zfzData.data.zfzErr,'k.')  
    plot(expt.upstream.zfzData.fit.psiN,expt.upstream.zfzData.fit.zfz,'k','linewidth',2,'HandleVisibility','off')
    for i = 1:length(run_path)
        plot(profs_jxa{i}.psiN(2:end),profs_jxa{i}.zfzc6(2:end),'linewidth',lineWidth,'color',colors(i,:))
    end    
    legend(['expt',myleg])
    axis([xlim(1),xlim(2),ylim(1),ylim(2)])        
        
end



if PLOT_DOWNSTREAM_R
    Rsep = profs_odiv{1}.Rsep;
    
    xlimR = [1.4,1.8];
    ylim = [0,4];
    figure; hold on; box on; grid on;set(gcf,'color','w'); set(gca,'fontsize',fontSize);
    xlabel('R_{div}','fontsize',fontSize)
    ylabel('q_\perp (MW/m^2)','fontsize',fontSize)      
    errorbar(expt.downstream.R,expt.downstream.qperp_MWm2,expt.downstream.qperp_err_MWm2,'k.','linewidth',1)
    plot(LP.hf_r_rsep+Rsep,medflt1d_jl(LP.hf,5),'rx','HandleVisibility','off')
    for i = 1:length(run_path)
        plot([profs_idiv{i}.r1d_left(2:end);NaN;profs_odiv{i}.r1d_left(2:end)],[profs_idiv{i}.qperp_MWm2(2:end),NaN,profs_odiv{i}.qperp_MWm2(2:end)],'linewidth',lineWidth,'color',colors(i,:),'linewidth',2)
    end       
    legend(['expt',myleg],'location','northwest')
    axis([xlimR(1),xlimR(2),ylim(1),ylim(2)])        
        
    
%     Rtest = linspace(1.372,1.768,100);
%     Ztest = -1.25*ones(size(Rtest));
%     alpha_deg = calc_Bangle_g(g,Rtest,Ztest);
%     cs = sqrt(LP.te*2*1.6e-19/2/1.67e-27);
%     qq = 8*LP.te*1.6e-19.*LP.dens*1e6.*cs;
%     aTest = interp1(Rtest,alpha_deg,LP.dens_r_rsep+Rsep);
%     plot(LP.te_r_rsep+Rsep,medflt1d_jl(qq,5)./1e6.*sind(abs(aTest)),'b.')
    
    ylim = [0,100];
    figure; hold on; box on; grid on;set(gcf,'color','w'); set(gca,'fontsize',fontSize);
    xlabel('R_{div} (m)','fontsize',fontSize)
    ylabel('T_e (eV)','fontsize',fontSize)      
    plot(LP.te_r_rsep+Rsep,medflt1d_jl(LP.te,5),'k.')
    for i = 1:length(run_path)
        plot([profs_idiv{i}.r1d_left(2:end);NaN;profs_odiv{i}.r1d_left(2:end)],[profs_idiv{i}.te_eV(2:end),NaN,profs_odiv{i}.te_eV(2:end)],'linewidth',lineWidth,'color',colors(i,:),'linewidth',2)
    end       
    legend(['expt',myleg],'location','northwest')
    axis([xlimR(1),xlimR(2),ylim(1),ylim(2)])   
    
    

    ylim = [0,12];
    figure; hold on; box on; grid on;set(gcf,'color','w'); set(gca,'fontsize',fontSize);
    xlabel('R_{div} (m)','fontsize',fontSize)
    ylabel('n_e (10^{19} m^{-3})','fontsize',fontSize)      

    plot(LP.dens_r_rsep+Rsep,medflt1d_jl(LP.dens,5)./1e13,'k.')
    for i = 1:length(run_path)
        plot([profs_idiv{i}.r1d_left(2:end);NaN;profs_odiv{i}.r1d_left(2:end)],1e-19.*[profs_idiv{i}.ne_m3(2:end),NaN,profs_odiv{i}.ne_m3(2:end)],'linewidth',lineWidth,'color',colors(i,:),'linewidth',2)
    end       
    legend(['expt',myleg],'location','northwest')
    axis([xlimR(1),xlimR(2),ylim(1),ylim(2)])  
    

end

