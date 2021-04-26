function plot_2DTS_fit_data(rcens,zcens,lcens)
% clearvars;


switch 1
    case 1
        fname_out_part = '2MW_attach';
        fpath = 'C:\Work\DIII-D\APS 2016 - H mode 2 and 5MW\2dts_data\';
        fname ='dts_fitted1d_osp_multishot_156855156856156857156858_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';
        fname_fits = 'dts_2d_dataandfit_multishot_156855156856156857156858_remapto156855at4527msEFIT04_tsrev-1_ELMfilter_chisq10.mat';
    case 2
        fname_out_part = '2MW_detach';
        fname ='C:\Work\DIII-D\APS 2016 - H mode 2 and 5MW\2dts_data\dts_fitted1d_osp_multishot_156859156861156862156863_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';
    case 3
        fname_out_part = '5MW_attach';
        fname ='C:\Work\DIII-D\APS 2016 - H mode 2 and 5MW\2dts_data\dts_fitted1d_osp_multishot_156866156867156868156869_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';
    case 4
        fname_out_part = '5MW_detach';
        fname ='C:\Work\DIII-D\APS 2016 - H mode 2 and 5MW\2dts_data\dts_fitted1d_osp_multishot_156871156872156873156874_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';        
end

ts2d = read_2dts_dat_file(fullfile(fpath,fname));
f = load(fullfile(fpath,fname_fits));


f = f.dtsfit_all;
r1 = f.fit.rmap(:,1);
z1 = f.fit.zmap(1,:);
[r2,z2] = meshgrid(r1,z1);
d2 = f.fit.dens;
t2 = f.fit.temp;
% [d2d] = griddata(r1,z1,d2,r2,z2);

dcens = interp2(r2,z2,d2,rcens,zcens,'linear',NaN);
tcens = interp2(r2,z2,t2,rcens,zcens,'linear',NaN);
a=1