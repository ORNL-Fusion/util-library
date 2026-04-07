clearvars;

% config_name = {'D1-1','D3-6','D4-1','D4-2','D4-3','D4-4'};
% config_name = {'D3-6','D1-1','D1-2','D1-3','D1-4','D1-5','D1-6'};
% config_name = {'D3-6','D1-1','D1-2','D2-2'};
config_name = {'D3-6'};
verbose = 1;
PLOT_COILS = 1;

nTest = 200;
zTest = linspace(-2,10,nTest);
rTest = 0*ones(size(zTest));

% geo = get_MPEX_geometry;

%%
figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14)
for i = 1:length(config_name)
    % [fil,cur] = setup_MPEX_coils(config_name{i},0);
    [coil,current,fil,cur] = build_MPEX_coils_jackson(config_name{i},verbose);
    bfield.coil = coil;
    bfield.current = current;
    bfield.type = 'MPEX';

    if PLOT_COILS && i == 1
        [rcoil,zcoil] = get_coil_cross_sections(fil);
        newfig = 1;
        plot_coil_cross_section(rcoil,zcoil,newfig,cur);
    end

    [Br,Bz] = bfield_circular_coils_vectorized(bfield.coil,bfield.current,rTest,zTest);
    Btot = sqrt(Br.^2 + Bz.^2);
    plot(zTest,Btot,'DisplayName',config_name{i},'LineWidth',2)
end

legend(config_name,'Location','best')
xlabel('Z (m)')
xlabel('R (m)')
ylabel('|B| on-axis (T)')

% plot(geo.z,geo.r,'r-','LineWidth',2)