clearvars;


config_name = {'D3-6'};
verbose = 1;



PLOT_COILS = 1;


geo = get_MPEX_geometry;

%%
for i = 1:length(config_name)
    figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14)
    [coil,current,fil,cur] = build_MPEX_coils_jackson(config_name{i},verbose);
    bfield.coil = coil;
    bfield.current = current;
    bfield.type = 'MPEX';

    if PLOT_COILS
        [rcoil,zcoil] = get_coil_cross_sections(fil);
        plot_coil_cross_section(rcoil,zcoil,0,cur);
    end

    
end

xlabel('Z (m)')
xlabel('R (m)')
plot(geo.vessel.z,geo.vessel.r,'r-','LineWidth',2)
axis([-4,9,0,0.5])

