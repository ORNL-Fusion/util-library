clearvars;

% config_name = {'D1-1','D3-6','D4-1','D4-2','D4-3','D4-4'};
% config_name = {'D3-6','D1-1','D1-2','D1-3','D1-4','D1-5','D1-6'};
config_name = {'D3-6','D2-2','D1-1'};
% config_name = {'D3-6'};
verbose = 1;
PLOT_COILS = 0;
field_models = {'jackson','hybrid'};
primary_model = 'hybrid';
hybridSimplifyCoils = [];
hybridNTurns = 3;
hybridNLayers = 3;

nTest = 200;
zTest = linspace(-2,10,nTest);
rTest = 0*ones(size(zTest));

% geo = get_MPEX_geometry;

%%
figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14)
ax = gca;
colors = lines(length(config_name));
for i = 1:length(config_name)
    model_data = struct();

    for imodel = 1:length(field_models)
        model_name = field_models{imodel};
        switch model_name
            case 'jackson'
                [coil,current,fil,cur] = build_MPEX_coils_jackson(config_name{i},verbose);
            case 'hybrid'
                [coil,current,fil,cur] = build_MPEX_coils_jackson_hybrid( ...
                    config_name{i},hybridSimplifyCoils,hybridNTurns,hybridNLayers,verbose);
            otherwise
                error('Unsupported field model: %s',model_name)
        end

        bfield.coil = coil;
        bfield.current = current;
        bfield.type = 'MPEX';

        [Br,Bz] = bfield_circular_coils_vectorized(bfield.coil,bfield.current,rTest,zTest);
        model_data.(model_name).Btot = sqrt(Br.^2 + Bz.^2);
        model_data.(model_name).fil = fil;
        model_data.(model_name).cur = cur;
    end

    if PLOT_COILS && i == 1
        [rcoil,zcoil] = get_coil_cross_sections(model_data.(primary_model).fil);
        newfig = 1;
        plot_coil_cross_section(rcoil,zcoil,newfig,model_data.(primary_model).cur);
    end

    for imodel = 1:length(field_models)
        model_name = field_models{imodel};
        if strcmp(model_name,primary_model)
            line_style = '-';
        else
            line_style = '--';
        end
        plot(ax,zTest,model_data.(model_name).Btot, ...
            'Color',colors(i,:), ...
            'LineStyle',line_style, ...
            'LineWidth',2, ...
            'DisplayName',sprintf('%s (%s)',config_name{i},model_name))
    end
end

legend('Location','best')
xlabel('Z (m)')
ylabel('|B| on-axis (T)')

% plot(geo.z,geo.r,'r-','LineWidth',2)
