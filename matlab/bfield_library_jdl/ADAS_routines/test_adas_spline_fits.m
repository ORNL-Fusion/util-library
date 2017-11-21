clearvars;

species_string = 'h';
[scd,acd,plt,prb,ccd] = setup_rates(species_string);

ne = 1e19;
ntry = 10000;
Te = linspace(1,1000,ntry);

% Te = 1;

tic1 = tic;
for i = 1:ntry
    sv_iz(i) = calc_sv_ee_internal(Te(i),ne,scd,1);
end
toc(tic1)

% sv_iz = calc_sv_ee_internal(Te,ne,scd,0)

% [xx,yy,zz] = prepareSurfaceData(scd.te_log10,scd.ne_log10,scd.coeff_log10);
% [fit_data,gof] = fit( [xx,yy],zz, 'cubicinterp' );
% tic1 = tic;
% for i = 1:ntry
%     1e-6*10^feval(fit_data,log10(Te(i)),log10(ne/1e6));
% end
% toc(tic1)


[tt,nn] = ndgrid(scd.te_log10,scd.ne_log10);
F = griddedInterpolant(tt,nn,scd.coeff_log10.','spline');
tic1 = tic;
for i = 1:ntry
    sv_iz_gd(i) = 1e-6*10^F(log10(Te(i)),log10(ne/1e6));
end
toc(tic1)

figure; hold on; box on;
plot(Te,sv_iz)
plot(Te,sv_iz_gd)

function [scd,acd,plt,prb,ccd] = setup_rates(species_string)
    adas_path = 'C:/Work/ADAS/adf11_all';
    %     species_string = 'h';
    
    fname = fullfile(adas_path,strcat('scd96/scd96_',species_string,'.dat'));  % Effective ionization coefficients (cm^-3/s)
    scd = read_adas_adf11_file(fname);
    fname = fullfile(adas_path,strcat('acd96/acd96_',species_string,'.dat'));  % Effective recombination coefficients (cm^-3/s)
    acd = read_adas_adf11_file(fname);
    fname = fullfile(adas_path,strcat('plt96/plt96_',species_string,'.dat'));  % Radiated power IZ (W/cm^3)
    plt = read_adas_adf11_file(fname);
    fname = fullfile(adas_path,strcat('prb96/prb96_',species_string,'.dat'));  % Radiated power RC (W/cm^3)
    prb = read_adas_adf11_file(fname);
    fname = fullfile(adas_path,strcat('ccd96/ccd96_',species_string,'.dat'));  % Charge exchange rate (cm^-3/s)
    if isfile(fname)
        ccd = read_adas_adf11_file(fname);
    else
        ccd = [];
    end
end

function sv_iz = calc_sv_ee_internal(Te,ne,rate,method)
    % Te in eV, ne in m^-3
    % Calculate <sigma*v>_iz
    global perform
    persistent iwarn
    if isempty(iwarn)
        iwarn = 1;
    end
    if perform.LEVEL > 0
        perform.ncall_sv_iz = perform.ncall_sv_iz + 1;
    end
    if method == 1
        sv_iz = 1e-6*10^(log_interp_adas_rate_coefficient_internal(log10(Te),log10(ne./1e6),rate.te_log10,rate.ne_log10,rate.coeff_log10));
    else
        if iwarn == 1
            warning('Using simple form for sv, only appropriate for H ionization!')
            iwarn = 0;
        end
        sv_iz = 3e-16*Te.^2./(3+0.01*Te.^2);
    end
    if any(sv_iz) < 0
        error('Negative sv!!')
    end
    
end




function sigma_v = log_interp_adas_rate_coefficient_internal(Te,ne,te_adas,ne_adas,coeff_adas)
    % Te in eV, ne in cm^-3
    % coefficient cm^3/s
    % above all log10!
    
%     my_method = 'makima';
        my_method = 'spline';
    % my_method = 'linear';
    
    if ne > max(ne_adas)
        warning('Some densities are out of bounds!')
        ne = max(ne_adas);
    end
    extrap_val = NaN;
    sigma_v = interp2(te_adas,ne_adas,coeff_adas,Te,ne,my_method,extrap_val);
end