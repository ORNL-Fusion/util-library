clearvars;

% current_C_vals = linspace(0,3300,10);
current_C_vals = -200;

debug_plots = 1;

skimmer = 1;
geo = get_Proto_geometry(debug_plots,debug_plots,skimmer);

config = 'standard';
helicon_current = 0;
current_A = 3300; 
current_B = 0;
% current_C = 0;

verbose = 1;


for ic = 1:length(current_C_vals)
    fprintf('Working on ic %d of %d\n',ic,length(current_C_vals))
    current_C = current_C_vals(ic);
    
%     [coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config,verbose,current_C);    
%     bfield.coil = coil;
%     bfield.current = current;
%     bfield.type = 'just_coils';
    
    shot_tmp{1} = helicon_current;
    shot_tmp{2} = current_A;
    shot_tmp{3} = current_B;
    shot_tmp{4} = config;
    shot_tmp{5} = skimmer;
    shot_tmp{6} = current_C;
    f = find_lcfs(shot_tmp,0,10);
    lcfs_info{ic}.f = f;
    if debug_plots
        plot(f.z,f.r)
%         plot(f.hit_rz(2),f.hit_rz(1),'.')
    end
end