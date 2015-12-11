clearvars;
% time_str indicates the "000001" of files like 2d.000001.003.dat,
% and can also be "final".

%
%  Setup
%

% run_path = 'C:\Work\C2_example\example_jdl\run01\';  grid_dir = [run_path,'ingrid\'];
run_path = 'C:\Work\C2_SOLPS_benchmark\c2_benchmark_v0\case01\'; grid_dir = [run_path,'../ingrid\'];

ndomain = 6;
time_str = 'final';


%
%  Examples
% 

%  --> Load c2 output
grid = read_c2_grid(grid_dir,ndomain);
data = read_c2_data(run_path,ndomain,time_str);

% --> plot grid
plot_c2_grid(grid);

% --> plot output 2d
% data_plot = data.Te2d; mytitle = 'T_e [eV]';
% plot_2d_c2(grid,data_plot,mytitle)

data_plot = data.np2d; mytitle = 'n_e [10^{19} m^-^3]';
plot_2d_c2(grid,data_plot,mytitle)


% --> Example script to plot Te, Ti, np at OMP
plot_omp_profiles_c2(grid,data)