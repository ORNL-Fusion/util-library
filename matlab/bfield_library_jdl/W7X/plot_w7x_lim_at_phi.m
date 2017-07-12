function plot_w7x_lim_at_phi(phi_rad,newfig)
if nargin < 2
newfig = 0;
end
run_info.run_path = 'C:\Work\util-library\matlab\bfield_library_jdl\W7X\geo';
lim = load_all_limiter_files(run_info);
if newfig
    figure; hold on; box on;
end
plot_emc3_plates_at_phi(lim,phi_rad);


