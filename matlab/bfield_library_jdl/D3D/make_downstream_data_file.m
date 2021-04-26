clearvars;

run_path = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX';

outFileName = 'downstream_data_174306.mat';

% From plot_ir_profile
IRfileName = 'IR_profile_174306_twin_3100_4000_elm_80_95.mat';

IR = load(fullfile(run_path,IRfileName)); IR = IR.IR;




DownStream.HeatFlux.IR = IR;

save(fullfile(run_path,outFileName),'DownStream')