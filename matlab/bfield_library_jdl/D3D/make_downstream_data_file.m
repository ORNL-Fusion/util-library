clearvars;

run_path = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX';
outFileName = 'downstream_data_174306.mat';

IRFileName = 'IR_profile_174306_twin_3100_4000_elm_80_95.mat';
LPFileName = 'LP_174306_processed.mat';

% From plot_ir_profile
IR = load(fullfile(run_path,IRFileName)); IR = IR.IR;
DownStream.HeatFlux.IR = IR;

% From process_LP
LP = load(fullfile(run_path,LPFileName)); LP = LP.LP;
DownStream.LP = LP;


save(fullfile(run_path,outFileName),'DownStream')