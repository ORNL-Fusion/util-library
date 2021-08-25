clearvars;

run_path = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX';


% outFileName = 'downstream_data_174306.mat'; IRFileName = 'IR_profile_174306_twin_3200_4000_elm_80_95.mat'; LPFileName = 'LP_174306_processed.mat'; TSFileName = 'dts_osp_174306.mat';
outFileName = 'downstream_data_174310.mat'; IRFileName = 'IR_profile_174310_twin_3200_4000_elm_80_95.mat'; LPFileName = 'LP_174310_processed.mat'; TSFileName = 'dts_osp_174310.mat';


% From plot_ir_profile
IR = load(fullfile(run_path,IRFileName)); IR = IR.IR;
DownStream.HeatFlux.IR = IR;

% From process_LP
LP = load(fullfile(run_path,LPFileName)); LP = LP.LP;
DownStream.LP = LP;

TS = load(fullfile(run_path,TSFileName)); TS = TS.ts2d;
DownStream.TS = TS;

save(fullfile(run_path,outFileName),'DownStream')