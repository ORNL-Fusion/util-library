function mex_path = compile_calc_psi_mex()
% Compile calc_psi_mex.c in place and return the mex path.

src = fullfile(fileparts(mfilename('fullpath')),'calc_psi_mex.c');
outdir = fileparts(src);

mex('-outdir',outdir,src);

mex_path = fullfile(outdir,['calc_psi_mex.',mexext]);
