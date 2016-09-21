clearvars;


run_path = 'C:\Work\DIII-D\164723\VMEC_XPAND\3059\';
fname = fullfile(run_path,'xpand_164723_3059.dat');

gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
g = readg_g3d(gfile_name);

field = read_xpand_field_file(fname);
field.g = g;


R = [1.8];
Z = [.5];
P = [2*pi];
field_choice = 1
nowarn = 1;
[Br,Bz,Bphi]=bfield_xpand(R,Z,P,field,nowarn,field_choice)
