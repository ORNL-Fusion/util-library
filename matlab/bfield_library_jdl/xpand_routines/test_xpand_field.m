clearvars;


run_path = 'C:\Work\DIII-D\164723\VMEC_XPAND\3059\';
fname = fullfile(run_path,'xpand_164723_3059.dat');

field = read_xpand_field_file(fname);


R = 2.1;
Z = 0.05;
P = 0.1;
field_choice = 0
nowarn = 1;
[Br,Bz,Bphi]=bfield_xpand(R,Z,P,field,nowarn,field_choice)
