% function stress_test_find_xpt
clearvars;


fnames{1} = 'C:\Work\DIII-D\160884\g160884.03400_522';
fnames{2} = 'C:\Work\DIII-D\164723\g164723.03059_410';
fnames{3} = 'C:\Work\DIII-D\165274\kinetic\g165274.02120';
fnames{4} = 'C:\Work_archive\gfile_ARCHIVE\NSTX\g140508.00403';
fnames{5} = 'C:\Work_archive\gfile_ARCHIVE\NSTX\g135183.00433';
fnames{6} = 'C:\Work_archive\CMOD\PSI2016\EXP_DATA\g1150625014.01009_983';
fnames{7} = 'C:\Work_archive\JET\gfiles_corrected\g_p82806_t54.4586';
fnames{8} = 'C:\Work\DIII-D\164723\g164723.03059_C1';
fnames{9} = 'C:\Users\jjl\Dropbox (ORNL)\PSIscidac\WEST\west_54034_10p2s_mag\west_54034_10p2s_mag.eqdsk';
fnames{10} = 'C:\Users\jjl\Dropbox (ORNL)\PSIscidac\ITER\400s.eqdsk';
fnames{11} = 'C:\Users\jjl\Dropbox (ORNL)\PSIscidac\ITER\Baseline2008-li0.70.x4.eqdsk';

for i = 1:length(fnames)
    gfile_name = fnames{i};
    fprintf('%s\n',fnames{i})
    g = readg_g3d(gfile_name,1);
    plot_gfile(g);
end