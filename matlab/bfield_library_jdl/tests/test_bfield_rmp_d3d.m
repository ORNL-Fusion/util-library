clearvars;

TEST = 2; % 1 =icoils, 2 = ccoils

gfile_name = 'C:\Work\DIII-D\160884\efits\g160884.03014_251';
g = readg_g3d(gfile_name);

if TEST == 1
    taper = [ -3740.0000000000000        3863.0000000000000       -3720.0000000000000        3855.0000000000000       -3718.0000000000000        3858.0000000000000        3862.0000000000000       -3791.0000000000000        3884.0000000000000       -3854.0000000000000        3923.0000000000000       -3847.0000000000000];
    rmp = build_d3d_icoils_jl(taper);
elseif TEST == 2
    taper = [ -3740 3863  -3720  3855 -3718   3858];
    rmp = build_d3d_ccoils_jl(taper);
end

Rtest = [1.8,1.8];
Ztest = [-0.5,0.5];
phitest = [0,0.1];
ntest = length(Rtest);

fprintf('Test call to geq_bicub\n')
b = bfield_geq_bicub(g,Rtest,Ztest);
for i = 1:ntest
    fprintf('Br, Bz, Bphi = [%e,%e,%e]\n',b.br(i),b.bz(i),b.bphi(i))
end

fprintf('Test call to BS xyz\n')
Xtest = Rtest.*cos(phitest);
Ytest = Rtest.*sin(phitest);
[Bx,By,Bz] = bfield_bs_jdl(Xtest,Ytest,Ztest,rmp.coil,rmp.current);
for i = 1:ntest
    fprintf('Bx, By, Bz = [%e,%e,%e]\n',Bx(i),By(i),Bz(i))
end

fprintf('Test call to BS cyl\n')
[Br,Bphi,Bz] = bfield_bs_cyl(Rtest,phitest,Ztest,rmp.coil,rmp.current);
for i = 1:ntest
    fprintf('Br, Bphi, Bz = [%e,%e,%e]\n',Br(i),Bphi(i),Bz(i))
end
