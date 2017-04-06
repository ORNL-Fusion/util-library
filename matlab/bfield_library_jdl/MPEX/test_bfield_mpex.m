clearvars;


% Set up magnetic field
helicon_current = 440;
current_A = 6400;
current_B = 6400;
current_C = [];
config = 'flat';
verbose = 1;


[coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config,verbose,current_C);
bfield.coil = coil;
bfield.current = current;
bfield.type = 'just_coils';

[coil2,current2] = build_Proto_coils_jackson(helicon_current,current_A,current_B,config,verbose,current_C);
bfield2.coil = coil2;
bfield2.current = current2;
bfield2.type = 'MPEX';

nowarn = 0;
npts = 100;
Z = linspace(0,5,npts);
R = linspace(0.001,.15,npts);
phi_radian = linspace(-pi,pi,npts);

% phi_radian = 0;
% R = 0.05;
% Z = 1.05;
tic;
[Bout,ierr] = bfield_general_rzphi(R,Z,phi_radian,bfield,nowarn);
toc
tic;
[Bout2,ierr2] = bfield_general_rzphi(R,Z,phi_radian,bfield2,nowarn);
toc 

if length(R) <= 10
    for i = 1:length(R)
        fprintf('Evaluating B at R,Z,phi = [%f,%f,%f]\n',R(i),Z(i),phi_radian(i))        
        fprintf('   FILAMENTS:    [Br,Bz,Bphi] = [%e,%e,%e]\n',Bout.br(i),Bout.bz(i),Bout.bphi(i))
        fprintf('   ANALYTIC :    [Br,Bz,Bphi] = [%e,%e,%e]\n',Bout2.br(i),Bout2.bz(i),Bout2.bphi(i))
        fprintf('   Rel Error: %f%%\n',sqrt( (Bout.br(i)-Bout2.br(i))^2 + (Bout.bphi(i)-Bout2.bphi(i))^2 + (Bout.bz(i)-Bout2.bz(i))^2)./sqrt( Bout2.br(i)^2 + Bout2.bz(i)^2)*100)
    end
end

