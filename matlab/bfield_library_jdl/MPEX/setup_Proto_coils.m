function [nturns,nlayers,rr1,rr2,cl,z0,cur] = setup_Proto_coils(helicon_current,current_A,current_B,config,verbose,current_C)
% helicon_current = coils 3-4
% nturns = number of turns per coil layer in each coil (see note below)
% nlayers = number of layers in each coil (see note below)
% rr1 = inner radius of each coil
% rr2 = outer radius of each coil
% cl = axial length of each coil
% z0 = (axial) starting position of each coil relative to the start of the magnetic field profile at z=0. The coil extends from        this location in the direction of larger z
% cur = the current in each WINDING (amps)
ncoils = 12;
omat = ones(1,ncoils);
nturns  = 8*omat;
nlayers = 5*omat;
rr1 = 0.1221*omat;
rr2 = 0.1785*omat;
cl = 0.0979*omat;
z0  = [0.9392    1.2492    1.5792    1.8152    2.1412    2.3392    2.8952    3.1712    3.3692    3.6852    3.9992    4.3172];
if nargin == 0
%     fprintf('Just returning coil geometry\n')
    cur = [];
    return;
end
if nargin < 6
    current_C = [];
end

if nargin < 3
    error('Must specify currents!')
end
if nargin < 4
    config = 'standard';
    warning('Assuming standard configuration!')
end
if nargin < 5
    verbose = 1;
end


cur = 0*omat;
cur(3:4) = helicon_current;
cur(10:12) = current_B;
switch config
    case 'standard'
        cur([1:2,5:9]) = current_A;
    case 'focus'
        cur([1,5:9]) = current_A;
    case 'flat'
        cur([1,6:9]) = current_A;        
    otherwise
        error(['Did not recognize field configuration: ',config])
end
if ~isempty(current_C)
    cur(2) = current_C;
%     warning(['SETTING CURRENT_C = ',num2str(current_C)]);
end

%     cur = [4000      4000      200       200       4000      4000      4000      0         4000      4000      4000      4000];
%          1         2         3         4         5         6         7         8         9         10        11        12
%     cur = [3300      0         120       120       0         3300      3300      3300      3300      3300      3300      3300  ];
%     cur = [3300      0         210       210       0         3300      3300      3300      3300      0           0         0  ];  % shot 7445

if verbose
    fprintf('--------------------------------------\n')
    fprintf('Coil configuration is %s\n',config)
    fprintf('Using helicon current (coils 3:4) of %f\n',cur(3))
    fprintf('Using current_A of %f\n',current_A)
    fprintf('Using current_B of %f\n',current_B)    
    if ~isempty(current_C)
        fprintf('Using current_C of %f\n',current_C)    
    end
    fprintf('--------------------------------------\n')
end
end