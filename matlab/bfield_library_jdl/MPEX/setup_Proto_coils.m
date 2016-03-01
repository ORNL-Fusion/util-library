function [nturns,nlayers,rr1,rr2,cl,z0,cur] = setup_Proto_coils(helicon_current,current_A,current_B,config,use_split)
% helicon_current = coils 3-4
if nargin < 4
    config = 'standard';
    warning('Assuming standard configuration!'
end
if nargin < 5
    use_split = 0;
end

% nturns = number of turns per coil layer in each coil (see note below)
% nlayers = number of layers in each coil (see note below)
% rr1 = inner radius of each coil
% rr2 = outer radius of each coil
% cl = axial length of each coil
% z0 = (axial) starting position of each coil relative to the start of the magnetic field profile at z=0. The coil extends from        this location in the direction of larger z
% cur = the current magnitudes in each coil  (see note below)


if use_split
    warning('Using split coils: should just be for validation!')
    nturns  = [8     8     8     8     3     3     8     8     3     3     8     8     8     3     3];
    nlayers = [5     5     5     5     7     7     5     5     7     7     5     5     5     7     7];
    rr1 = [0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221    0.1221];
    rr2 = [0.1785    0.1785    0.1785    0.1785    0.2100    0.2100    0.1785    0.1785    0.2100    0.2100    0.1785    0.1785    0.1785    0.2100    0.2100];
    cl  = [0.0979    0.0979    0.0979    0.0979    0.0417    0.0417    0.0979    0.0979    0.0417    0.0417    0.0979    0.0979    0.0979    0.0417    0.0417];
    z0  = [0.9392    1.2492    1.5792    1.8152    2.1412    2.2079    2.3392    2.8952    3.1712    3.2379    3.3692    3.6852    3.9992    4.3172    4.3839];
    cur = [4000        4000         200         200        4000        4000        4000        4000           0           0        4000        4000        4000        4000        4000];
    %                                                       s           s                                     s           s                                             s            s
else
    ncoils = 12;
    omat = ones(1,ncoils);
    nturns  = 8*omat;
    nlayers = 5*omat;
    rr1 = 0.1221*omat;
    rr2 = 0.1785*omat;
    cl = 0.0979*omat;    
    z0  = [0.9392    1.2492    1.5792    1.8152    2.1412    2.3392    2.8952    3.1712    3.3692    3.6852    3.9992    4.3172];
%     cur = [4000      4000      200       200       4000      4000      4000      0         4000      4000      4000      4000];
%          1         2         3         4         5         6         7         8         9         10        11        12
%     cur = [3300      0         120       120       0         3300      3300      3300      3300      3300      3300      3300  ];
    cur = [3300      0         210       210       0         3300      3300      3300      3300      0           0         0  ];  % shot 7445
    if ~isempty(helicon_current)
        cur(3:4) = helicon_current;
    end
    fprintf('--------------------------------------\n')
    fprintf('Using helicon current (coils 3:4) of %f, %f\n',cur(3),cur(4))
    fprintf('--------------------------------------\n')
end