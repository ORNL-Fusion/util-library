function [nturns,nlayers,rr1,rr2,cl,z0,cur] = setup_Proto_coils(helicon_current,current_A,current_B,config,verbose,current_C)
% helicon_current = coils 3-4
% cur = the current in each WINDING (amps)
[nturns,nlayers,rr1,rr2,cl,z0] = define_proto_coil_filaments;
if nargin == 0
    fprintf('Just returning coil geometry\n')
    error('Replace with call to define_proto_coil_filaments')
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

ncoils = length(nturns);
omat = ones(1,ncoils);
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