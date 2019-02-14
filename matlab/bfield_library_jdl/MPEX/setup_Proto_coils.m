function [fil,cur] = setup_Proto_coils(current_in,config,verbose)
% [fil,cur] = setup_Proto_coils(current_in,config,verbose)
% Sets the array cur for each winding of the coil set.
%
% current_in can be specified in two ways
%
% If current_in is equal in length to the number of coils in define_proto_coil_filaments (12)
% then the values will be set as the WINDING CURRENT for each coil [Amps].
%
% If length(current_in) < 12 then config must also be specified.
% In this case the elements will be taken as WINDING CURRENTS [Amps] for
% the coils as defined by config. The elements of current_in will be interpreted as
% current_in = [helicon_current, current_A, current_B, current_C] where current_C is
% optional. That is, current_in can be length 3 or 4, see note below.
% The 12 coil winding currents will be set as follows:
%
% config == 'standard'
%   current_in = [A,A,H,H,A,A,A,A,A,B ,B ,B ]
%              1,2,3,4,5,6,7,8,9,10,11,12
% config == 'focus'
%   current_in = [A,0,H,H,A,A,A,A,A,B ,B ,B ]
%              1,2,3,4,5,6,7,8,9,10,11,12
% config == 'flat'
%   current_in = [A,0,H,H,0,A,A,A,A,B ,B ,B ]
%              1,2,3,4,5,6,7,8,9,10,11,12
%
% If current_C is specified, then coil 2 winding current_in will be set to current_C.
%
% verbose ~= 0 will print the coil currents to screen.
% config can be [] when current_in is a 12 element array.

% helicon_current = coils 3-4
% The rest depend on config
% cur = the current_in in each WINDING (amps)
% J.D. Lore
fil = define_proto_coil_filaments;

if nargin == 0
    error('Must specify inputs')
end
if nargin < 2
    config = [];
end
if nargin < 3
    verbose = 0;
end
if ~any(length(current_in) == [3,4,fil.ncoils])
    error('Current array must have length 3, 4, or fil.ncoils = %d\n',fil.ncoils)
end
if any(length(current_in) == [3,4])
    if isempty(config)
        error('If full current array is unset a configuration must be specified')
    end
end

cur = zeros(1,fil.ncoils);
if length(current_in) == fil.ncoils
    cur = current_in;
else
    helicon_current = current_in(1);
    current_A = current_in(2);
    current_B = current_in(3);
    if length(current_in) == 4
        current_C = current_in(4);
    else
        current_C = [];
    end
    
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
    end
end

if verbose
    fprintf('------------------------------------------------------------------------------------------------\n')
    fprintf('Coil currents set to: \n')
    fprintf('%8d',1:fil.ncoils)
    fprintf('\n')
    fprintf('%8.1f',cur)
    fprintf('\n')
    if length(current_in) ~= fil.ncoils        
        fprintf('Coil configuration is %s\n',config)
        fprintf('Using helicon current (coils 3:4) of %f\n',cur(3))
        fprintf('Using current_A of %f\n',current_A)
        fprintf('Using current_B of %f\n',current_B)
        if ~isempty(current_C)
            fprintf('Using current_C of %f\n',current_C)
        end
    end
    fprintf('------------------------------------------------------------------------------------------------\n')
end