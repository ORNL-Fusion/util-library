function [fil,cur] = setup_MPEX_coils(current_in,verbose)
% [fil,cur] = setup_Proto_coils(current_in,config,verbose)
% Sets the array cur for each winding of the coil set.
%
% current_in can be specified in two ways
%
% If current_in is equal in length to the number of coils in define_MPEX_coil_filaments (21)
% then the values will be set as the WINDING CURRENT for each coil [Amps].
%
% If length(current_in) == 5 then the elements of current_in will be interpreted as
% current_in = [helicon_current, current_ECH, current_ICH, current_transport, current_target]
%
% verbose ~= 0 will print the coil currents to screen.

% current_in = the current_in in each WINDING (amps)
% J.D. Lore
fil = define_MPEX_coil_filaments;

if nargin == 0
    error('Must specify inputs')
end
if nargin < 2
    verbose = 0;
end
if ~any(length(current_in) == [5,fil.ncoils])
    error('Current array must have length 5 or fil.ncoils = %d\n',fil.ncoils)
end

cur = zeros(1,fil.ncoils);
if length(current_in) == fil.ncoils
    cur = current_in;
else
    helicon_current = current_in(1);
    current_ECH = current_in(2);
    current_ICH = current_in(3);
    current_transport = current_in(4);
    current_target = current_in(5);
    
    cur(1:6) = helicon_current;
    cur(7:10) = current_ECH;
    cur(11:15) = current_ICH;
    cur(16:18) = current_transport;
    cur(19:21) = current_target;
end

if verbose
    fprintf('------------------------------------------------------------------------------------------------\n')
    fprintf('Coil currents set to: \n')
    fprintf('%8d',1:fil.ncoils)
    fprintf('\n')
    fprintf('%8.1f',cur)
    fprintf('\n')
    if length(current_in) ~= fil.ncoils        
        fprintf('Using helicon current %f\n',helicon_current)
        fprintf('Using current_ECH of %f\n',current_ECH)
        fprintf('Using current_ICH of %f\n',current_ICH)
        fprintf('Using current_transport of %f\n',current_transport)
        fprintf('Using current_target of %f\n',current_target)
    end
    fprintf('------------------------------------------------------------------------------------------------\n')
end
