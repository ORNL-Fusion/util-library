function [coil,current] = build_Proto_coils(current_in,config,verbose)
% [coil,current] = build_Proto_coils(current_in,config,verbose)
% Builds the coil set and corresponding current array for use in Biot-Savart 
% filament routines, e.g., bfield type 'just_coils'. Coils that have a very small
% current magnitude are excluded to speed up later calculations.
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
% J.D. Lore
if nargin == 0
    error('Must specify inputs')
end
if nargin < 2 
    config = [];
end
if nargin < 3
    verbose = 0;
end
debug_plots = 0;

[fil,cur] = setup_Proto_coils(current_in,config,verbose);

ntheta_per_wind = 101;  % >= 100 seems to give reasonable results for B components, fewer may be ok for field line following

ibuild = 0;
for i=1:fil.ncoils
    if abs(cur(i)) > 1e-8
        [coil0,current0] = build_circular_coil(fil.rr1(i),fil.rr2(i),fil.z0(i),fil.cl(i),fil.nturns(i),fil.nlayers(i),fil.cur(i),ntheta_per_wind);
        if ibuild == 0
            coil = coil0;
            current = current0;            
        else
            coil = [coil;coil0];
            current = [current;current0];                        
        end
        ibuild = ibuild + 1;
    end
end
if verbose
    fprintf('Built %d out of %d coils. Coils with zero current were excluded.\n',ibuild,length(cur))
end
if debug_plots
    figure; hold on; box on;
    % plot3(coil(:,1),coil(:,2),coil(:,3))
    for i = 1:length(current)-1
        if abs(current(i)) < 1e-8
            plot3(coil(i:i+1,1),coil(i:i+1,2),coil(i:i+1,3),'r')
        else
            plot3(coil(i:i+1,1),coil(i:i+1,2),coil(i:i+1,3),'b')
        end
    end    
end