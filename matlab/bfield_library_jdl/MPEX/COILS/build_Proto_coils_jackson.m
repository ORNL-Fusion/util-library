function [coil,current] = build_Proto_coils_jackson(current_in,config,verbose)
% [coil,current] = build_Proto_coils_jackson(current_in,config,verbose)
% Builds the coil set and corresponding current array for use in analytical 
% evaluation of the magnetic field and vector potential from circular coils
% corresponding to bfield type 'MPEX'.
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
    
[fil,cur] = setup_Proto_coils(current_in,config,verbose);
nwind = fil.nturns.*fil.nlayers;
nwind_tot = sum(fil.nturns.*fil.nlayers);
current = zeros(1,nwind_tot);
coil.rwind = zeros(1,nwind_tot);
coil.zwind = zeros(1,nwind_tot);
for i = 1:fil.ncoils
    coil_an = build_circular_coil_jackson(fil.rr1(i),fil.rr2(i),fil.z0(i),fil.cl(i),fil.nturns(i),fil.nlayers(i));
    i1 = 1 + sum(nwind(1:i-1));
    i2 = sum(nwind(1:i));
    coil.rwind(i1:i2) = coil_an.rwind;
    coil.zwind(i1:i2) = coil_an.zwind;
    current(i1:i2) = cur(i);
end
