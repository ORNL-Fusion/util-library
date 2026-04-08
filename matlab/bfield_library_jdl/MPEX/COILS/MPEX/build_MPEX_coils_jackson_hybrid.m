function [Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson_hybrid(current_in,simplifyCoils,nturnsSimple,nlayersSimple,verbose)
% [Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson_hybrid(current_in,simplifyCoils,nturnsSimple,nlayersSimple,verbose)
% Builds the MPEX coil set for analytical circular-coil field evaluation.
% Selected coils can be simplified to a coarse equivalent winding set.
%
% simplifyCoils can be specified as:
%   []                    use the default simplification mask
%   logical(1,ncoils)     true for coils to simplify
%   numeric index array   coil indices to simplify
%
% The default is to simplify coils whose center radius exceeds 0.4 m.
% nturnsSimple and nlayersSimple use the usual MPEX definitions:
% turns are axial, layers are radial. The defaults are 3 turns and 2 layers.
% J.D. Lore

if nargin == 0
    error('Must specify inputs')
end
if nargin < 2
    simplifyCoils = [];
end
if nargin < 3
    nturnsSimple = 3;
end
if nargin < 4
    nlayersSimple = 2;
end
if nargin < 5
    verbose = 0;
end

[CoilGeometry,currentPerWinding] = setup_MPEX_coils(current_in,verbose);
simplifyMask = get_simplify_mask(CoilGeometry,simplifyCoils);
nsimple = nturnsSimple*nlayersSimple;

if nturnsSimple <= 0 || nlayersSimple <= 0
    error('nturnsSimple and nlayersSimple must be positive')
end
if mod(nturnsSimple,1) ~= 0 || mod(nlayersSimple,1) ~= 0
    error('nturnsSimple and nlayersSimple must be integers')
end

nwind = CoilGeometry.nturns.*CoilGeometry.nlayers;
nwind_tot = sum(nwind(~simplifyMask)) + nsimple*sum(simplifyMask);
windingCurrent = zeros(1,nwind_tot);
Coil.rwind = zeros(1,nwind_tot);
Coil.zwind = zeros(1,nwind_tot);

i0 = 1;
for i = 1:CoilGeometry.ncoils
    if simplifyMask(i)
        coil_an = build_circular_coil_jackson(CoilGeometry.rr1(i),CoilGeometry.rr2(i),CoilGeometry.z0(i),CoilGeometry.cl(i),nturnsSimple,nlayersSimple);
        i1 = i0;
        i2 = i0 + nsimple - 1;
        Coil.rwind(i1:i2) = coil_an.rwind;
        Coil.zwind(i1:i2) = coil_an.zwind;
        windingCurrent(i1:i2) = currentPerWinding(i)*nwind(i)/nsimple;
        i0 = i2 + 1;
    else
        coil_an = build_circular_coil_jackson(CoilGeometry.rr1(i),CoilGeometry.rr2(i),CoilGeometry.z0(i),CoilGeometry.cl(i),CoilGeometry.nturns(i),CoilGeometry.nlayers(i));
        i1 = i0;
        i2 = i0 + nwind(i) - 1;
        Coil.rwind(i1:i2) = coil_an.rwind;
        Coil.zwind(i1:i2) = coil_an.zwind;
        windingCurrent(i1:i2) = currentPerWinding(i);
        i0 = i2 + 1;
    end
end

if verbose
    fprintf('Hybrid builder simplified %d of %d coils.\n',sum(simplifyMask),CoilGeometry.ncoils)
    fprintf('Simplified coils use %d turns and %d layers.\n',nturnsSimple,nlayersSimple)
    if any(simplifyMask)
        fprintf('Simplified coil indices: ')
        fprintf('%d ',find(simplifyMask))
        fprintf('\n')
    end
end


function simplifyMask = get_simplify_mask(CoilGeometry,simplifyCoils)
ncoils = CoilGeometry.ncoils;

if isempty(simplifyCoils)
    rmid = (CoilGeometry.rr1 + CoilGeometry.rr2)/2;
    simplifyMask = rmid > 0.4;
elseif islogical(simplifyCoils)
    if length(simplifyCoils) ~= ncoils
        error('Logical simplifyCoils array must have length CoilGeometry.ncoils = %d',ncoils)
    end
    simplifyMask = reshape(simplifyCoils,1,ncoils);
elseif isnumeric(simplifyCoils)
    if any(mod(simplifyCoils,1) ~= 0)
        error('Numeric simplifyCoils entries must be integer coil indices')
    end
    if any(simplifyCoils < 1) || any(simplifyCoils > ncoils)
        error('simplifyCoils indices must be between 1 and CoilGeometry.ncoils = %d',ncoils)
    end
    simplifyMask = false(1,ncoils);
    simplifyMask(simplifyCoils) = true;
else
    error('simplifyCoils must be empty, a logical array, or a numeric index array')
end
