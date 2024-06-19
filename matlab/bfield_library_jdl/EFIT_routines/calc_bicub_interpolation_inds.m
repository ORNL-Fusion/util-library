function [ir,iz,index,ierr] = calc_bicub_interpolation_inds(g,R,Z,nowarn)
% Find the position in the psirz grid and return indices for interpolation

if iscolumn(R)
    R = R';
    Z = Z';
end

% if any(isnan(R)) || any(isnan(Z))
%     error('R or Z is nan')
% end

% Add a small offset to avoid error with first grid point
epsilon = 1e-10;

ir = floor( (R - g.r(1) + epsilon)/g.dR ) + 1;  
iz = floor( (Z - g.z(1) + epsilon)/g.dZ ) + 1;

% Correct for first interior index if it lands on the last grid point
% ir(ir == g.mw - 1) = floor((R(ir == g.mw - 1) - g.r(1) - epsilon) / g.dR) + 1;
% iz(iz == g.mh - 1) = floor((Z(iz == g.mh - 1) - g.z(1) - epsilon) / g.dZ) + 1;

ir(ir == g.mw) = floor((R(ir == g.mw) - g.r(1) - epsilon) / g.dR) + 1;
iz(iz == g.mh) = floor((Z(iz == g.mh) - g.z(1) - epsilon) / g.dZ) + 1;

% check for points off grid, no derivatives on boundary cells
ierr = false(size(R));
% ierr(ir < 2 | ir > g.mw - 1) = true;
% ierr(iz < 2 | iz > g.mh - 1) = true;

ierr(ir < 1 | ir > g.mw) = true;
ierr(iz < 1 | iz > g.mh) = true;

ir(ierr) = nan;
iz(ierr) = nan;

index = ir + g.mw*(iz-1);

if ~nowarn
    for i = find(ierr)
        warning('Point off grid: [R,Z] = [%.2f,%.2f], [Rmin,Zmin] = [%.2f,%.2f], [Rmax,Zmax] = [%.2f,%.2f]\n',R(ierr(i)),Z(ierr(i)),g.r(2),g.r(g.mw-1),g.z(2),g.z(g.mh-1))
    end
end

ierr(isnan(R) | isnan(Z)) = true;