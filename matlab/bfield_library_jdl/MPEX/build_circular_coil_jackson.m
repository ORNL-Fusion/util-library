function coil = build_circular_coil_jackson(r1,r2,z1,dz,nturns,nlayers)
% layers in radius, turns in axial
% coils centered around r=0
% windings ordered from min Z, min R; R increases fastest
%  r2 -------------------
%     |     |     |     |
%     |  4  |  5  |  6  |        % two layers
%     |     |     |     |
%     -------------------
%     |     |     |     |
%     |  1  |  2  |  3  |        % three turns
%     |     |     |     |
%  r1 -------------------
%    z1                  z1+dz


fw = dz/nturns;
fh = (r2-r1)/nlayers;
zwind = repmat(z1+fw/2:fw:z1+dz-fw/2,1,nlayers);
rwind = repmat(r1+fh/2:fh:r2-fh/2,1,nturns);

coil.rwind = rwind;
coil.zwind = zwind;