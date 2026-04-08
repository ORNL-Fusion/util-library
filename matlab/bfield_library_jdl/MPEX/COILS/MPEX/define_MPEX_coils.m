function CoilGeometry = define_MPEX_coils
% This function relies on the currently external private repo
% MPEX-modeling-data where CoilSetup_MPEX.m lives
%
% nturns = number of turns per coil layer in each coil (see note below)
% nlayers = number of layers in each coil 
% rr1 = inner radius of each coil
% rr2 = outer radius of each coil
% cl = axial length of each coil
% z0 = (axial) starting position of each coil relative to z=0. The coil extends from        
% this location in the direction of larger z
%
% layers in radius, turns in axial
% coils centered around r=0
%
% Windings are constructed downstream in build_circular_coil*. The
% numbering used there is shown schematically below
%  r2 -------------------
%     |     |     |     |
%     |  2  |  4  |  6  |        % two layers
%     |     |     |     |
%     -------------------
%     |     |     |     |
%     |  1  |  3  |  5  |        % three turns
%     |     |     |     |
%  r1 -------------------
%    z1                  z1+dz
%
%
% JD Lore

data = CoilSetup_MPEX;

CoilGeometry.ncoils = size(data,2);
CoilGeometry.nturns  = data(7,:);
CoilGeometry.nlayers = data(8,:);
CoilGeometry.nwind = CoilGeometry.nturns.*CoilGeometry.nlayers;

CoilGeometry.z0  = data(3,:)./1000;
CoilGeometry.cl  = data(4,:)./1000;
CoilGeometry.rr1 = data(5,:)./1000;
CoilGeometry.rr2 = data(6,:)./1000;

thick = CoilGeometry.rr2 - CoilGeometry.rr1;

CoilGeometry.area = CoilGeometry.cl.*thick;



