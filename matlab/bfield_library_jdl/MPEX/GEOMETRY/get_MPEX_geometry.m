function Geo = get_MPEX_geometry
% This function relies on the currently external private repo
% MPEX-modeling-data where PFCSetup_MPEX.m lives
% J.D. Lore

GeoRaw = PFCSetup_MPEX;

%% Vessel
% Trim vessel data around R = 0 and add closing points on R = 0
check = GeoRaw.VesselData.r < 0;
GeoRaw.VesselData.r(check) = [];
GeoRaw.VesselData.z(check) = [];

% Remove looped point if it exists
tol = 1e-4;
if abs(GeoRaw.VesselData.r(1) - GeoRaw.VesselData.r(end)) < tol && ...
   abs(GeoRaw.VesselData.z(1) - GeoRaw.VesselData.z(end)) < tol
    GeoRaw.VesselData.r(end) = [];
    GeoRaw.VesselData.z(end) = [];
end

% Close the vessel contour on axis for polygon operations
GeoRaw.VesselData.r = [0; GeoRaw.VesselData.r(:); 0];
GeoRaw.VesselData.z = [GeoRaw.VesselData.z(1); GeoRaw.VesselData.z(:); GeoRaw.VesselData.z(end)];


Geo.Vessel.r = GeoRaw.VesselData.r;
Geo.Vessel.z = GeoRaw.VesselData.z;

%% Target
if ~(length(GeoRaw.TargetData.r) == 2 && length(GeoRaw.TargetData.z) == 2)
    error("Unexpected size for TargetData, expected 2 points")
end
Geo.Target.z = GeoRaw.TargetData.z;
Geo.Target.r = GeoRaw.TargetData.r;

% Trim to r > 0
Geo.Target.r(Geo.Target.r < 0) = 0;

%% Target Assembly Shield
if ~(length(GeoRaw.TargetAssemblyShieldData.r) == 2 && length(GeoRaw.TargetAssemblyShieldData.z) == 2)
    error("Unexpected size for TargetAssemblyShieldData, expected 2 points")
end

Geo.TargetAssemblyShield.z = GeoRaw.TargetAssemblyShieldData.z;
Geo.TargetAssemblyShield.r = GeoRaw.TargetAssemblyShieldData.r;

% Trim to r > 0
Geo.TargetAssemblyShield.r(Geo.TargetAssemblyShield.r < 0) = 0;


