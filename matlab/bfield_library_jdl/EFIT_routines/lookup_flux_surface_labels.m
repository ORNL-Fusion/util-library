function Labels = lookup_flux_surface_labels(g,psiN)
% Map normalized poloidal flux to cached closed-core radial labels.
% Usage: Labels = lookup_flux_surface_labels(g,psiN)
% Build the table first with g = postprocess_gfile_geometry(g). psiN may
% be an array; Labels is a struct array of the same size. This psiN-only
% query assumes the closed-core branch. Use calc_flux_surface_labels for
% spatial measurements requiring a branch check.
% Returns psiN, x, phiTor, rTor, rSaarelma, dpsiNdx, dPhiTorDpsiN,
% drTorDx, drSaarelmaDx, b0Reference, r0Reference, toroidalFieldSign,
% zOMP, the coordinate-definition strings, iErr and errorId.
% phiTor is a positive flux magnitude in Wb; x and the radii are in m.
% Interior queries must lie strictly inside the sampled table interval.
% Exact psiN=0 and 1 also return endpoint labels, with derivatives NaN.
% At the axis x is defined only if zOMP=zmaxis. Unsampled interior gaps
% remain unavailable. rhoTor/rhoSaarelma normalize to rTorSep/rSaarelmaSep;
% xTor/xSaarelma subtract those radii. Check separatrixIErr before using
% normalization or shifts, even when iErr=0 for the interior geometry.

Geo = lookup_flux_surface_table(g,psiN);
Table = g.fluxSurfaceGeometry;
for iPoint = 1:numel(Geo)
    if psiN(iPoint) == 0
        Geo(iPoint).phiTor = 0;
        Geo(iPoint).rTor = 0;
        Geo(iPoint).rSaarelma = 0;
        Geo(iPoint).iErr = 0;
        Geo(iPoint).errorId = '';
        if Table.zOMP == g.zmaxis
            Axis = compute_flux_surface_geometry(g,g.rmaxis,g.zmaxis,Table.refinement,Table.zOMP);
            Geo(iPoint).x = Axis.x;
        end
    elseif psiN(iPoint) == 1
        if Table.separatrixIErr == 0
            Geo(iPoint).phiTor = Table.phiTorSep;
            Geo(iPoint).rTor = Table.rTorSep;
            Geo(iPoint).rSaarelma = Table.rSaarelmaSep;
            Geo(iPoint).x = 0;
            Geo(iPoint).iErr = 0;
            Geo(iPoint).errorId = '';
        else
            Geo(iPoint).errorId = Table.separatrixErrorId;
        end
    end
    Geo(iPoint).rhoTor = Geo(iPoint).rTor/Table.rTorSep;
    Geo(iPoint).rhoSaarelma = Geo(iPoint).rSaarelma/Table.rSaarelmaSep;
    Geo(iPoint).xTor = Geo(iPoint).rTor-Table.rTorSep;
    Geo(iPoint).xSaarelma = Geo(iPoint).rSaarelma-Table.rSaarelmaSep;
end
Labels = select_flux_surface_labels(Geo);
end
