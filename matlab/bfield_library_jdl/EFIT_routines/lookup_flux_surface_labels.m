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
% Queries must lie strictly inside the sampled table interval. The axis,
% separatrix and unsampled regions are not extrapolated. iErr = 1 marks
% unavailable labels, with InvalidPoint, OutsideTable or InvalidInterpolation.

Geo = lookup_flux_surface_table(g,psiN);
Labels = select_flux_surface_labels(Geo);
end
