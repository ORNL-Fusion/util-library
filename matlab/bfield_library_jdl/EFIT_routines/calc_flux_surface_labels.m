function Labels = calc_flux_surface_labels(g,R,Z)
% Map spatial measurements to cached closed-core radial labels.
% Usage: Labels = calc_flux_surface_labels(g,R,Z)
% Build the table first with g = postprocess_gfile_geometry(g). R and Z
% must have matching sizes, or one may be scalar. Labels retains that shape.
% Returns the same fields as lookup_flux_surface_labels, using the same
% interpolation after evaluating psiN(R,Z). A stored-contour check rejects
% private-flux and other points outside the closed-core branch.
% Check Labels.iErr and Labels.errorId before using a returned label.
% Geometry values come from the existing table. Ambiguous spatial points
% require a query-level contour for branch selection; no cache is saved.

Geo = lookup_flux_surface_geometry(g,R,Z);
Labels = select_flux_surface_labels(Geo);
end
