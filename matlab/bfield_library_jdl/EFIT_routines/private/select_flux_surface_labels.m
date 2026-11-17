function Labels = select_flux_surface_labels(Geo)
% Preserve the common radial-label fields and per-point lookup status.
fields = {'phiTorSep','rTorSep','rSaarelmaSep','separatrixRelativeChange','separatrixRefinement','separatrixIErr','separatrixErrorId','rhoTor','rhoSaarelma','xTor','xSaarelma','psiN','x','phiTor','rTor','rSaarelma','dpsiNdx','dPhiTorDpsiN','drTorDx','drSaarelmaDx','b0Reference','r0Reference','toroidalFieldSign','zOMP','xDefinition','phiTorDefinition','rTorDefinition','rSaarelmaDefinition','iErr','errorId'};
Labels = rmfield(Geo,setdiff(fieldnames(Geo),fields));
end
