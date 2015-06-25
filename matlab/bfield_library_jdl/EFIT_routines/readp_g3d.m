function p = readp_g3d(filename)

% filename = 'C:\Work\DIII-D\160884\efits\p160884.03014_251';

fid = fopen(filename,'r');

counter = 1; counter_max = 1000;
while ~feof(fid)
    line = fgetl(fid);
    if ~isempty(strfind(line,'ION SPECIES'))
        p.nions = sscanf(line,'%d');
        p.ion_data_NZA = reshape(fscanf(fid,'%f %f %f\n',p.nions*3),3,p.nions).';
    else
        
        p.npts = sscanf(line,'%d');
        if counter > 1
            if p.npts ~= npts_last
                error('Different value for npts, have to account for this')
            end            
        end
        npts_last = p.npts;
        temp = reshape(fscanf(fid,'%f %f %f\n',p.npts*3),3,p.npts).';
        if ~isempty(strfind(line,'psinorm ne(10^20/m^3) dne/dpsiN'))
            p.ne_psin = temp(:,1); p.ne = temp(:,2); p.dne_dpsin = temp(:,3); p.ne_units = '10^20/m^3';
        elseif ~isempty(strfind(line,'psinorm te(KeV) dte/dpsiN'))
            p.te_psin = temp(:,1); p.te = temp(:,2); p.dte_dpsin = temp(:,3); p.te_units = 'keV';
        elseif ~isempty(strfind(line,'psinorm ni(10^20/m^3) dni/dpsiN'))
            p.ni_psin = temp(:,1); p.ni = temp(:,2); p.dni_dpsin = temp(:,3); p.ni_units = '10^20/m^3';
        elseif ~isempty(strfind(line,'psinorm ti(KeV) dti/dpsiN'))
            p.ti_psin = temp(:,1); p.ti = temp(:,2); p.dti_dpsin = temp(:,3); p.ti_units = 'keV';
        elseif ~isempty(strfind(line,'psinorm nb(10^20/m^3) dnb/dpsiN'))
            p.nb_psin = temp(:,1); p.nb = temp(:,2); p.dnb_dpsin = temp(:,3); p.nb_units = '10^20/m^3';
        elseif ~isempty(strfind(line,'psinorm pb(kPa) dpb/dpsiN'))
            p.pb_psin = temp(:,1); p.pb = temp(:,2); p.dpb_dpsin = temp(:,3); p.pb_units = 'kPa';
        elseif ~isempty(strfind(line,'psinorm ptot(kPa) dptot/dpsiN'))
            p.ptot_psin = temp(:,1); p.ptot = temp(:,2); p.dptot_dpsin = temp(:,3); p.ptot_units = 'kPa';
        elseif ~isempty(strfind(line,'psinorm omeg(kRad/s) domeg/dpsiN'))
            p.omeg_psin = temp(:,1); p.omeg = temp(:,2); p.domeg_dpsin = temp(:,3); p.omeg_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm omegp(kRad/s) domegp/dpsiN'))
            p.omegp_psin = temp(:,1); p.omegp = temp(:,2); p.domegp_dpsin = temp(:,3); p.omegp_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm omgvb(kRad/s) domgvb/dpsiN'))
            p.omgvb_psin = temp(:,1); p.omgvb = temp(:,2); p.domgvb_dpsin = temp(:,3); p.omgvb_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm omgpp(kRad/s) domgpp/dpsiN'))
            p.omgpp_psin = temp(:,1); p.omgpp = temp(:,2); p.domgpp_dpsin = temp(:,3); p.omgpp_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm omgeb(kRad/s) domgeb/dpsiN'))
            p.omgeb_psin = temp(:,1); p.omgeb = temp(:,2); p.domgeb_dpsin = temp(:,3); p.omgeb_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm ommvb(kRad/s) dommvb/dpsiN'))
            p.ommvb_psin = temp(:,1); p.ommvb = temp(:,2); p.dommvb_dpsin = temp(:,3); p.ommvb_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm ommpp(kRad/s) dommpp/dpsiN'))
            p.ommpp_psin = temp(:,1); p.ommpp = temp(:,2); p.dommpp_dpsin = temp(:,3); p.ommpp_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm omevb(kRad/s) domevb/dpsiN'))
            p.omevb_psin = temp(:,1); p.omevb = temp(:,2); p.domevb_dpsin = temp(:,3); p.omevb_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm omepp(kRad/s) domepp/dpsiN'))
            p.omepp_psin = temp(:,1); p.omepp = temp(:,2); p.domepp_dpsin = temp(:,3); p.omepp_units = 'kRad/s';
        elseif ~isempty(strfind(line,'psinorm er(kV/m) der/dpsiN'))
            p.er_psin = temp(:,1); p.er = temp(:,2); p.der_dpsin = temp(:,3); p.er_units = 'kV/m';
        elseif ~isempty(strfind(line,'psinorm kpol(km/s/T) dkpol/dpsiN'))
            p.kpol_psin = temp(:,1); p.kpol = temp(:,2); p.dkpol_dpsin = temp(:,3); p.kpol_units = 'km/s/T';        
        elseif ~isempty(strfind(line,'psinorm nz1(10^20/m^3) dnz1/dpsiN'))
            p.nz1_psin = temp(:,1); p.nz1 = temp(:,2); p.dnz1_dpsin = temp(:,3); p.nz1_units = '10^20/m^3';
        elseif ~isempty(strfind(line,'psinorm vtor1(km/s) dvtor1/dpsiN'))
            p.vtor1_psin = temp(:,1); p.vtor1 = temp(:,2); p.vtor1_dpsin = temp(:,3); p.vtor1_units = 'km/s';
        elseif ~isempty(strfind(line,'psinorm vpol1(km/s) dvpol1/dpsiN'))
            p.vpol1_psin = temp(:,1); p.vpol1 = temp(:,2); p.dvpol1_dpsin = temp(:,3); p.vpol1_units = 'km/s';            
        end
        
    end
    
    counter = counter + 1;
    if counter > counter_max
        error('Hit counter_max in readp_g3d')
    end
end
    
fclose(fid);


