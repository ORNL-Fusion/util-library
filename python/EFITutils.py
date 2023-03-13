import numpy as np


    
def find_xpt(g):
    print('hi')
    
    # Make a first guess based on the boundary info. 
    # This won't work if the configuration is limited, so set 
    # a tolerance to see if the guess is sane.
    firstGuessTol = 1e-3  
    
    b = calc_bfield(g,R,Z)
    #[bpx,ix] = min(bp);    

# Get magnetic field components from gfile, units of Tesla    
# also returns psi and psiN
def calc_bfield(g,R,Z):

    inds = calc_interpolation_inds(g,R,Z)
    ir = inds['ir']
    iz = inds['iz']
    dr = inds['dr']
    dz = inds['dz']
    c = g['psi_bicub_coeffs_inv']    

    psi = calc_psi(g,R,Z,inds=inds)
    psiN = (psi - g['ssimag'])/(g['ssibry'] - g['ssimag'])
    
    dpsidr = ( c['c10'][ir,iz]          + 2*c['c20'][ir,iz]*dr          + 3*c['c30'][ir,iz]*dr*dr    
             + c['c11'][ir,iz]*dz       + 2*c['c21'][ir,iz]*dr*dz       + 3*c['c31'][ir,iz]*dr*dr*dz 
             + c['c12'][ir,iz]*dz*dz    + 2*c['c22'][ir,iz]*dr*dz*dz    + 3*c['c32'][ir,iz]*dr*dr*dz*dz 
             + c['c13'][ir,iz]*dz*dz*dz + 2*c['c23'][ir,iz]*dr*dz*dz*dz + 3*c['c33'][ir,iz]*dr*dr*dz*dz*dz)/g['dR']    

    dpsidz = (   c['c01'][ir,iz]       +   c['c11'][ir,iz]*dr       +   c['c21'][ir,iz]*dr*dr       +   c['c31'][ir,iz]*dr*dr*dr    
             + 2*c['c02'][ir,iz]*dz    + 2*c['c12'][ir,iz]*dr*dz    + 2*c['c22'][ir,iz]*dr*dr*dz    + 2*c['c32'][ir,iz]*dr*dr*dr*dz    
             + 3*c['c03'][ir,iz]*dz*dz + 3*c['c13'][ir,iz]*dr*dz*dz + 3*c['c23'][ir,iz]*dr*dr*dz*dz + 3*c['c33'][ir,iz]*dr*dr*dr*dz*dz)/g['dZ']

    Br = -dpsidz/R
    Bz =  dpsidr/R    
    
    fpol = np.polyval(g['fpol_coeffs'],psiN)
    if psiN <= 1:
        Bt = fpol/R
    else:
        Bt = g['bcentr']*g['rzero']/R
    
    Bpol = np.sqrt(Br**2 + Bz**2)
    Btot = np.sqrt(Br**2 + Bz**2 + Bt**2)
    
    return {'Br':Br,'Bz':Bz,'Bt':Bt,'Bpol':Bpol,'Btor':Btot,'psi':psi,'psiN':psiN}

# Get toroidal magnetic field from gfile
# units of Tesla
def calc_Btor(g,R,Z):
    psiN = calc_psiN(g,R,Z)
    fpol = np.polyval(g['fpol_coeffs'],psiN)
    if psiN <= 1:
        Bt = fpol/R
    else:
        Bt = g['bcentr']*g['rzero']/R
    return Bt


# Get psiN at a point from a gfile
# psiN = (psi(R,Z) - psi(axis))/(psi(boundary) - psi(axis))
# where psi = 2*pi*isign(Ip)*psirz(R,Z), and psirz is the value in the gfile
def calc_psiN(g,R,Z):    
    psi = calc_psi(g,R,Z)
    psiN = (psi - g['ssimag'])/(g['ssibry'] - g['ssimag'])
    return psiN
    
# Get psi at a point from a gfile
# psi = 2*pi*sign(Ip)*psirz(R,Z), where psirz is the value in the gfile
# Can accept inds dict if this was already computed
def calc_psi(g,R,Z,inds=None):

    if inds == None:
        inds = calc_interpolation_inds(g,R,Z)
        
    ir = inds['ir']
    iz = inds['iz']
    dr = inds['dr']
    dz = inds['dz']
    c = g['psi_bicub_coeffs_inv']

    psi = (   c['c00'][ir,iz]          + c['c10'][ir,iz]*dr          + c['c20'][ir,iz]*dr*dr          + c['c30'][ir,iz]*dr*dr*dr 
            + c['c01'][ir,iz]*dz       + c['c11'][ir,iz]*dr*dz       + c['c21'][ir,iz]*dr*dr*dz       + c['c31'][ir,iz]*dr*dr*dr*dz
            + c['c02'][ir,iz]*dz*dz    + c['c12'][ir,iz]*dr*dz*dz    + c['c22'][ir,iz]*dr*dr*dz*dz    + c['c32'][ir,iz]*dr*dr*dr*dz*dz 
            + c['c03'][ir,iz]*dz*dz*dz + c['c13'][ir,iz]*dr*dz*dz*dz + c['c23'][ir,iz]*dr*dr*dz*dz*dz + c['c33'][ir,iz]*dr*dr*dr*dz*dz*dz)

    psi = psi*g['ip_sign']
    return psi

# Helper routine to find relative position in table.
def calc_interpolation_inds(g,R,Z):
    ir = np.floor( (R - g['R'][0])/g['dR'] ).astype(int)
    iz = np.floor( (Z - g['Z'][0])/g['dZ'] ).astype(int)        

    # check for points off grid  (account for derivatives at edge)
    # Grid is already reduced by one on upper side, which is why these are asymmetric
    if (ir <= 1) or (ir >= g['mw'] - 1):
        print('Point off grid in R')

    if (iz <= 1) or (iz >= g['mh'] - 1):
        print('Point off grid in Z')                 

    dr = (R - g['R'][ir])/g['dR']
    dz = (Z - g['Z'][iz])/g['dZ']    
    return {'ir':ir,'iz':iz,'dr':dr,'dz':dz}

#
# Reads a gfile and adds additional quantities for psi and
# bfield interpolation
# Use readg_g3d_simple if you just want the gfile contents
# Note, this applies 2*pi to the psirz grid from the gfile
# JDL
def readg_g3d(filename):

    g = readg_g3d_simple(filename)

    nR = g['mw']
    nZ = g['mh']

    g['dR'] = np.array(g['xdim']/(nR - 1))
    g['dZ'] = np.array(g['zdim']/(nZ - 1))
    g['R'] = np.array([g['rgrid1'] + g['dR']*i for i in range(nR)])
    g['Z'] = np.array([g['zmid'] - 0.5*g['zdim'] + g['dZ']*i for i in range(nZ)])
    g['pn'] = np.array([i/(nZ-1) for i in range(nZ)])

    if g['cpasma'] == []:
        g['ip_sign'] = 1
    else:
        g['ip_sign'] = -np.sign(g['cpasma'])

    psi = g['ip_sign']*g['psirz']

    # Compute bicubic array inverse for divB = 0 interpolation
    # Note that values at edges are not meaningful, points evaluate
    dpsidr = (np.roll(psi,(-1,0),axis=(0,1))-np.roll(psi,(1,0),axis=(0,1)))/(2*g['dR'])
    dpsidz = (np.roll(psi,(0,-1),axis=(0,1))-np.roll(psi,(0,1),axis=(0,1)))/(2*g['dZ'])
    d2psidrdz = (
        (  np.roll(psi,(-1,-1),axis=(0,1)) 
         - np.roll(psi,( 1,-1),axis=(0,1)) 
         - np.roll(psi,(-1, 1),axis=(0,1)) 
         + np.roll(psi,( 1, 1),axis=(0,1))
        )/(4*g['dZ']*g['dR']))

    psi_bicub_coeffs = {}
    psi_bicub_coeffs['c00'] = psi[0:nR-1,0:nZ-1]
    psi_bicub_coeffs['c10'] = dpsidr[0:nR-1,0:nZ-1]*g['dR']
    psi_bicub_coeffs['c20'] = -3*psi[0:nR-1,0:nZ-1] + 3*psi[1:nR,0:nZ-1] - 2*dpsidr[0:nR-1,0:nZ-1]*g['dR'] - dpsidr[1:nR,0:nZ-1]*g['dR']
    psi_bicub_coeffs['c30'] =  2*psi[0:nR-1,0:nZ-1] - 2*psi[1:nR,0:nZ-1] +   dpsidr[0:nR-1,0:nZ-1]*g['dR'] + dpsidr[1:nR,0:nZ-1]*g['dR']

    psi_bicub_coeffs['c01'] = dpsidz[0:nR-1,0:nZ-1]*g['dZ']
    psi_bicub_coeffs['c11'] = d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ']
    psi_bicub_coeffs['c21'] = (-3*dpsidz[0:nR-1,0:nZ-1]*g['dZ'] + 3*dpsidz[1:nR,0:nZ-1]*g['dZ'] 
                               - 2*d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ'] - d2psidrdz[1:nR,0:nZ-1]*g['dR']*g['dZ'])
    psi_bicub_coeffs['c31'] =  (2*dpsidz[0:nR-1,0:nZ-1]*g['dZ'] - 2*dpsidz[1:nR,0:nZ-1]*g['dZ'] 
                                + d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ'] + d2psidrdz[1:nR,0:nZ-1]*g['dR']*g['dZ'])

    psi_bicub_coeffs['c02'] = -3*psi[0:nR-1,0:nZ-1] + 3*psi[0:nR-1,1:nZ] - 2*dpsidz[0:nR-1,0:nZ-1]*g['dZ'] - dpsidz[0:nR-1,1:nZ]*g['dZ']
    psi_bicub_coeffs['c12'] = (-3*dpsidr[0:nR-1,0:nZ-1]*g['dR'] + 3*dpsidr[0:nR-1,1:nZ]*g['dR'] 
                               - 2*d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ'] - d2psidrdz[0:nR-1,1:nZ]*g['dR']*g['dZ'])
    psi_bicub_coeffs['c22'] = ( 9*psi[0:nR-1,0:nZ-1] - 9*psi[1:nR,0:nZ-1] - 9*psi[0:nR-1,1:nZ] + 9*psi[1:nR,1:nZ]
        + 6*dpsidr[0:nR-1,0:nZ-1]*g['dR'] + 3*dpsidr[1:nR,0:nZ-1]*g['dR'] - 6*dpsidr[0:nR-1,1:nZ]*g['dR'] - 3*dpsidr[1:nR,1:nZ]*g['dR'] 
        + 6*dpsidz[0:nR-1,0:nZ-1]*g['dZ'] - 6*dpsidz[1:nR,0:nZ-1]*g['dZ'] + 3*dpsidz[0:nR-1,1:nZ]*g['dZ'] - 3*dpsidz[1:nR,1:nZ]*g['dZ'] 
        + 4*d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ'] + 2*d2psidrdz[1:nR,0:nZ-1]*g['dR']*g['dZ'] 
        + 2*d2psidrdz[0:nR-1,1:nZ]*g['dR']*g['dZ'] + d2psidrdz[1:nR,1:nZ]*g['dR']*g['dZ'])
    psi_bicub_coeffs['c32'] = ( -6*psi[0:nR-1,0:nZ-1] + 6*psi[1:nR,0:nZ-1] + 6*psi[0:nR-1,1:nZ] - 6*psi[1:nR,1:nZ] 
        - 3*dpsidr[0:nR-1,0:nZ-1]*g['dR'] - 3*dpsidr[1:nR,0:nZ-1]*g['dR'] + 3*dpsidr[0:nR-1,1:nZ]*g['dR'] + 3*dpsidr[1:nR,1:nZ]*g['dR'] 
        - 4*dpsidz[0:nR-1,0:nZ-1]*g['dZ'] + 4*dpsidz[1:nR,0:nZ-1]*g['dZ'] - 2*dpsidz[0:nR-1,1:nZ]*g['dZ'] + 2*dpsidz[1:nR,1:nZ]*g['dZ'] 
        - 2*d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ'] - 2*d2psidrdz[1:nR,0:nZ-1]*g['dR']*g['dZ'] 
        - d2psidrdz[0:nR-1,1:nZ]*g['dR']*g['dZ'] - d2psidrdz[1:nR,1:nZ]*g['dR']*g['dZ'])

    psi_bicub_coeffs['c03'] = 2*psi[0:nR-1,0:nZ-1] - 2*psi[0:nR-1,1:nZ] + dpsidz[0:nR-1,0:nZ-1]*g['dZ'] + dpsidz[0:nR-1,1:nZ]*g['dZ']
    psi_bicub_coeffs['c13'] = (2*dpsidr[0:nR-1,0:nZ-1]*g['dR'] - 2*dpsidr[0:nR-1,1:nZ]*g['dR'] 
                               + d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ'] + d2psidrdz[0:nR-1,1:nZ]*g['dR']*g['dZ'])
    psi_bicub_coeffs['c23'] = (-6*psi[0:nR-1,0:nZ-1] + 6*psi[1:nR,0:nZ-1] + 6*psi[0:nR-1,1:nZ] - 6*psi[1:nR,1:nZ] 
        - 4*dpsidr[0:nR-1,0:nZ-1]*g['dR'] - 2*dpsidr[1:nR,0:nZ-1]*g['dR'] + 4*dpsidr[0:nR-1,1:nZ]*g['dR'] + 2*dpsidr[1:nR,1:nZ]*g['dR'] 
        - 3*dpsidz[0:nR-1,0:nZ-1]*g['dZ'] + 3*dpsidz[1:nR,0:nZ-1]*g['dZ'] - 3*dpsidz[0:nR-1,1:nZ]*g['dZ'] + 3*dpsidz[1:nR,1:nZ]*g['dZ'] 
        - 2*d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ'] - d2psidrdz[1:nR,0:nZ-1]*g['dR']*g['dZ'] 
        - 2*d2psidrdz[0:nR-1,1:nZ]*g['dR']*g['dZ'] - d2psidrdz[1:nR,1:nZ]*g['dR']*g['dZ'])
    psi_bicub_coeffs['c33'] =  (4*psi[0:nR-1,0:nZ-1] - 4*psi[1:nR,0:nZ-1] - 4*psi[0:nR-1,1:nZ] + 4*psi[1:nR,1:nZ] 
        + 2*dpsidr[0:nR-1,0:nZ-1]*g['dR'] + 2*dpsidr[1:nR,0:nZ-1]*g['dR'] - 2*dpsidr[0:nR-1,1:nZ]*g['dR'] - 2*dpsidr[1:nR,1:nZ]*g['dR'] 
        + 2*dpsidz[0:nR-1,0:nZ-1]*g['dZ'] - 2*dpsidz[1:nR,0:nZ-1]*g['dZ'] + 2*dpsidz[0:nR-1,1:nZ]*g['dZ'] - 2*dpsidz[1:nR,1:nZ]*g['dZ'] 
        + d2psidrdz[0:nR-1,0:nZ-1]*g['dR']*g['dZ'] + d2psidrdz[1:nR,0:nZ-1]*g['dR']*g['dZ'] 
        + d2psidrdz[0:nR-1,1:nZ]*g['dR']*g['dZ'] + d2psidrdz[1:nR,1:nZ]*g['dR']*g['dZ'])

    g['fpol_coeffs'] = np.polyfit(g['pn'],g['fpol'],7)
    g['psi_bicub_coeffs_inv'] = psi_bicub_coeffs    
    
    return g
    
# Simple reading of gfile
#
# Note: 
#   F = R*Btor
#  Jtor(Amp/m2) = R*P'(psi) + FF'(psi)/R
# JDL
def readg_g3d_simple(filename):
    f = open(filename, "r")
    lines = f.readlines()        
    g = {}
    
    ## Line 1
    # This should follow the formatting, but there are many files that do not.
    # All we really need here are mw and mh, so there are a few attempts to get
    # this from files that do not follow the original format.
    # format (6a8,3i4)

    line = lines[0]
    g['line0'] = line    
    if len(line) < 60:
        # strict formatting not followed
        # Try getting two integers from the end
        splitline = line.split()
        g['mw'] = int(splitline[-2])
        g['mh'] = int(splitline[-1])        
    else:
        # Try formatted read
        g['ecase'] = line[0:8]
        g['mw'] = int(line[52:56])
        g['mh'] = int(line[56:60])
            
    fieldwidth = 16
    
    ## Line 2
    # Format 5e16.9
    line = lines[1];
    g['xdim']   = float(line[0*fieldwidth:(0+1)*fieldwidth])  # rdim. Horizontal dim in meter of comp. box
    g['zdim']   = float(line[1*fieldwidth:(1+1)*fieldwidth])  # zdim  Vertical dim in meter of comp. box
    g['rzero']  = float(line[2*fieldwidth:(2+1)*fieldwidth])  # rcentr R in meter of toroidal magnetic field bcentr
    g['rgrid1'] = float(line[3*fieldwidth:(3+1)*fieldwidth])  # rleft Min R in meter of comp box
    g['zmid']   = float(line[4*fieldwidth:(4+1)*fieldwidth])  # zmid  Z in center of grid in m

    ## Line 3
    # Format 5e16.9
    line = lines[2];
    g['rmaxis'] = float(line[0*fieldwidth:(0+1)*fieldwidth]) # R of magnetic axis in meter
    g['zmaxis'] = float(line[1*fieldwidth:(1+1)*fieldwidth]) # Z of magnetic axis in meter
    g['ssimag'] = float(line[2*fieldwidth:(2+1)*fieldwidth]) # simag Poloidal flux at mag. axis in Weber/rad
    g['ssibry'] = float(line[3*fieldwidth:(3+1)*fieldwidth]) # sibry Poloidal flux at the plasma boundary in Weber/rad
    g['bcentr'] = float(line[4*fieldwidth:(4+1)*fieldwidth]) # Vacuum toroidal mag field in Tesla at RCENTR        

    ## Line 4
    # Some files have empty lines for 4 and 5
    # Format 5e16.9
    line = lines[3];
    if line.strip() == "":
        g['cpasma'] = []
    else:
        g['cpasma'] = float(line[0*fieldwidth:(0+1)*fieldwidth]) # Plasma current in Ampere
    # remaining data is redundant or meaningless (simag,xdum,rmaxis,xdum)

    # Line 5
    # Format 5e16.9
    line = lines[4];
    # all data is redundant or meaningless (zmaxis,xdum,sibry,xdum,xdum)

    ## Read profiles 
    iline = 5     
    
    # Poloidal current function in m-T. F = R*Bt
    array = []      
    icount = 0
    while icount < g['mw']:
        line = lines[iline]
        nfields = len(line) // fieldwidth
        icount2 = 0
        while icount2 < nfields:
            array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
            icount2 += 1
            icount += 1
        iline += 1
    g['fpol'] = np.array(array)

    # Plasma pressure in Nt/m^2 on uniform flux grid
    array = []      
    icount = 0
    while icount < g['mw']:
        line = lines[iline]
        nfields = len(line) // fieldwidth
        icount2 = 0
        while icount2 < nfields:
            array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
            icount2 += 1
            icount += 1
        iline += 1
    g['pres'] = np.array(array)    

    # FF'(psi) in (mT)^2/(Wb/rad) on uniform flux grid
    array = []      
    icount = 0    
    while icount < g['mw']:
        line = lines[iline]
        nfields = len(line) // fieldwidth
        icount2 = 0
        while icount2 < nfields:
            array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
            icount2 += 1
            icount += 1
        iline += 1
    g['ffprim'] = np.array(array)    

    # P'(psi) in (Nt/m^2)/(Wb/rad) on uniform flux grid
    array = []      
    icount = 0
    while icount < g['mw']:
        line = lines[iline]
        nfields = len(line) // fieldwidth
        icount2 = 0
        while icount2 < nfields:
            array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
            icount2 += 1
            icount += 1
        iline += 1
    g['pprime'] = np.array(array)    
    
    # % Poloidal flux in Weber/rad on the rectangular grid points
    array = []      
    icount = 0
    while icount < g['mw']*g['mh']:
        line = lines[iline]
        nfields = len(line) // fieldwidth
        icount2 = 0
        while icount2 < nfields:
            array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
            icount2 += 1
            icount += 1
        iline += 1
    g['psirz'] = np.reshape(np.array(array),(g['mw'],g['mh']),order='F')

    # q values on uniform flux grid from axis to boundary
    array = []      
    icount = 0
    while icount < g['mw']:
        line = lines[iline]
        nfields = len(line) // fieldwidth
        icount2 = 0
        while icount2 < nfields:
            array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
            icount2 += 1
            icount += 1
        iline += 1
    g['qpsi'] = np.array(array)  
    
    if np.equal(g['qpsi'],0).any(): 
        g['qpsi'] = []

    # file ends here in some cases
    if not g['qpsi'] == []:
        # Format 2i5
        line = lines[iline]
        splitline = line.split()
        g['nbdry'] = int(splitline[0])  # nbbbs number of boundary pts
        g['limitr'] = int(splitline[1]) # number of limiter pts              

        # rbbs, zbbs
        array = []      
        icount = 0
        while icount < 2*g['nbdry']:                        
            line = lines[iline]
            nfields = len(line) // fieldwidth
            icount2 = 0
            while icount2 < nfields:
                array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
                icount2 += 1
                icount += 1
            iline += 1
        g['bdry'] = np.reshape(np.array(array),(2,g['nbdry']),order='F')
        g['rbdry'] = g['bdry'][0,:]
        g['zbdry'] = g['bdry'][1,:]
        
        # rlim, zlim
        array = []      
        icount = 0
        while icount < 2*g['limitr']:            
            line = lines[iline]
            nfields = len(line) // fieldwidth
            icount2 = 0
            while icount2 < nfields:
                array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
                icount2 += 1
                icount += 1
            iline += 1
        g['lim'] = np.reshape(np.array(array),(2,g['limitr']),order='F')
        if g['lim'].max() > 100:
            printf('Warning: lim seems to be in [cm], converting to [m]')    
            g['lim'] = g['lim']/100;
        g['rlim'] = g['lim'][0,:]
        g['zlim'] = g['lim'][1,:]            
    else:
        g['nbdry'] = 0
        g['limitr'] = 0
        g['lim'] = []
        g['rlim'] = []
        g['zlim'] = []
        g['bdry'] = []
        g['rbdry'] = []
        g['zbdry'] = []
        
    return g

