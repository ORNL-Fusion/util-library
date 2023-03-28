import numpy as np

#
# R/Z1, R/Z2 ... start and end points used to calculate surface normal.
# Surface normal is (RZ2-RZ1,0)x(0,0,1) for (R,Z,phi), with phi coming out
# of the page. This means if wall elements are anticlockwise then the
# surface normal will point "out", away from the axis.
#
# npts is number of points to evaluate along this segment. 
# If npts == 1 the midpoint is used.
#
# P.alpha is the angle of incidence relative to the surface, i.e., the 
# incident flux is qdep = qprl*sin(P.alpha)
#
# P.beta is the angle in the poloidal plane between B and the surface
# normal
#
def calc_Bangle_g(g,R1,Z1,R2,Z2,npts=3):
    if npts < 1:
        print('Error: npts must be at least 1')
        return None
    
    if npts == 1:
        Reval = 0.5*(R2 + R1)
        Zeval = 0.5*(Z2 + Z1)
    else:
        Reval = np.linspace(R1,R2,npts)
        Zeval = np.linspace(Z1,Z2,npts)
    
    # This is normalized cross product of v21 = P2 - P1 with a unit vector in the toroidal direction
    # since only a single segment is allowed here, this is constant across the segment
    vNorm = np.array((Zeval[-1] - Z1,-Reval[-1] + R1,0))/np.sqrt((Reval[-1] - R1)**2 + (Zeval[-1] - Z1)**2)
    alphaDeg = np.full((npts,),np.nan)
    betaDeg  = np.full((npts,),np.nan)
    for i in range(npts):
        b = calc_bfield(g,Reval[i],Zeval[i])
        bVec = np.array((b['Br'],b['Bz'],b['Bt']))
        bNorm = bVec/np.linalg.norm(bVec)
        if b['ierr']:
            continue
            
        alphaDeg[i] = np.arcsin(np.dot(vNorm,bNorm))*180/np.pi
        
        bVec[2] = 0
        bNorm = bVec/np.linalg.norm(bVec)
        betaDeg[i] = np.arccos(np.dot(vNorm,bNorm))*180/np.pi
    
    return {'alphaDeg':alphaDeg,'betaDeg':betaDeg}    


def make_target_shape(g,targetAngle,rStart,zStart,forwardLength,backwardLength):
    
    nCircle = 180
    angle = np.zeros((nCircle,))
    r = np.zeros((nCircle,))
    z = np.zeros((nCircle,))
    
    length = forwardLength
    npts = 20
    r[0] = rStart
    z[0] = zStart
    for i in range(npts-1):
        rTest = np.full((npts,),np.nan)
        zTest = np.full((npts,),np.nan)
        for j in range(nCircle):
            b = calc_bfield(g,r[i],z[i])
            if b['ierr']:
                break
            angleOffset = np.arctan2(b['Bz'],b['Br']) + np.pi/2
            rTest[j] = length/npts*np.cos((j-1)*np.pi/nCircle + angleOffset) + r[i]
            zTest[j] = length/npts*np.sin((j-1)*np.pi/nCircle + angleOffset) + z[i]
            angInfo = calc_Bangle_g(g,r[i],z[i],rTest[j],zTest[j])
            angle[j+1] = angInfo['alphaDeg'][1]
        
    
        indMin = np.argmin(
    

# -------------------------------------------------------------------------------------------------------------------------

def move_L_on_C(L,rline,zline):

    dL = np.zeros((rline.size,1))
    dL[1:] = np.sqrt( (rline[0:-1] - rline[1:])**2 + (zline[0:-1] - zline[1:])**2)
    Ltot = np.nansum(dL)
    
    if Ltot < L:
        print('Error: Requested length is less than total line length')
        return None
    
    sumL = np.cumsum(dL)
    if L < 1e-10:
        ind = 0
    else:        
        ind = np.argmax(sumL >= L) - 1
    
    f = (L - sumL[ind])/dL[ind+1]
    R_L = f*(rline[ind+1] - rline[ind]) + rline[ind]
    Z_L = f*(zline[ind+1] - zline[ind]) + zline[ind]
    if f > 0.5:
        icurve_near_L = ind + 1
        err_near_L = (1-f)*dL[ind+1]
    else:
        icurve_near_L = ind
        err_near_L = f*dL[ind+1]
                
    return {'icurve_near_L':icurve_near_L,'err_near_L':err_near_L,'R_L':R_L,'Z_L':Z_L}
# -------------------------------------------------------------------------------------------------------------------------
def follow_fieldlines_rzphi_dl(g,Rstart,Zstart,phistart,dl,nsteps):
    # Number of equations for each ode system
    Neq = 3
    # Number of simultaneous systems to solve
    Nsys = 1
    
    if isinstance(nsteps,float):
        if nsteps.is_integer():
            nsteps = int(nsteps)            
    if not isinstance(nsteps,int):
        print('Error: nsteps must be an integer')
        return None
        
    y = np.empty((Neq*Nsys,))
    y[0::Neq] = Rstart
    y[1::Neq] = phistart
    y[2::Neq] = Zstart

    x = 0
    dx = dl
    
    this = _rk4_fixed_step_integrate_dl(y,x,dx,nsteps,g)
    
    return {'r':this['yout'][:,0::Neq],'phi':this['yout'][:,1::Neq],'z':this['yout'][:,2::Neq],
            'L':this['xout'],'ierr':this['ierr'],
            'iLastGood':this['iLastGood']}
    
def _rk4_fixed_step_integrate_dl(y0,x0,dx,nsteps,g):
    yout = np.full((nsteps+1,y0.size),np.nan)
    xout = np.full((nsteps+1,),np.nan)
    yout[0,] = y0
    xout[0] = x0
    
    y = y0
    x = x0
    for i in range(nsteps):
        dydx = _fl_derivs_dl_gfile(y,g)
        if any(np.isnan(dydx)):
            ierr = 0
            iLastGood = i
            return {'yout':yout,'xout':xout,'ierr':ierr,'iLastGood':iLastGood}
        
        ytmp = _rk4_core_dl(y,dydx,dx,g)
        if any(np.isnan(ytmp)):
            ierr = 0
            iLastGood = i
            return {'yout':yout,'xout':xout,'ierr':ierr,'iLastGood':iLastGood}
        
        x = x + dx
        
        yout[i+1,] = ytmp
        xout[i+1] = x
        y = ytmp
        
    ierr = 0
    iLastGood = nsteps + 1
    return {'yout':yout,'xout':xout,'ierr':ierr,'iLastGood':iLastGood}
                

def _rk4_core_dl(y,dydx0,dx,g):    
    d1 = dx*dydx0
    dydx = _fl_derivs_dl_gfile(y + d1/2,g)
    if any(np.isnan(dydx)):
        return np.full(y.size,np.nan)
    d2 = dx*dydx
    dydx = _fl_derivs_dl_gfile(y + d2/2,g)
    if any(np.isnan(dydx)):
        return np.full(y.size,np.nan)
    d3 = dx*dydx
    dydx = _fl_derivs_dl_gfile(y + d3,g)
    if any(np.isnan(dydx)):
        return np.full(y.size,np.nan)
    d4 = dx*dydx
    
    return y + (d1 + 2*d2 + 2*d3 + d4)/6
    
    
def _fl_derivs_dl_gfile(RPZ,g):
    b = calc_bfield(g,RPZ[0::3],RPZ[2::3])
    df = np.empty((RPZ.size,))
    df[0::3] = b['Br']/b['Btot'] # dR/dl
    df[1::3] = b['Bt']/(RPZ[0::3]*b['Btot']) # dphi/dl
    df[2::3] = b['Bz']/b['Btot'] # dZ/dl
    return df
    
# -------------------------------------------------------------------------------------------------------------------------
# Helper routine to refine psiN for plotting
# -------------------------------------------------------------------------------------------------------------------------
def _refine_psi(g,r,z,fac):
    r2 = np.linspace(g['R'][0],g['R'][-1],fac*g['mw'])
    z2 = np.linspace(g['Z'][0],g['Z'][-1],fac*g['mh'])
    p2 = np.empty((r2.size,z2.size))
    for i in range(r2.size):
        p2[:,i] = calc_psiN(g,r2,z2[i]*np.ones_like(r2))
    return {'r':r2,'z':z2,'psiN':p2}

# -------------------------------------------------------------------------------------------------------------------------
# Find xpt(s) by searching for min(Bpol)
# JDL
# -------------------------------------------------------------------------------------------------------------------------
def find_xpt(g):

    # Tolerance used on refinement
    BpolTol = 1e-6
    
    # Make a first guess based on the boundary info. 
    # This won't work if the configuration is limited, so set 
    # a tolerance to see if the guess is sane.
    # If tolerance exceeded then set a larger search range
    firstGuessBpolTol = 1e-3      
    b = calc_bfield(g,g['rbdry'],g['zbdry'])
    ix = np.argmin(b['Bpol'])
    xpt1 = {'rx':g['rbdry'][ix],'zx':g['zbdry'][ix],'bpx':b['Bpol'][ix]}
    print('Rough 1st x-point:',xpt1['rx'],xpt1['zx'],xpt1['bpx'])
    
    # Setup box for refined search
    dr = (g['R'][-1] - g['R'][0])
    dz = (g['Z'][-1] - g['Z'][0])
    if xpt1['bpx'] <= firstGuessBpolTol:
        dr = dr*.05
        dz = dz*.05
    else:    
        dr = dr*.15
        dz = dz*.10    

    xpt1 = _refine_xpt(g,xpt1,dr,dz,BpolTol)
    print('Refined 1st x-point:',xpt1['rx'],xpt1['zx'],xpt1['bpx'])
    
    # Todo: make sure not off grid, b= None?
    # Guess up/down symmetric for 2nd x-point
    xpt2 = xpt1.copy()
    xpt2['zx'] = -xpt2['zx']
    b = calc_bfield(g,xpt2['rx'],xpt2['zx'])
    xpt2['bpx'] = b['Bpol']
    xpt2 = _refine_xpt(g,xpt2,dr,dz,BpolTol)
    print('Refined 2nd x-point:',xpt2['rx'],xpt2['zx'],xpt2['bpx'])

    return {'xpt1':xpt1,'xpt2':xpt2}

# -------------------------------------------------------------------------------------------------------------------------
# helper routine to refine xpt guess
# -------------------------------------------------------------------------------------------------------------------------
def _refine_xpt(g,xptStart,drStart,dzStart,tol):
                                
    nIterMax = 15    
    nGrid = 100  # grid dimension (square)
    
    Rmin_eval = g['R'][2] + 1e-3
    Zmin_eval = g['Z'][2] + 1e-3
    Rmax_eval = g['R'][-2] - 1e-3
    Zmax_eval = g['Z'][-2] - 1e-3    
    rGrid = np.empty((nGrid,nGrid))
    zGrid = np.empty((nGrid,nGrid))
    bpolGrid = np.empty((nGrid,nGrid))

    # initialize
    xpt = xptStart.copy()
    dr = drStart
    dz = dzStart
                    
    nIter = 0
    while (xpt['bpx'] > tol) and (nIter < nIterMax):      
        rTest = np.linspace(max((Rmin_eval,xpt['rx']-dr)),min((Rmax_eval,xpt['rx']+dr)),nGrid)
        zTest = np.linspace(max((Zmin_eval,xpt['zx']-dz)),min((Zmax_eval,xpt['zx']+dz)),nGrid)
        for i in range(nGrid):            
            ztmp = zTest[i]*np.ones_like(rTest)
            b = calc_bfield(g,rTest,ztmp)
            rGrid[i,:] = rTest
            zGrid[i,:] = ztmp
            bpolGrid[i,:] = b['Bpol']
        
        # find new minimum and set as error
        ix,jx = np.unravel_index(bpolGrid.argmin(), bpolGrid.shape)
        
        xpt['bpx'] = bpolGrid[ix,jx]
        
        rxLast = xpt['rx']
        zxLast = xpt['zx']          
        dr = np.abs(rGrid[ix,jx] - rxLast)
        dz = np.abs(zGrid[ix,jx] - zxLast)
        
        xpt['rx'] = rGrid[ix,jx]
        xpt['zx'] = zGrid[ix,jx]                
        
        #print(nIter,ix,jx,rx,zx,bpx,dr,dz,err > BpolTol)
        
        nIter += 1        
        
        if nIter >= nIterMax:
            print('Warning: max iteration exceeded for 1st x-point')            
        
    return xpt
    
# -------------------------------------------------------------------------------------------------------------------------
# Get magnetic field components from gfile, units of Tesla    
# also returns psi and psiN
# -------------------------------------------------------------------------------------------------------------------------
def calc_bfield(g,R,Z):

    inds = _calc_interpolation_inds(g,R,Z)  

    # this applies g.ip_sign
    psi = calc_psi(g,R,Z,inds=inds)
    psiDerivs = calc_psi_derivs(g,R,Z,inds=inds)
    # so have to apply it again because ssi quantities will be flipped (if ip_sign = -1)
    psiN = (g['ip_sign']*psi - g['ssimag'])/(g['ssibry'] - g['ssimag'])        

    Br = np.where(inds['ierr'] == True, np.nan, -psiDerivs['dpsidz']/R)
    Bz = np.where(inds['ierr'] == True, np.nan,  psiDerivs['dpsidr']/R)
    
    fpol = np.polyval(g['fpol_coeffs'],psiN)
    Bt  = np.where(psiN <= 1, fpol/R, g['bcentr']*g['rzero']/R)
    
    Bpol = np.sqrt(Br**2 + Bz**2)
    Btot = np.sqrt(Br**2 + Bz**2 + Bt**2)
    
    return {'Br':Br,'Bz':Bz,'Bt':Bt,'Bpol':Bpol,'Btot':Btot,'psi':psi,'psiN':psiN,'ierr':inds['ierr']}

# -------------------------------------------------------------------------------------------------------------------------
# Get psiN at a point from a gfile
# psiN = (psi(R,Z) - psi(axis))/(psi(boundary) - psi(axis))
# -------------------------------------------------------------------------------------------------------------------------
def calc_psiN(g,R,Z):
    # this applies g.ip_sign
    psi = calc_psi(g,R,Z)
    # so have to apply it again because ssi quantities will be flipped (if ip_sign = -1)
    psiN = (g['ip_sign']*psi - g['ssimag'])/(g['ssibry'] - g['ssimag'])    
    return psiN

# -------------------------------------------------------------------------------------------------------------------------
# Get psi at a point from a gfile
# psi = -sign(Ip)*psirz(R,Z), where psirz is the value in the gfile
# -------------------------------------------------------------------------------------------------------------------------
def calc_psi(g,R,Z,inds=None):
    
    if inds is None:
        inds = _calc_interpolation_inds(g,R,Z)        
            
    ir = inds['ir']
    iz = inds['iz']
    dr = (R - g['R'][ir])/g['dR']
    dz = (Z - g['Z'][iz])/g['dZ']
    
    c = g['psi_bicub_coeffs_inv']

    psi = g['ip_sign']*(
              c['c00'][ir,iz]          + c['c10'][ir,iz]*dr          + c['c20'][ir,iz]*dr*dr          + c['c30'][ir,iz]*dr*dr*dr 
            + c['c01'][ir,iz]*dz       + c['c11'][ir,iz]*dr*dz       + c['c21'][ir,iz]*dr*dr*dz       + c['c31'][ir,iz]*dr*dr*dr*dz
            + c['c02'][ir,iz]*dz*dz    + c['c12'][ir,iz]*dr*dz*dz    + c['c22'][ir,iz]*dr*dr*dz*dz    + c['c32'][ir,iz]*dr*dr*dr*dz*dz 
            + c['c03'][ir,iz]*dz*dz*dz + c['c13'][ir,iz]*dr*dz*dz*dz + c['c23'][ir,iz]*dr*dr*dz*dz*dz + c['c33'][ir,iz]*dr*dr*dr*dz*dz*dz)

    psi = np.where(inds['ierr'] == True, np.nan, psi)
    return psi

# -------------------------------------------------------------------------------------------------------------------------
# Get psi derivatives at a point from a gfile
# derivatives computed from psi = -sign(Ip)*psirz(R,Z), where psirz is the value in the gfile
# -------------------------------------------------------------------------------------------------------------------------
def calc_psi_derivs(g,R,Z,inds=None):
    
    if inds is None:
        inds = _calc_interpolation_inds(g,R,Z)
            
    ir = inds['ir']
    iz = inds['iz']
    dr = (R - g['R'][ir])/g['dR']
    dz = (Z - g['Z'][iz])/g['dZ']
    
    c = g['psi_bicub_coeffs_inv']
    
    dpsidr = g['ip_sign']*(
               c['c10'][ir,iz]          + 2*c['c20'][ir,iz]*dr          + 3*c['c30'][ir,iz]*dr*dr    
             + c['c11'][ir,iz]*dz       + 2*c['c21'][ir,iz]*dr*dz       + 3*c['c31'][ir,iz]*dr*dr*dz 
             + c['c12'][ir,iz]*dz*dz    + 2*c['c22'][ir,iz]*dr*dz*dz    + 3*c['c32'][ir,iz]*dr*dr*dz*dz 
             + c['c13'][ir,iz]*dz*dz*dz + 2*c['c23'][ir,iz]*dr*dz*dz*dz + 3*c['c33'][ir,iz]*dr*dr*dz*dz*dz)/g['dR']    

    dpsidz = g['ip_sign']*(
                 c['c01'][ir,iz]       +   c['c11'][ir,iz]*dr       +   c['c21'][ir,iz]*dr*dr       +   c['c31'][ir,iz]*dr*dr*dr    
             + 2*c['c02'][ir,iz]*dz    + 2*c['c12'][ir,iz]*dr*dz    + 2*c['c22'][ir,iz]*dr*dr*dz    + 2*c['c32'][ir,iz]*dr*dr*dr*dz    
             + 3*c['c03'][ir,iz]*dz*dz + 3*c['c13'][ir,iz]*dr*dz*dz + 3*c['c23'][ir,iz]*dr*dr*dz*dz + 3*c['c33'][ir,iz]*dr*dr*dr*dz*dz)/g['dZ']    

    dpsidr = np.where(inds['ierr'] == True, np.nan, dpsidr)
    dpsidz = np.where(inds['ierr'] == True, np.nan, dpsidz)
    return {'dpsidr':dpsidr,'dpsidz':dpsidz}

# -------------------------------------------------------------------------------------------------------------------------
# Helper routine to find position in table.
# Points off valid grid region will have ierr = True, ir = iz = 0
# -------------------------------------------------------------------------------------------------------------------------
def _calc_interpolation_inds(g,R,Z):
    # Add a small offset to avoid error with first grid point
    ir = np.asarray(np.floor( (R - g['R'][0] + 1e-10)/g['dR'] ).astype(int))
    iz = np.asarray(np.floor( (Z - g['Z'][0] + 1e-10)/g['dZ'] ).astype(int))

    # check for points off grid, no derivatives on boundary cells
    ierr = np.where((ir < 1) | (ir > g['mw'] - 1) |
                    (iz < 1) | (iz > g['mh'] - 1), True, False)
    ir[ierr == True] = 0
    iz[ierr == True] = 0
    
    return {'ierr':ierr,'ir':ir,'iz':iz}

# -------------------------------------------------------------------------------------------------------------------------
#
# Reads a gfile and adds additional quantities for psi and
# bfield interpolation
# Use readg_g3d_simple if you just want the gfile contents
# JDL
# -------------------------------------------------------------------------------------------------------------------------
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

    g['fpol_coeffs'] = np.polyfit(g['pn'],g['fpol'],7)
    g['psi_bicub_coeffs_inv'] = _get_psi_bicub_coeffs_inv(g)
    
    return g

# -------------------------------------------------------------------------------------------------------------------------
# Simple reading of gfile
#
# Note: 
#   F = R*Btor
#  Jtor(Amp/m2) = R*P'(psi) + FF'(psi)/R
# J.D. Lore
# -------------------------------------------------------------------------------------------------------------------------
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


# -------------------------------------------------------------------------------------------------------------------------
# Calculate bicubic matrix cij such that bicubic interpolation 
#     z(t,u) = sum(i=0:3, sum(j=0:3, cij*ti*uj))
# where
#     t = (x - xi)/(xi+1 - xi)
#     u = (y - yj)/(yj+1 - yj)
# are the relative points on the rectilinear grid. 
# 
# To evaluate z inside of a grid cell, dz/dx, dz/dy, and dz/dxdy are 
# required at each vertex. We use central differences, which will preserve 
# div(B) = 0. The derivatives cannot be calculated for the vertices on the
# boundary, which means that z cannot be evaluated for the boundary cells.
#
# Note: ip_sign = -sign(g.Ip) is ****NOT**** applied to the psirz grid here
#
# Since the bicubic matrix  can be inverted analytically, we use the 
# inverted matrix coefficients directly. 
# 
# J.D. Lore
# -------------------------------------------------------------------------------------------------------------------------
def _get_psi_bicub_coeffs_inv(g):

    nR = g['mw']
    nZ = g['mh']
    dR = g['dR']
    dZ = g['dZ']
    
    psi = g['psirz']

    # Initialize full arrays, only interior points will be evaluated
    dpsidr = np.full([nR,nZ],np.nan)
    dpsidz = np.full([nR,nZ],np.nan)
    d2psidrdz = np.full([nR,nZ],np.nan)
    
    pbci = {}
    pbci['c00'] = np.full([nR,nZ],np.nan)
    pbci['c10'] = np.full([nR,nZ],np.nan)
    pbci['c20'] = np.full([nR,nZ],np.nan)
    pbci['c30'] = np.full([nR,nZ],np.nan)
    pbci['c01'] = np.full([nR,nZ],np.nan)
    pbci['c11'] = np.full([nR,nZ],np.nan)
    pbci['c21'] = np.full([nR,nZ],np.nan)
    pbci['c31'] = np.full([nR,nZ],np.nan)
    pbci['c02'] = np.full([nR,nZ],np.nan)
    pbci['c12'] = np.full([nR,nZ],np.nan)
    pbci['c22'] = np.full([nR,nZ],np.nan)
    pbci['c32'] = np.full([nR,nZ],np.nan)
    pbci['c03'] = np.full([nR,nZ],np.nan)
    pbci['c13'] = np.full([nR,nZ],np.nan)
    pbci['c23'] = np.full([nR,nZ],np.nan)
    pbci['c33'] = np.full([nR,nZ],np.nan)
    
    # Calculate central derivatives on interior points. The 1D derivatives
    # would be valid on the boundaries in the orthogonal direction, but we 
    # keep them as nan for simplicity.

    dpsidr[1:nR-1,1:nZ-1] = (psi[2:nR,1:nZ-1] - psi[0:nR-2,1:nZ-1])/(2*dR);
    dpsidz[1:nR-1,1:nZ-1] = (psi[1:nR-1,2:nZ] - psi[1:nR-1,0:nZ-2])/(2*dZ);
    d2psidrdz[1:nR-1,1:nZ-1] = ( psi[2:nR,2:nZ] - psi[0:nR-2,2:nZ]
                   - psi[2:nR,0:nZ-2] + psi[0:nR-2,0:nZ-2] )/(4*dR*dZ);
    
    # Evaluate coefficients on interior points
    pbci['c00'][1:nR-1,1:nZ-1] = psi[1:nR-1,1:nZ-1]
    pbci['c10'][1:nR-1,1:nZ-1] = dpsidr[1:nR-1,1:nZ-1]*dR
    pbci['c20'][1:nR-1,1:nZ-1] = -3*psi[1:nR-1,1:nZ-1] + 3*psi[2:nR,1:nZ-1] - 2*dpsidr[1:nR-1,1:nZ-1]*dR - dpsidr[2:nR,1:nZ-1]*dR
    pbci['c30'][1:nR-1,1:nZ-1] =  2*psi[1:nR-1,1:nZ-1] - 2*psi[2:nR,1:nZ-1] +   dpsidr[1:nR-1,1:nZ-1]*dR + dpsidr[2:nR,1:nZ-1]*dR

    pbci['c01'][1:nR-1,1:nZ-1] = dpsidz[1:nR-1,1:nZ-1]*dZ
    pbci['c11'][1:nR-1,1:nZ-1] = d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ
    pbci['c21'][1:nR-1,1:nZ-1] = -3*dpsidz[1:nR-1,1:nZ-1]*dZ + 3*dpsidz[2:nR,1:nZ-1]*dZ - 2*d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ - d2psidrdz[2:nR,1:nZ-1]*dR*dZ
    pbci['c31'][1:nR-1,1:nZ-1] =  2*dpsidz[1:nR-1,1:nZ-1]*dZ - 2*dpsidz[2:nR,1:nZ-1]*dZ +   d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ + d2psidrdz[2:nR,1:nZ-1]*dR*dZ

    pbci['c02'][1:nR-1,1:nZ-1] = -3*psi[1:nR-1,1:nZ-1] + 3*psi[1:nR-1,2:nZ] - 2*dpsidz[1:nR-1,1:nZ-1]*dZ - dpsidz[1:nR-1,2:nZ]*dZ
    pbci['c12'][1:nR-1,1:nZ-1] = -3*dpsidr[1:nR-1,1:nZ-1]*dR + 3*dpsidr[1:nR-1,2:nZ]*dR - 2*d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ - d2psidrdz[1:nR-1,2:nZ]*dR*dZ
    pbci['c22'][1:nR-1,1:nZ-1] = ( 9*psi[1:nR-1,1:nZ-1] - 9*psi[2:nR,1:nZ-1] - 9*psi[1:nR-1,2:nZ] + 9*psi[2:nR,2:nZ]
        + 6*dpsidr[1:nR-1,1:nZ-1]*dR + 3*dpsidr[2:nR,1:nZ-1]*dR - 6*dpsidr[1:nR-1,2:nZ]*dR - 3*dpsidr[2:nR,2:nZ]*dR 
        + 6*dpsidz[1:nR-1,1:nZ-1]*dZ - 6*dpsidz[2:nR,1:nZ-1]*dZ + 3*dpsidz[1:nR-1,2:nZ]*dZ - 3*dpsidz[2:nR,2:nZ]*dZ 
        + 4*d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ + 2*d2psidrdz[2:nR,1:nZ-1]*dR*dZ 
        + 2*d2psidrdz[1:nR-1,2:nZ]*dR*dZ + d2psidrdz[2:nR,2:nZ]*dR*dZ)
    pbci['c32'][1:nR-1,1:nZ-1] = ( -6*psi[1:nR-1,1:nZ-1] + 6*psi[2:nR,1:nZ-1] + 6*psi[1:nR-1,2:nZ] - 6*psi[2:nR,2:nZ] 
        - 3*dpsidr[1:nR-1,1:nZ-1]*dR - 3*dpsidr[2:nR,1:nZ-1]*dR + 3*dpsidr[1:nR-1,2:nZ]*dR + 3*dpsidr[2:nR,2:nZ]*dR 
        - 4*dpsidz[1:nR-1,1:nZ-1]*dZ + 4*dpsidz[2:nR,1:nZ-1]*dZ - 2*dpsidz[1:nR-1,2:nZ]*dZ + 2*dpsidz[2:nR,2:nZ]*dZ 
        - 2*d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ - 2*d2psidrdz[2:nR,1:nZ-1]*dR*dZ - d2psidrdz[1:nR-1,2:nZ]*dR*dZ - d2psidrdz[2:nR,2:nZ]*dR*dZ)

    pbci['c03'][1:nR-1,1:nZ-1] = 2*psi[1:nR-1,1:nZ-1] - 2*psi[1:nR-1,2:nZ] + dpsidz[1:nR-1,1:nZ-1]*dZ + dpsidz[1:nR-1,2:nZ]*dZ
    pbci['c13'][1:nR-1,1:nZ-1] = 2*dpsidr[1:nR-1,1:nZ-1]*dR - 2*dpsidr[1:nR-1,2:nZ]*dR + d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ + d2psidrdz[1:nR-1,2:nZ]*dR*dZ
    pbci['c23'][1:nR-1,1:nZ-1] = (-6*psi[1:nR-1,1:nZ-1] + 6*psi[2:nR,1:nZ-1] + 6*psi[1:nR-1,2:nZ] - 6*psi[2:nR,2:nZ] 
        - 4*dpsidr[1:nR-1,1:nZ-1]*dR - 2*dpsidr[2:nR,1:nZ-1]*dR + 4*dpsidr[1:nR-1,2:nZ]*dR + 2*dpsidr[2:nR,2:nZ]*dR 
        - 3*dpsidz[1:nR-1,1:nZ-1]*dZ + 3*dpsidz[2:nR,1:nZ-1]*dZ - 3*dpsidz[1:nR-1,2:nZ]*dZ + 3*dpsidz[2:nR,2:nZ]*dZ 
        - 2*d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ - d2psidrdz[2:nR,1:nZ-1]*dR*dZ - 2*d2psidrdz[1:nR-1,2:nZ]*dR*dZ - d2psidrdz[2:nR,2:nZ]*dR*dZ)
    pbci['c33'][1:nR-1,1:nZ-1] =  (4*psi[1:nR-1,1:nZ-1] - 4*psi[2:nR,1:nZ-1] - 4*psi[1:nR-1,2:nZ] + 4*psi[2:nR,2:nZ] 
        + 2*dpsidr[1:nR-1,1:nZ-1]*dR + 2*dpsidr[2:nR,1:nZ-1]*dR - 2*dpsidr[1:nR-1,2:nZ]*dR - 2*dpsidr[2:nR,2:nZ]*dR 
        + 2*dpsidz[1:nR-1,1:nZ-1]*dZ - 2*dpsidz[2:nR,1:nZ-1]*dZ + 2*dpsidz[1:nR-1,2:nZ]*dZ - 2*dpsidz[2:nR,2:nZ]*dZ 
        + d2psidrdz[1:nR-1,1:nZ-1]*dR*dZ + d2psidrdz[2:nR,1:nZ-1]*dR*dZ + d2psidrdz[1:nR-1,2:nZ]*dR*dZ + d2psidrdz[2:nR,2:nZ]*dR*dZ)
    
    return pbci