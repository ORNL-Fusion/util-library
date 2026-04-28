# -------------------------------------------------------------------------------------------------------------------------
# Collection of bfield routines using EFIT/geqdsk as the bfield method. Adapted from general bfield routines
# so a lot of inefficiency currently. 
#
# JDL
# -------------------------------------------------------------------------------------------------------------------------

import argparse
import copy

import numpy as np

class geq():
    def __init__(self, gfile_name="gfile"):
        self.gfile_name = gfile_name
        pass
    
    def setup_bfield(self):        
        self.read_gfile(self.gfile_name)
        self.find_xpt()
        
    def plot(self,ax):
        psiFine = _refine_psi(self.g,self.g['R'][1:-1],self.g['Z'][1:-1],2)
        ax.plot(self.xpt_info['xpt1']['rx'],self.xpt_info['xpt1']['zx'],'rx')
        ax.contour(psiFine['r'],psiFine['z'],np.transpose(psiFine['psiN']),[1.0])        
        
    def read_gfile(self, gfile_name="gfile"):    
        self.g = readg_g3d(gfile_name)        
        
    # -------------------------------------------------------------------------------------------------------------------------
    # Find xpt(s) by searching for min(Bpol)
    # JDL
    # -------------------------------------------------------------------------------------------------------------------------
    def find_xpt(self):

        verbose = False
        
        # Tolerance used on refinement
        BpolTol = 1e-6

        # Make a first guess based on the boundary info. 
        # This won't work if the configuration is limited, so set 
        # a tolerance to see if the guess is sane.
        # If tolerance exceeded then set a larger search range
        firstGuessBpolTol = 1e-3      
        b = calc_bfield(self.g,self.g['rbdry'],self.g['zbdry'])
        ix = np.argmin(b['Bpol'])
        xpt1 = {'rx':self.g['rbdry'][ix],'zx':self.g['zbdry'][ix],'bpx':b['Bpol'][ix]}
        if verbose:
            print('Rough 1st x-point:',xpt1['rx'],xpt1['zx'],xpt1['bpx'])

        # Setup box for refined search
        dr = (self.g['R'][-1] - self.g['R'][0])
        dz = (self.g['Z'][-1] - self.g['Z'][0])
        if xpt1['bpx'] <= firstGuessBpolTol:
            dr = dr*.05
            dz = dz*.05
        else:    
            dr = dr*.15
            dz = dz*.10    

        xpt1 = _refine_xpt(self.g,xpt1,dr,dz,BpolTol)
        if verbose:
            print('Refined 1st x-point:',xpt1['rx'],xpt1['zx'],xpt1['bpx'])

        # Todo: make sure not off grid, b= None?
        # Guess up/down symmetric for 2nd x-point
        xpt2 = xpt1.copy()
        xpt2['zx'] = -xpt2['zx']
        b = calc_bfield(self.g,xpt2['rx'],xpt2['zx'])
        xpt2['bpx'] = b['Bpol']
        xpt2 = _refine_xpt(self.g,xpt2,dr,dz,BpolTol)
        if verbose:
            print('Refined 2nd x-point:',xpt2['rx'],xpt2['zx'],xpt2['bpx'])

        self.xpt_info = {'xpt1':xpt1,'xpt2':xpt2}
        return self

# -------------------------------------------------------------------------------------------------------------------------
# Calculate angle between magnetic field vector and a line segment in the RZ plane.
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
# P.beta is the angle in the poloidal plane between B and the surface normal
#
# JDL
# -------------------------------------------------------------------------------------------------------------------------
def calc_Bangle_g(g,R1,Z1,R2,Z2,npts=3):
    if npts < 1:
        print('Error: npts must be at least 1')
        return None
    
    if npts == 1:
        Reval = np.array([0.5*(R2 + R1)])
        Zeval = np.array([0.5*(Z2 + Z1)])
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


# -------------------------------------------------------------------------------------------------------------------------
# Driver for fieldline following with parallel distance as integration variable
#
# JDL
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
    
#    if R.size == 1:
#        dpsidr = np.array([dpsidr])
#        dpsidz = np.array([dpsidz])
#        
#    dpsidr[inds['ierr']] = np.nan
#    dpsidz[inds['ierr']] = np.nan
    

    dpsidr = np.where(inds['ierr'] == True, np.nan, dpsidr)
    dpsidz = np.where(inds['ierr'] == True, np.nan, dpsidz)
    return {'dpsidr':dpsidr,'dpsidz':dpsidz}

# -------------------------------------------------------------------------------------------------------------------------
#
# Reads a gfile and adds additional quantities for psi and
# bfield interpolation
# Use readg_g3d_simple if you just want the gfile contents
# JDL
# -------------------------------------------------------------------------------------------------------------------------
def readg_g3d(filename):

    g = readg_g3d_simple(filename)
    g = postprocess_gfile(g)
    
    return g


def postprocess_gfile(g):
    nR = g['mw']
    nZ = g['mh']

    g['dR'] = np.array(g['xdim']/(nR - 1))
    g['dZ'] = np.array(g['zdim']/(nZ - 1))
    g['R'] = np.array([g['rgrid1'] + g['dR']*i for i in range(nR)])
    g['Z'] = np.array([g['zmid'] - 0.5*g['zdim'] + g['dZ']*i for i in range(nZ)])
    g['pn'] = np.array([i/(nR-1) for i in range(nR)])

    if g['cpasma'] == []:
        g['ip_sign'] = 1
    else:
        g['ip_sign'] = -np.sign(g['cpasma'])

    g['fpol_coeffs'] = np.polyfit(g['pn'],g['fpol'],7)
    g['psi_bicub_coeffs_inv'] = _get_psi_bicub_coeffs_inv(g)
    
    return g


def refine_gfile(g, fac=2):
    if fac < 1:
        raise ValueError('refinement factor must be at least 1')

    gNew = copy.deepcopy(g)
    gNew['mw'] = (g['mw'] - 1)*fac + 1
    gNew['mh'] = (g['mh'] - 1)*fac + 1

    nR = gNew['mw']
    nZ = gNew['mh']
    gNew['dR'] = np.array(gNew['xdim']/(nR - 1))
    gNew['dZ'] = np.array(gNew['zdim']/(nZ - 1))
    gNew['R'] = np.array([gNew['rgrid1'] + gNew['dR']*i for i in range(nR)])
    gNew['Z'] = np.array([gNew['zmid'] - 0.5*gNew['zdim'] + gNew['dZ']*i for i in range(nZ)])

    if gNew['cpasma'] == []:
        gNew['ip_sign'] = 1
    else:
        gNew['ip_sign'] = -np.sign(gNew['cpasma'])

    gNew['psirz'] = np.empty((nR, nZ))
    for i, zval in enumerate(gNew['Z']):
        gNew['psirz'][:, i] = gNew['ip_sign']*calc_psi(
            g,
            gNew['R'],
            zval*np.ones_like(gNew['R']),
        )

    x_old = np.arange(g['mw'])
    x_new = np.linspace(0, g['mw'] - 1, nR)
    for field in ['fpol', 'pres', 'ffprim', 'pprime', 'qpsi']:
        if field in g and np.asarray(g[field]).size > 0:
            gNew[field] = np.interp(x_new, x_old, np.asarray(g[field]))

    return postprocess_gfile(gNew)


def _write_gfile_values(fid, vals):
    vals = np.asarray(vals).reshape(-1, order='F')
    for i in range(0, vals.size, 5):
        fid.write(''.join([f'{v:16.9e}' for v in vals[i:i+5]]) + '\n')


def write_gfile(g, filename):
    xdum = 0.0
    line0 = g.get('line0', '')
    txt = line0[:29]
    header = f'{txt:<48}{0:4d}{g["mw"]:4d}{g["mh"]:4d}'

    if np.asarray(g['lim']).shape != (2, g['limitr']):
        raise ValueError('lim shape and limitr do not match')
    if np.asarray(g['bdry']).shape != (2, g['nbdry']):
        raise ValueError('bdry shape and nbdry do not match')

    with open(filename, 'w') as fid:
        fid.write(header + '\n')
        _write_gfile_values(fid, [g['xdim'], g['zdim'], g['rzero'], g['rgrid1'], g['zmid']])
        _write_gfile_values(fid, [g['rmaxis'], g['zmaxis'], g['ssimag'], g['ssibry'], g['bcentr']])
        _write_gfile_values(fid, [g['cpasma'], g['ssimag'], xdum, g['rmaxis'], xdum])
        _write_gfile_values(fid, [g['zmaxis'], xdum, g['ssibry'], xdum, xdum])
        _write_gfile_values(fid, g['fpol'])
        _write_gfile_values(fid, g['pres'])
        _write_gfile_values(fid, g['ffprim'])
        _write_gfile_values(fid, g['pprime'])
        _write_gfile_values(fid, g['psirz'])
        _write_gfile_values(fid, g['qpsi'])
        fid.write(f'{g["nbdry"]:5d}{g["limitr"]:5d}\n')
        _write_gfile_values(fid, g['bdry'])
        _write_gfile_values(fid, g['lim'])

    return filename


def double_gfile_resolution(gfile_name, output_file):
    print(f'Reading {gfile_name}')
    g = readg_g3d(gfile_name)
    print(f'Input resolution: nr={g["mw"]}, nz={g["mh"]}')

    gNew = refine_gfile(g, fac=2)
    print(f'Output resolution: nr={gNew["mw"]}, nz={gNew["mh"]}')
    print(f'Writing {output_file}')
    write_gfile(gNew, output_file)

    return output_file


def double_gfile_resolution_main(argv=None):
    parser = argparse.ArgumentParser(
        description='Double the R and Z resolution of an EFIT gfile/geqdsk.'
    )
    parser.add_argument('gfile', help='Input EFIT gfile/geqdsk path')
    parser.add_argument(
        'output',
        help='Output EFIT gfile/geqdsk path',
    )

    args = parser.parse_args(argv)
    double_gfile_resolution(args.gfile, args.output)
    return 0


def write_gfile_vessel_ogr(gfile_name, output_file='vvfile.ogr'):
    print(f'Reading {gfile_name}')
    g = readg_g3d_simple(gfile_name)

    lim = np.asarray(g['lim'])
    if lim.size == 0:
        raise ValueError(f'No limiter data found in {gfile_name}')
    if lim.shape[0] != 2:
        raise ValueError('Expected limiter data with shape (2, nlim)')

    poly = {'r': lim[0, :], 'z': lim[1, :]}
    poly = close_ogr_polygon(poly, block_id=1, verbose=True)
    lim = np.vstack((poly['r'], poly['z']))

    print(f'Limiter points: nlim={lim.shape[1]}')
    print(f'Writing {output_file}')
    with open(output_file, 'w') as fid:
        for rlim, zlim in lim.T:
            fid.write(f'{rlim*1000.0:16.8e} {zlim*1000.0:16.8e}\n')

    return output_file


def gfile_vessel_main(argv=None):
    parser = argparse.ArgumentParser(
        description='Write vvfile.ogr limiter coordinates from an EFIT gfile/geqdsk.'
    )
    parser.add_argument('gfile', help='Input EFIT gfile/geqdsk path')
    parser.add_argument(
        '-o', '--output',
        default='vvfile.ogr',
        help='Output OGR vessel file path (default: vvfile.ogr)',
    )

    args = parser.parse_args(argv)
    write_gfile_vessel_ogr(args.gfile, args.output)
    return 0


def plot_ogr_files(filenames, label_elements=False):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    for filename in filenames:
        polygons = close_ogr_polygons(order_ogr_polygons(read_ogr(filename)), verbose=True)
        _plot_ogr_polygons(ax, polygons, filename=filename, label_elements=label_elements)

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
    ax.grid(True)
    fig.tight_layout()
    plt.show()


def plot_ogr_main(argv=None):
    parser = argparse.ArgumentParser(
        description='Plot one or more two-column OGR coordinate files.'
    )
    parser.add_argument('files', nargs='+', help='Input OGR coordinate files')
    parser.add_argument(
        '--label_elements',
        action='store_true',
        help='Label adjacent OGR elements as block.element',
    )

    args = parser.parse_args(argv)
    plot_ogr_files(args.files, label_elements=args.label_elements)
    return 0


def read_ogr(filename):
    polygons = []
    coords = []

    def flush_coords():
        if coords:
            arr = np.asarray(coords, dtype=float)*0.001
            polygons.append({'r': arr[:, 0], 'z': arr[:, 1]})
            coords.clear()

    with open(filename, 'r') as fid:
        for line in fid:
            if len(line.strip()) == 0:
                flush_coords()
                continue

            vals = line.split()
            if len(vals) < 2:
                raise ValueError(f'Invalid OGR coordinate line in {filename}: {line.rstrip()}')
            coords.append([float(vals[0]), float(vals[1])])

    flush_coords()
    return polygons


def write_ogr(filename, polygons):
    with open(filename, 'w') as fid:
        for i, poly in enumerate(polygons):
            r = np.asarray(poly['r'])
            z = np.asarray(poly['z'])
            if r.shape != z.shape:
                raise ValueError('OGR r and z arrays must have matching shapes')

            for ri, zi in zip(r, z):
                fid.write(f'{ri*1000.0:12.4f} {zi*1000.0:12.4f}\n')

            if i < len(polygons) - 1:
                fid.write('\n')


def _point_line_distance(points, start, end):
    points = np.asarray(points)
    start = np.asarray(start)
    end = np.asarray(end)
    seg = end - start
    seg_len2 = np.dot(seg, seg)

    if seg_len2 == 0:
        return np.linalg.norm(points - start, axis=1)

    rel = points - start
    t = np.clip(rel @ seg / seg_len2, 0.0, 1.0)
    proj = start + t[:, None]*seg
    return np.linalg.norm(points - proj, axis=1)


def _douglas_peucker_keep(points, tolerance):
    n = points.shape[0]
    keep = np.zeros(n, dtype=bool)
    keep[0] = True
    keep[-1] = True

    def recurse(i0, i1):
        if i1 <= i0 + 1:
            return

        interior = points[i0 + 1:i1]
        dist = _point_line_distance(interior, points[i0], points[i1])
        imax = int(np.argmax(dist))
        dmax = dist[imax]

        if dmax > tolerance:
            ikeep = i0 + 1 + imax
            keep[ikeep] = True
            recurse(i0, ikeep)
            recurse(ikeep, i1)

    recurse(0, n - 1)
    return keep


def _is_closed_polyline(points, tolerance):
    if points.shape[0] < 3:
        return False
    return np.linalg.norm(points[0] - points[-1]) <= tolerance


def simplify_ogr_polygon(poly, tolerance):
    r = np.asarray(poly['r'])
    z = np.asarray(poly['z'])
    points = np.column_stack((r, z))

    if points.shape[0] <= 2:
        return {'r': r.copy(), 'z': z.copy()}

    closed = _is_closed_polyline(points, tolerance=1e-12)
    work = points[:-1] if closed else points

    if work.shape[0] <= 2:
        simplified = work
    elif closed:
        i0 = int(np.argmin(work[:, 0]))
        rolled = np.roll(work, -i0, axis=0)
        rolled_closed = np.vstack((rolled, rolled[0]))
        keep = _douglas_peucker_keep(rolled_closed, tolerance)
        simplified = rolled_closed[keep]
    else:
        keep = _douglas_peucker_keep(work, tolerance)
        simplified = work[keep]

    return {'r': simplified[:, 0], 'z': simplified[:, 1]}


def _max_distance_to_polyline(points, reduced_points):
    if points.shape[0] == 0:
        return 0.0
    if reduced_points.shape[0] == 1:
        return float(np.max(np.linalg.norm(points - reduced_points[0], axis=1)))

    best = np.full(points.shape[0], np.inf)
    for i in range(reduced_points.shape[0] - 1):
        dist = _point_line_distance(points, reduced_points[i], reduced_points[i + 1])
        best = np.minimum(best, dist)
    return float(np.max(best))


def simplify_ogr_file(input_file, output_file, tolerance):
    if tolerance < 0:
        raise ValueError('tolerance must be non-negative')

    print(f'Reading {input_file}')
    polygons = close_ogr_polygons(order_ogr_polygons(read_ogr(input_file)), verbose=True)
    if not polygons:
        raise ValueError(f'No OGR polygons found in {input_file}')

    reduced = []
    total_in = 0
    total_out = 0
    max_error = 0.0

    for i, poly in enumerate(polygons, start=1):
        points = np.column_stack((poly['r'], poly['z']))
        new_poly = simplify_ogr_polygon(poly, tolerance)
        new_points = np.column_stack((new_poly['r'], new_poly['z']))
        err = _max_distance_to_polyline(points, new_points)

        total_in += points.shape[0]
        total_out += new_points.shape[0]
        max_error = max(max_error, err)
        reduced.append(close_ogr_polygon(new_poly, block_id=i, verbose=True))

        print(
            f'  block {i}: points {points.shape[0]} -> {new_points.shape[0]}, '
            f'max_error={err:.6e} m'
        )

    print(f'Total points: {total_in} -> {total_out}')
    print(f'Max error: {max_error:.6e} m')
    print(f'Writing {output_file}')
    write_ogr(output_file, reduced)
    return output_file


def simply_ogr_file_main(argv=None):
    parser = argparse.ArgumentParser(
        description='Simplify/reduce OGR polylines using a distance tolerance.'
    )
    parser.add_argument('input_ogr', help='Input OGR coordinate file')
    parser.add_argument('output_ogr', help='Output simplified OGR coordinate file')
    parser.add_argument(
        'tolerance',
        type=float,
        help='Maximum allowed point-to-simplified-polyline distance in meters',
    )

    args = parser.parse_args(argv)
    simplify_ogr_file(args.input_ogr, args.output_ogr, args.tolerance)
    return 0


def _point_key(point, ndigits=12):
    return tuple(np.round(point, ndigits))


def _order_edges(points, edges):
    key_to_point = {}
    adjacency = {}
    edge_keys = []

    for edge_id, (i0, i1) in enumerate(edges):
        k0 = _point_key(points[i0])
        k1 = _point_key(points[i1])
        key_to_point.setdefault(k0, points[i0])
        key_to_point.setdefault(k1, points[i1])
        edge_keys.append((k0, k1))
        adjacency.setdefault(k0, []).append((k1, edge_id))
        adjacency.setdefault(k1, []).append((k0, edge_id))

    unused = set(range(len(edges)))
    ordered_blocks = []

    while unused:
        endpoints = [
            key for key, nbrs in adjacency.items()
            if len([edge_id for _, edge_id in nbrs if edge_id in unused]) == 1
        ]
        if endpoints:
            start = sorted(endpoints)[0]
        else:
            candidate_edges = [edge_keys[i] for i in unused]
            start = sorted([key for edge in candidate_edges for key in edge])[0]

        chain = [start]
        current = start
        previous = None

        while True:
            choices = [
                (nbr, edge_id) for nbr, edge_id in adjacency[current]
                if edge_id in unused and nbr != previous
            ]
            if not choices:
                choices = [
                    (nbr, edge_id) for nbr, edge_id in adjacency[current]
                    if edge_id in unused
                ]
            if not choices:
                break

            nbr, edge_id = choices[0]
            unused.remove(edge_id)
            chain.append(nbr)
            previous, current = current, nbr

        arr = np.asarray([key_to_point[key] for key in chain])
        ordered_blocks.append({'r': arr[:, 0], 'z': arr[:, 1]})

    return ordered_blocks


def _edges_from_ogr_block(points):
    if points.shape[0] < 2:
        return []

    adjacent_edges = [(i, i + 1) for i in range(points.shape[0] - 1)]
    if points.shape[0] % 2 == 0:
        pair_edges = [(i, i + 1) for i in range(0, points.shape[0], 2)]
        pair_components, pair_max_degree = _edge_graph_stats(points, pair_edges)
        adjacent_components, adjacent_max_degree = _edge_graph_stats(points, adjacent_edges)

        if (
            pair_components < adjacent_components
            or (pair_components == adjacent_components and pair_max_degree < adjacent_max_degree)
        ):
            return pair_edges

    return adjacent_edges


def _edge_graph_stats(points, edges):
    adjacency = {}
    for i0, i1 in edges:
        k0 = _point_key(points[i0])
        k1 = _point_key(points[i1])
        adjacency.setdefault(k0, set()).add(k1)
        adjacency.setdefault(k1, set()).add(k0)

    unseen = set(adjacency)
    components = 0
    while unseen:
        components += 1
        stack = [unseen.pop()]
        while stack:
            key = stack.pop()
            for nbr in adjacency[key]:
                if nbr in unseen:
                    unseen.remove(nbr)
                    stack.append(nbr)

    max_degree = max((len(nbrs) for nbrs in adjacency.values()), default=0)
    return components, max_degree


def order_ogr_file(input_file, output_file):
    print(f'Reading {input_file}')
    polygons = read_ogr(input_file)
    if not polygons:
        raise ValueError(f'No OGR polygons found in {input_file}')

    ordered = close_ogr_polygons(order_ogr_polygons(polygons, verbose=True), verbose=True)
    print(f'Writing {output_file}')
    write_ogr(output_file, ordered)
    return output_file


def order_ogr_polygons(polygons, verbose=False):
    ordered = []
    for i, poly in enumerate(polygons, start=1):
        points = np.column_stack((poly['r'], poly['z']))
        edges = _edges_from_ogr_block(points)
        if not edges:
            ordered.append({'r': points[:, 0], 'z': points[:, 1]})
            if verbose:
                print(f'  block {i}: points={points.shape[0]}, no edges')
            continue

        new_blocks = _order_edges(points, edges)
        ordered.extend(new_blocks)
        if verbose:
            counts = ', '.join(str(block['r'].size) for block in new_blocks)
            print(
                f'  block {i}: points={points.shape[0]}, edges={len(edges)}, '
                f'ordered_blocks={len(new_blocks)} ({counts} points)'
            )

    return ordered


def _points_are_closed(r, z, tolerance=1e-12):
    if r.size < 2:
        return False
    return np.hypot(r[0] - r[-1], z[0] - z[-1]) <= tolerance


def close_ogr_polygon(poly, block_id=None, verbose=False, tolerance=1e-12):
    r = np.asarray(poly['r'])
    z = np.asarray(poly['z'])
    if r.size == 0:
        return {'r': r.copy(), 'z': z.copy()}

    if _points_are_closed(r, z, tolerance=tolerance):
        return {'r': r.copy(), 'z': z.copy()}

    if verbose:
        if block_id is None:
            print('Adding element to close polygon')
        else:
            print(f'Adding element to close polygon for block {block_id}')

    return {
        'r': np.concatenate((r, r[:1])),
        'z': np.concatenate((z, z[:1])),
    }


def close_ogr_polygons(polygons, verbose=False, tolerance=1e-12):
    return [
        close_ogr_polygon(poly, block_id=i, verbose=verbose, tolerance=tolerance)
        for i, poly in enumerate(polygons, start=1)
    ]


def order_ogr_main(argv=None):
    parser = argparse.ArgumentParser(
        description='Order OGR points so each block follows connected polygon/polyline geometry.'
    )
    parser.add_argument('input_ogr', help='Input OGR coordinate file')
    parser.add_argument('output_ogr', help='Output ordered OGR coordinate file')

    args = parser.parse_args(argv)
    order_ogr_file(args.input_ogr, args.output_ogr)
    return 0


def _parse_ogr_element_label(label):
    parts = str(label).split('.')
    if len(parts) == 1:
        return 1, int(parts[0])
    if len(parts) == 2:
        return int(parts[0]), int(parts[1])
    raise ValueError('Element label must be N or block.N')


def refine_ogr_element(input_file, output_file, element_label, nsections=2):
    if nsections < 2:
        raise ValueError('nsections must be at least 2')

    block_id, element_id = _parse_ogr_element_label(element_label)
    print(f'Reading {input_file}')
    polygons = close_ogr_polygons(order_ogr_polygons(read_ogr(input_file)), verbose=True)

    if block_id < 1 or block_id > len(polygons):
        raise ValueError(f'Block {block_id} not found in {input_file}')

    poly = polygons[block_id - 1]
    r = np.asarray(poly['r'])
    z = np.asarray(poly['z'])
    n_elements = r.size - 1

    if element_id < 1 or element_id > n_elements:
        raise ValueError(
            f'Element {block_id}.{element_id} not found; block has {n_elements} elements'
        )

    i0 = element_id - 1
    p0 = np.array([r[i0], z[i0]])
    p1 = np.array([r[i0 + 1], z[i0 + 1]])
    t = np.linspace(0.0, 1.0, nsections + 1)
    new_points = p0[None, :] + t[:, None]*(p1 - p0)
    inserted = new_points[1:-1]

    new_r = np.concatenate((r[:i0 + 1], inserted[:, 0], r[i0 + 1:]))
    new_z = np.concatenate((z[:i0 + 1], inserted[:, 1], z[i0 + 1:]))
    polygons[block_id - 1] = {'r': new_r, 'z': new_z}
    polygons = close_ogr_polygons(order_ogr_polygons(polygons), verbose=True)

    print(
        f'Refined element {block_id}.{element_id}: '
        f'1 section -> {nsections} sections'
    )
    print(f'Block {block_id}: points {r.size} -> {new_r.size}')
    print(f'Writing {output_file}')
    write_ogr(output_file, polygons)
    return output_file


def refine_ogr_elements_main(argv=None):
    parser = argparse.ArgumentParser(
        description='Split one labeled OGR element into multiple equal sections.'
    )
    parser.add_argument('input_ogr', help='Input OGR coordinate file')
    parser.add_argument('output_ogr', help='Output refined OGR coordinate file')
    parser.add_argument('element', help='Element label, either N or block.N')
    parser.add_argument(
        'nsections',
        nargs='?',
        type=int,
        default=2,
        help='Number of sections to split the element into (default: 2)',
    )

    args = parser.parse_args(argv)
    refine_ogr_element(args.input_ogr, args.output_ogr, args.element, args.nsections)
    return 0


def _plot_ogr_labels(ax, r, z, block_id, label_segments=True, label_points=False):
    if label_points:
        for i, (ri, zi) in enumerate(zip(r, z), start=1):
            ax.text(ri, zi, f'{block_id}.p{i}', fontsize=7, ha='center', va='center')

    if label_segments:
        for i in range(len(r) - 1):
            rm = 0.5*(r[i] + r[i + 1])
            zm = 0.5*(z[i] + z[i + 1])
            ax.text(rm, zm, f'{block_id}.{i + 1}', fontsize=8, ha='center', va='center')


def _plot_ogr_polygons(
    ax,
    polygons,
    filename=None,
    label_elements=False,
    label_points=False,
    color=None,
):
    for block_id, poly in enumerate(polygons, start=1):
        r = poly['r']
        z = poly['z']
        if filename is None:
            label = f'block {block_id}'
        else:
            label = filename if len(polygons) == 1 else f'{filename}:{block_id}'

        ax.plot(r, z, '-o', label=label, color=color)
        if label_elements or label_points:
            _plot_ogr_labels(
                ax,
                r,
                z,
                block_id,
                label_segments=label_elements,
                label_points=label_points,
            )


def _plot_structures(ax, structures):
    for structure_id, points in enumerate(structures, start=1):
        r = points[:, 0]
        z = points[:, 1]
        line, = ax.plot(r, z, '-o', label=f'structure {structure_id}')

        if points.shape[0] > 0:
            imid = points.shape[0]//2
            ax.text(
                r[imid],
                z[imid],
                f'S{structure_id}',
                color=line.get_color(),
                fontsize=10,
                fontweight='bold',
                ha='center',
                va='center',
            )

        for i in range(points.shape[0] - 1):
            rm = 0.5*(r[i] + r[i + 1])
            zm = 0.5*(z[i] + z[i + 1])
            ax.text(
                rm,
                zm,
                f'{structure_id}.{i + 1}',
                color=line.get_color(),
                fontsize=8,
                ha='center',
                va='center',
            )


def plot_ogr_for_structure_selection(
    ogr_file,
    output_png=None,
    show=True,
    label_segments=True,
    label_points=False,
):
    import matplotlib.pyplot as plt

    polygons = close_ogr_polygons(order_ogr_polygons(read_ogr(ogr_file)), verbose=True)
    if not polygons:
        raise ValueError(f'No OGR polygons found in {ogr_file}')

    print(f'Reading {ogr_file}')
    print(f'OGR blocks: {len(polygons)}')

    fig, ax = plt.subplots()
    for block_id, poly in enumerate(polygons, start=1):
        r = poly['r']
        z = poly['z']
        nseg = max(r.size - 1, 0)
        print(f'  block {block_id}: points={r.size}, adjacent_segments={nseg}')
    _plot_ogr_polygons(
        ax,
        polygons,
        label_elements=label_segments,
        label_points=label_points,
    )

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(ogr_file)
    ax.grid(True)
    ax.legend()

    if output_png:
        print(f'Writing {output_png}')
        fig.savefig(output_png, dpi=200, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

    return polygons


def plot_ogr_structure_conversion(
    ogr_file,
    polygons,
    structures,
    output_png=None,
    show=False,
    label_segments=True,
    label_points=False,
):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)
    ax_ogr, ax_struct = axes

    _plot_ogr_polygons(
        ax_ogr,
        polygons,
        filename=ogr_file,
        label_elements=label_segments,
        label_points=label_points,
    )
    ax_ogr.set_title('OGR elements')

    _plot_structures(ax_struct, structures)
    ax_struct.set_title('structure.dat groups')

    for ax in axes:
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('R (m)')
        ax.set_ylabel('Z (m)')
        ax.grid(True)
        ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))

    fig.tight_layout()

    if output_png:
        print(f'Writing {output_png}')
        fig.savefig(output_png, dpi=200, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)


def convert_ogr_to_structure_dat_main(argv=None):
    parser = argparse.ArgumentParser(
        description='Convert grouped OGR elements to a GOAT/SOLPS structure.dat file.'
    )
    parser.add_argument('ogr_file', help='Input OGR coordinate file')
    parser.add_argument('structure_dat', help='Output structure.dat file')
    parser.add_argument(
        'groups',
        nargs='*',
        help='Element groups, e.g. "72:74,1" "2:67" "68:71".',
    )
    parser.add_argument(
        '--group',
        action='append',
        help='Element group, e.g. "72:74,1" or "2:67". Repeat for each structure.',
    )
    parser.add_argument('--png', help='Save labeled plot to this PNG file')
    parser.add_argument('--show', action='store_true', help='Show the labeled plot window')
    parser.add_argument('--no-show', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument(
        '--label-points',
        action='store_true',
        help='Also label each point as block.pN',
    )
    parser.add_argument(
        '--no-label-segments',
        action='store_true',
        help='Do not label adjacent segments as block.N',
    )

    args = parser.parse_args(argv)
    group_texts = expand_ogr_group_arguments((args.group or []) + args.groups)
    if not group_texts:
        parser.error('at least one element group is required')

    groups = [parse_ogr_element_group(group) for group in group_texts]
    convert_ogr_to_structure_dat(
        args.ogr_file,
        args.structure_dat,
        groups,
        output_png=args.png,
        show=args.show and not args.no_show,
        label_segments=not args.no_label_segments,
        label_points=args.label_points,
    )
    return 0


def expand_ogr_group_arguments(texts):
    groups = []
    for text in texts:
        text = text.strip()
        if not text:
            continue

        if '[' not in text:
            groups.append(text)
            continue

        current = []
        depth = 0
        for char in text:
            if char == '[':
                if depth == 0:
                    current = []
                else:
                    current.append(char)
                depth += 1
            elif char == ']':
                depth -= 1
                if depth < 0:
                    raise ValueError(f'Unmatched "]" in group argument: {text}')
                if depth == 0:
                    groups.append(''.join(current).strip())
                else:
                    current.append(char)
            elif depth > 0:
                current.append(char)
            elif char in ', ':
                continue
            else:
                groups.append(text)
                break

        if depth != 0:
            raise ValueError(f'Unmatched "[" in group argument: {text}')

    return groups


def parse_ogr_element_group(text):
    text = text.strip()
    if text.startswith('[') and text.endswith(']'):
        text = text[1:-1]

    group = []
    for item in text.split(','):
        item = item.strip()
        if not item:
            continue

        if ':' in item:
            start_text, stop_text = item.split(':', 1)
            start = int(start_text)
            stop = int(stop_text)
            step = 1 if stop >= start else -1
            group.extend(range(start, stop + step, step))
        else:
            group.append(int(item))

    if not group:
        raise ValueError(f'Empty OGR element group: {text}')
    return group


def convert_ogr_to_structure_dat(
    ogr_file,
    structure_dat,
    groups,
    output_png=None,
    show=False,
    label_segments=True,
    label_points=False,
):
    polygons = close_ogr_polygons(order_ogr_polygons(read_ogr(ogr_file)), verbose=True)
    if len(polygons) != 1:
        raise ValueError('Element grouping currently expects exactly one ordered OGR block')

    poly = polygons[0]
    r = np.asarray(poly['r'])
    z = np.asarray(poly['z'])
    n_elements = r.size - 1

    print(f'Reading {ogr_file}')
    print(f'OGR blocks: {len(polygons)}')
    print(f'  block 1: points={r.size}, elements={n_elements}')

    structures = []
    element_counts = np.zeros(n_elements + 1, dtype=int)
    for i, group in enumerate(groups, start=1):
        for elem in group:
            if elem < 1 or elem > n_elements:
                raise ValueError(f'Element {elem} out of range 1:{n_elements}')
            element_counts[elem] += 1

        pts, nbreaks = points_from_element_group(r, z, group)
        structures.append(pts)
        print(
            f'  structure {i}: elements={format_ogr_element_group(group)}, '
            f'points={pts.shape[0]}'
        )
        if nbreaks:
            print(f'    warning: group contains {nbreaks} disconnected chain break(s)')

    repeated = [i for i in range(1, n_elements + 1) if element_counts[i] > 1]
    missing = [i for i in range(1, n_elements + 1) if element_counts[i] == 0]
    print(f'Elements included in more than one structure: {format_ogr_element_list(repeated)}')
    print(f'Elements not included in any structure: {format_ogr_element_list(missing)}')

    print(f'Writing {structure_dat} (open structures: negative point counts)')
    write_structure_dat(structure_dat, structures, closed=False)

    if output_png or show:
        plot_ogr_structure_conversion(
            ogr_file,
            polygons,
            structures,
            output_png=output_png,
            show=show,
            label_segments=label_segments,
            label_points=label_points,
        )

    return structure_dat


def format_ogr_element_group(group):
    return ','.join(str(i) for i in group)


def format_ogr_element_list(elements):
    if not elements:
        return 'none'
    return ','.join(str(i) for i in elements)


def points_from_element_group(r, z, group):
    points = []
    n_elements = r.size - 1
    nbreaks = 0

    for elem in group:
        if elem < 1 or elem > n_elements:
            raise ValueError(f'Element {elem} out of range 1:{n_elements}')
        p0 = np.array([r[elem - 1], z[elem - 1]])
        p1 = np.array([r[elem], z[elem]])

        if not points:
            points.append(p0)
            points.append(p1)
            continue

        last = points[-1]
        if np.linalg.norm(last - p0) <= 1e-12:
            points.append(p1)
        elif np.linalg.norm(last - p1) <= 1e-12:
            points.append(p0)
        else:
            nbreaks += 1
            points.append(p0)
            points.append(p1)

    return _dedupe_consecutive_points(np.asarray(points)), nbreaks


def _dedupe_consecutive_points(points, tolerance=1e-12):
    if points.size == 0:
        return points.reshape(0, 2)

    deduped = [points[0]]
    for point in points[1:]:
        if np.linalg.norm(point - deduped[-1]) > tolerance:
            deduped.append(point)
    return np.asarray(deduped)


def sanitize_structure_points(points, tolerance=1e-12):
    points = _dedupe_consecutive_points(np.asarray(points), tolerance=tolerance)
    if points.shape[0] >= 2 and np.linalg.norm(points[0] - points[-1]) <= tolerance:
        points = points[:-1]
    return points


def write_structure_dat(filename, structures, closed=False):
    with open(filename, 'w') as fid:
        fid.write(f'{len(structures):12d}\n')
        fid.write('$structures\n')

        for i, points in enumerate(structures, start=1):
            points = sanitize_structure_points(points)
            npoints = points.shape[0] if closed else -points.shape[0]
            fid.write(f'Structure    {i}\n')
            fid.write(f'{npoints:12d}\n')
            for r, z in points:
                fid.write(f'{r:20.16f} {z:20.16f}\n')

        fid.write('$end\n')

# -------------------------------------------------------------------------------------------------------------------------
# Helper routine to refine psiN for plotting
# JDL
# -------------------------------------------------------------------------------------------------------------------------
def _refine_psi(g,r,z,fac):
    r2 = np.linspace(g['R'][0],g['R'][-1],fac*g['mw'])
    z2 = np.linspace(g['Z'][0],g['Z'][-1],fac*g['mh'])
    p2 = np.empty((r2.size,z2.size))
    for i in range(r2.size):
        p2[:,i] = calc_psiN(g,r2,z2[i]*np.ones_like(r2))
    return {'r':r2,'z':z2,'psiN':p2}

# -------------------------------------------------------------------------------------------------------------------------
# helper routine to refine xpt guess
# JDL
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
# Helper routine to find position in table.
# Points off valid grid region will have ierr = True, ir = iz = 0
# JDL
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
    while icount < g['mw'] and iline < len(lines):
        line = lines[iline]
        nfields = len(line) // fieldwidth
        icount2 = 0
        while icount2 < nfields:
            array.append(float(line[icount2*fieldwidth:(icount2+1)*fieldwidth]))
            icount2 += 1
            icount += 1
        iline += 1
    g['qpsi'] = np.array(array)

    # file ends here in some cases
    if g['qpsi'].size == g['mw']:
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
# RK4 integration step with parallel distance as integration variable
#
# JDL
# -------------------------------------------------------------------------------------------------------------------------
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
                
# -------------------------------------------------------------------------------------------------------------------------
# RK4 core fieldline following with parallel distance as integration variable
#
# JDL
# -------------------------------------------------------------------------------------------------------------------------
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
    
# -------------------------------------------------------------------------------------------------------------------------
# Derivative function for fieldline following with parallel distance as integration variable
#
# JDL
# -------------------------------------------------------------------------------------------------------------------------
def _fl_derivs_dl_gfile(RPZ,g):
    b = calc_bfield(g,RPZ[0::3],RPZ[2::3])
    df = np.empty((RPZ.size,))
    df[0::3] = b['Br']/b['Btot'] # dR/dl
    df[1::3] = b['Bt']/(RPZ[0::3]*b['Btot']) # dphi/dl
    df[2::3] = b['Bz']/b['Btot'] # dZ/dl
    return df


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
