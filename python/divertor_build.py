import numpy as np 
import EFITutils as EFITutils
import GeoUtils as GeoUtils

class DivertorBuild():
    def __init__(self):
        self.dLegOuter = 0.7        # Move this distance along outer sep leg before starting target
        self.targetAngleOuter = 2   # Desired angle in deg
        self.SOLLengthOuter = 0.6   # Desired length of target in SOL direction
        self.PFRLengthOuter = 0.3   # Desired length of target in SOL direction

        self.dLegInner= 0.7         # Move this distance along inner sep leg before starting target
        self.targetAngleInner = -2
        self.SOLLengthInner = 0.6
        self.PFRLengthInner = 0.6  # (currently just needs to be large enough that intersection is found)        
        
        self.dL_ff = 0.01  # step size for RK4 fieldline following (m)
        self.LmaxOuter = 100 # Max distance to trace fieldline from xpoint to target on outer side
        self.LmaxInner = 100  # Max distance to trace fieldline from xpoint to target on inner side
        self.LmaxPFR   = 20  # Max distance to trace fieldline across PFR
        self.LmaxOW    = 40  # Max distance to trace fieldline for outer wall
        self.LmaxIW    = 60  # Max distance to trace fieldline for inner wall
                        
        pass
    
    def plot(self,ax):
        ax.plot(self.pfc['r'],self.pfc['z'])
    
        
    def make_divertors(self,geq,ax=None):
        
        if ax is None:
            plotit = False
        else:
            plotit = True
                    
        
        # Outer leg  ----------------------------------------------------------------------------------
        print("Moving ",self.dLegOuter," m along outer leg")
        dL = self.dL_ff
        Lmax = self.LmaxOuter
        Rstart = geq.xpt_info['xpt1']['rx'] + 1e-3
        Zstart = geq.xpt_info['xpt1']['zx']
        s = EFITutils.follow_fieldlines_rzphi_dl(geq.g,Rstart,Zstart,0,dL,abs(Lmax/dL))
        if plotit:
            ax.plot(s['r'],s['z'],'r-')        
        
        move = GeoUtils.move_L_on_C(self.dLegOuter,s['r'][0:s['iLastGood']+1],s['z'][0:s['iLastGood']+1] )
        ROSP = move['R_L']
        ZOSP = move['Z_L']
        if plotit:
            ax.plot(ROSP,ZOSP,'go')        

        print("Making outer target initial guess")
        OTf = make_target_shape_oneway(geq.g,self.targetAngleOuter,ROSP,ZOSP,self.SOLLengthOuter,1)
        OTb = make_target_shape_oneway(geq.g,self.targetAngleOuter,ROSP,ZOSP,self.PFRLengthOuter,-1)
        if plotit:
            ax.plot(OTf['r'],OTf['z'])
            ax.plot(OTb['r'],OTb['z'])
        
        OT = {}
        OT['r'] = np.concatenate((np.flip(OTb['r'][0:OTb['iLastGood']]),OTf['r'][0:OTf['iLastGood']+1]))
        OT['z'] = np.concatenate((np.flip(OTb['z'][0:OTb['iLastGood']]),OTf['z'][0:OTf['iLastGood']+1]))
        OT['s'] = np.concatenate((np.flip(OTb['s'][0:OTb['iLastGood']]),OTf['s'][0:OTf['iLastGood']+1]))
        
        # Inner leg  ----------------------------------------------------------------------------------
        print("Moving ",self.dLegInner," m along inner leg")
        Lmax = self.LmaxInner
        Rstart = geq.xpt_info['xpt1']['rx'] - 1e-3
        Zstart = geq.xpt_info['xpt1']['zx']
        s = EFITutils.follow_fieldlines_rzphi_dl(geq.g,Rstart,Zstart,0,-dL,abs(Lmax/dL))
        if plotit:
            ax.plot(s['r'],s['z'],'r-')
        
        move = GeoUtils.move_L_on_C(self.dLegInner,s['r'][0:s['iLastGood']+1],s['z'][0:s['iLastGood']+1] )
        RISP = move['R_L']
        ZISP = move['Z_L']
        if plotit:
            ax.plot(RISP,ZISP,'go')

        print("Making inner target initial guess")
        ITf = make_target_shape_oneway(geq.g,self.targetAngleInner,RISP,ZISP,self.SOLLengthInner,-1)
        ITb = make_target_shape_oneway(geq.g,self.targetAngleInner,RISP,ZISP,self.PFRLengthInner,1)        

        IT = {}
        IT['r'] = np.concatenate((np.flip(ITb['r'][0:ITb['iLastGood']]),ITf['r'][0:ITf['iLastGood']+1]))
        IT['z'] = np.concatenate((np.flip(ITb['z'][0:ITb['iLastGood']]),ITf['z'][0:ITf['iLastGood']+1]))
        IT['s'] = np.concatenate((np.flip(ITb['s'][0:ITb['iLastGood']]),ITf['s'][0:ITf['iLastGood']+1]))      
        if plotit:
            ax.plot(ITf['r'],ITf['z'])
            ax.plot(ITb['r'],ITb['z'])        
        
        self.IT = IT
        self.OT = OT                
        
    def connect_divertors(self,geq,ax=None):        
        if ax is None:
            plotit = False
        else:
            plotit = True
            
        print("Connecting divertors in PFR")
        # Move along PFR to inner leg
        Lmax = self.LmaxPFR
        dL = self.dL_ff
        Rstart = self.OT['r'][0]
        Zstart = self.OT['z'][0]
        if plotit:
            ax.plot(Rstart,Zstart,'go')
        s = EFITutils.follow_fieldlines_rzphi_dl(geq.g,Rstart,Zstart,0,-dL,abs(Lmax/dL));
        if plotit:
            ax.plot(s['r'],s['z'],'r-')

        # Find intersection
        ccInt = GeoUtils.int_curve_curve(s['r'][0:s['iLastGood']+1],s['z'][0:s['iLastGood']+1],self.IT['r'],self.IT['z'])
        if plotit: 
            ax.plot(ccInt['pInt'][0],ccInt['pInt'][1],'ro')
            
        # combine
        a = np.flip(self.IT['r'][ccInt['foundInd2'][0]+1:])
        b = np.flip(s['r'][0:ccInt['foundInd1'][0],0])
        c = self.OT['r']
        div = {}
        div['r'] = np.concatenate((a,[ccInt['pInt'][0]],b,c))

        a = np.flip(self.IT['z'][ccInt['foundInd2'][0]+1:])
        b = np.flip(s['z'][0:ccInt['foundInd1'][0],0])
        c = self.OT['z']
        div['z'] = np.concatenate((a,[ccInt['pInt'][1]],b,c))
                
        self.div = div
        
    def make_outer_wall(self,geq,ax=None):
        if ax is None:
            plotit = False
        else:
            plotit = True
            
        Lmax = self.LmaxOW
        dL = self.dL_ff
        print("Following line for outer wall")
        Rstart = self.OT['r'][-1]
        Zstart = self.OT['z'][-1]
        s = EFITutils.follow_fieldlines_rzphi_dl(geq.g,Rstart,Zstart,0,-dL,abs(Lmax/dL));
        if plotit:
            ax.plot(s['r'],s['z'],'r-')
        
        self.OW = {'r':s['r'][:,0],'z':s['z'][:,0]}
        
    def make_inner_wall(self,geq,ax=None):
        if ax is None:
            plotit = False
        else:
            plotit = True
            
        Lmax = self.LmaxIW
        dL = self.dL_ff
        print("Following line for inner wall")
        Rstart = self.IT['r'][-1]
        Zstart = self.IT['z'][-1]
        s = EFITutils.follow_fieldlines_rzphi_dl(geq.g,Rstart,Zstart,0,dL,abs(Lmax/dL));
        if plotit:
            ax.plot(s['r'],s['z'],'r-')
    
        self.IW = {'r':s['r'][:,0],'z':s['z'][:,0]}
            
    def combine_pfcs(self,ax=None):
        print('Combining pfc contours')
        a = np.flip(self.IW['r'][1:])
        b = self.div['r']
        c = self.OW['r'][1:]            
        r = np.concatenate((a,b,c))

        a = np.flip(self.IW['z'][1:])
        b = self.div['z']
        c = self.OW['z'][1:]            
        z = np.concatenate((a,b,c))


        self.pfc = {'r':r,'z':z}

    
# -------------------------------------------------------------------------------------------------------------------------
def make_target_shape_oneway(g,targetAngle,rStart,zStart,length,direction):
    # Angle is returned at midpoints of segments, so size is one smaller than r,z
    # direction = 1 : "forward".  
    # direction = -1: "backward". 
    
    
    if not any((direction == 1,direction == -1)):
        print("Error: direction must be 1 or -1")
        return None
    
    # Todo: allow targetAngle to be array
    
    # Output quantities and resolution
    nPts = 20
    angle = np.zeros((nPts-1,))
    r = np.full((nPts,),np.nan)
    z = np.full((nPts,),np.nan)
    s = np.full((nPts,),np.nan)
    r[0] = rStart
    z[0] = zStart
            
    
    # Make a semicircle around fieldline in RZ plane and search for target angle
    nCircle = 180 # Test resolution
    tol = 0.1  # Sanity check angle tolerance
    iLastGood = None
    for i in range(nPts-1):
        rTest = np.full((nCircle,),np.nan)
        zTest = np.full((nCircle,),np.nan)
        aTest = np.full((nCircle,),np.nan)
        b = EFITutils.calc_bfield(g,r[i],z[i])
        if b['ierr']:
            print('Break on b error: nans will be present')
            break        
        for j in range(nCircle):
            angleOffset = np.arctan2(direction*b['Bz'],direction*b['Br']) + np.pi/2
            rTest[j] = length/nPts*np.cos((j-1)*np.pi/nCircle + angleOffset) + r[i]
            zTest[j] = length/nPts*np.sin((j-1)*np.pi/nCircle + angleOffset) + z[i]
            angInfo = EFITutils.calc_Bangle_g(g,r[i],z[i],rTest[j],zTest[j],npts=1)
            aTest[j] = direction*angInfo['alphaDeg'][0]
                
        minTest = np.abs(aTest - targetAngle)
        indMin = np.argmin(minTest)
        thisMin = minTest[indMin]
        if thisMin > tol:
            print("Error: Failed sanity check. Maybe points went off gfile?")
            return None
        
        r[i+1] = rTest[indMin]
        z[i+1] = zTest[indMin]
        angle[i] = aTest[indMin]
        iLastGood = i + 1
        
    s[1:] = np.cumsum(np.sqrt(np.diff(r)**2 + np.diff(z)**2))
    
    return {'r':r,'z':z,'s':s,'angleDeg':angle,'iLastGood':iLastGood}

        
if __name__=="__main__":
    geq = EFITutils.geq('g000001.00001')
    geq.setup_bfield()

    fig = plt.figure(figsize=(8.5,11.0))
    ax = fig.add_subplot(111, aspect='equal')

    geq.plot(ax)
    this = None
    div = db.DivertorBuild()
    div.make_divertors(geq)
    div.connect_divertors(geq)
    div.make_outer_wall(geq)
    div.make_inner_wall(geq)
    div.combine_pfcs()
    div.plot(ax)