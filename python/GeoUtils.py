# -------------------------------------------------------------------------------------------------------------------------
# Collection of basic geometry routines, e.g., moving a distance along a line, finding intersections
#
# JDL
# -------------------------------------------------------------------------------------------------------------------------

import numpy as np 

# -------------------------------------------------------------------------------------------------------------------------
# Find intersection of two curves linearly by line segments.
#
# Curves are defined by array of points. 
# First curve is stepped and both curves linearly interpolated to find intersection with second curve.
# If first == true then checking is stopped when first intersection is found. 
#
# JDL
# -------------------------------------------------------------------------------------------------------------------------
def int_curve_curve(line1R,line1Z,line2R,line2Z,first=True):

    if line1R.size != line1Z.size:
        print('Error: line1 R and Z do not match')
        return {'ierr':True}
    if line2R.size != line2Z.size:
        print('Error: line2 R and Z do not match')
        return {'ierr':True}

    intCount = 0
    foundInd1 = []
    foundInd2 = []
    pInt = []
    
    for i in range(line1R.size - 1):
        p1 = np.array((line1R[i]  ,line1Z[i]  ))
        p2 = np.array((line1R[i+1],line1Z[i+1]))
        this = int_line_curve(p1,p2,line2R,line2Z,first=first)     
        if not this['ierr']:            
            ierr = False
            pInt.append(this['pInt'])
            foundInd1.append(i)
            foundInd2.extend(this['foundInd'])
            intCount += this['intCount']
            
            if first:
                done = True
                break
    
    if intCount == 0:
        ierr = True
        
    pInt = np.squeeze(np.asarray(pInt))
    
    return {'pInt':pInt,'ierr':ierr,'foundInd1':foundInd1,'foundInd2':foundInd2,'intCount':intCount}

    return None

# -------------------------------------------------------------------------------------------------------------------------
# Find intersection of a line segment with a curve.
#
# Curve is defined by array of points linearly interpolated to find intersection with line.
#
# JDL 
# -------------------------------------------------------------------------------------------------------------------------
def int_line_curve(p1,p2,lineR,lineZ,first=True):
    
    if p1.size != 2:
        print('Error: p1 must be a point')   
        return {'ierr':True}
    if p2.size != 2:
        print('Error: p2 must be a point')
        return {'ierr':True}
    
    intCount = 0
    foundInd = []
    pInt = []
    for i in range(lineR.size - 1):
        p3 = np.array((lineR[i],lineZ[i]))
        p4 = np.array((lineR[i+1],lineZ[i+1]))            
        this = int_two_lines(p1,p2,p3,p4)
        
        # Check for parallel line segments
        if this['ierr']: 
            continue
     
        if i == (lineR.size - 1):
            test = this['u1'] >= 0 and this['u1'] <= 1 and this['u2'] >= 0 and this['u2'] <= 1
        else:
            test = this['u1'] >= 0 and this['u1'] <= 1 and this['u2'] >= 0 and this['u2'] < 1
        
        if test:
            intCount += 1                    
            pInt.append(p1 + this['u1']*(p2 - p1) )
            ierr = False
            foundInd.append(i)
            if first == True:
                break    

    if intCount == 0:
        ierr = True
    
    pInt = np.asarray(pInt)
    
    return {'pInt':pInt,'ierr':ierr,'foundInd':foundInd,'intCount':intCount}

# -------------------------------------------------------------------------------------------------------------------------
# Calculates the intersection of lines p1 to p2, and p3 to p4.  Each
# point is an x-y pair.  Returns two values, which are the normalized distances
# to the intersection point along the line p1 to p2 (u1) and the line
# p3 to p4 (u2).
#
# Equations for line 1 (x1 = [x,y]), and line 2 (x2 = [x,y])
#   _   _        _  _
#   x1 = p1 + u1(p2-p1)
#   _   _        _  _
#   x2 = p3 + u2(p4-p3)
#
# Then solve equations for u1, u2
# 
# Note the distances can be greater than 1.0, i.e., these are not line segments!
# 
# JDL
# -------------------------------------------------------------------------------------------------------------------------
def int_two_lines(p1,p2,p3,p4):
    
    tol = 1e-15;
    
    denom = (p4[1]-p3[1])*(p2[0]-p1[0]) - (p4[0]-p3[0])*(p2[1]-p1[1])
    if np.abs(denom) < tol: # parallel lines
        u1 = np.nan
        u2 = np.nan
        pInt = np.nan
        ierr = True
    else:
        u1 = ((p4[0]-p3[0])*(p1[1]-p3[1]) - (p4[1]-p3[1])*(p1[0]-p3[0]))/denom
        u2 = ((p2[0]-p1[0])*(p1[1]-p3[1]) - (p2[1]-p1[1])*(p1[0]-p3[0]))/denom
        pInt = p1 + u1*(p2-p1)
        ierr = False
        
    return {'u1':u1,'u2':u2,'pInt':pInt,'ierr':ierr}
    
# -------------------------------------------------------------------------------------------------------------------------
# Move a distance L along a curve. 
#
# JDL
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