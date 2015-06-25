function [poss,vecrot,isign] = symmetrize(pos,nsymm,stellsymm)
%function [poss,vecrot,isign] = symmetrize(pos,nsymm,stellsymm)
%this function takes an input position in r,phi,z space
%which can be anywhere in the stellarator and returns
%a position poss which is equivalent but lies in the
%first period (in the case of no "stellarator symmetry"
%or in the first half field period if there is 
%"stellarator symmetry"
%inputs:
%pos:       [r,phi,z] anywhere in space
%nsymm:     toroidal symmetry number (number of field periods)
%stellsymm: 1 if "stellarator symmetry" (odd half period symmetry) exists
%           0 if "stellarator symmetry" does not exist
%           this input is optional, 0 is the default
%returns:
%
%poss:
%  if stellsym:
%             [r,phi,z] with 0<phi<pi/nsymm and with
%             z appropriately flipped if pi/nsymm<(new)phi<2pi/nsymm
%  else:
%             [r,phi,z] with 0<phi<2*pi/nsymm
%
%vecrot:      a matrix which will rotate
%             any vectors calculated at point poss to their
%             proper direction at the original point "pos".
%             Thus the new vector is isign*vecrot*b
%             gradb should be transformed as isign*vecrot*gradb*vecrot'
%
%isign:       a +1 or -1 which must multiply vectors or gradb matrices
%             calculated at point pos in order to be correct at
%             the point poss
  if nargin < 3
    stellsymm=0;
  end
  r=pos(1);
  phi=pos(2);
  z=pos(3);
  phi = mod(phi,2*pi/nsymm);
  
  if phi >= pi/nsymm & stellsymm
    phi = 2.*pi/nsymm - phi;
    z = -z;
    if nargout > 1
      isign=-1;
      rot = rotcx(-1,0); %a flip about x axis, reverse the vector sign
      addzrot=1;
    end
  else
    if nargout > 1
      isign=1;
      rot=eye(3);
      addzrot=0;
    end
  end
  poss=[r,phi,z];
  if nargout > 1
    phirot = (addzrot+floor(pos(2)*nsymm/pi/2))*2*pi/nsymm;
    vecrot = rotz(-phirot)*rot;
  end

