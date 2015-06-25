function r= rotcz( ct, st)
  %function r= rotcz( ct, st)
  %angle is given in terms of cos(theta), sin(theta)
  %returns rotation matrix
  %rotates coordinates about z axis, x into y
  %for positive theta this rotates the x and y axes counter
  %clockwise, or the point clockwise
       
  denom=sqrt(ct^2+st^2);
  if denom < 1.e-300
    error('rotation inputs too small');
  else
    ct=ct/denom;
    st=st/denom;
  end
  r=[ ct,st,0; -st,ct,0;0,0,1];

