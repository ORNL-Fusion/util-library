function r= rotcx( ct, st)
  %function r= rotcx( ct, st)
  %angle is given in terms of cos(theta), sin(theta)
  %returns rotation matrix
  %rotates coordinates about x axis, y into z
  r = circshift(rotcz(ct,st),[1,1]);
