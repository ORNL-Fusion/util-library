function r= rotz( theta)
	%returns rotation matrix
	%rotates coordinates about z axis, x into y
	ct=cos(theta);
	st=sin(theta);
	r=rotcz(ct,st);
