clearvars;


% Give two points to define line
p1 = [2.25,-1.0];
p2 = [2.25,-0.8];

% point 
p0 = [2.23,-0.92];



[dist,pu,CUTOFF,pu_unchecked,dist_unchecked] = distance_point_to_line_seg(p1,p2,p0)


% Give two points to define line
p13 = [2.25,-1.0,0.1];
p23 = [2.25,-0.8,-.10];

% point 
p03 = [2.23,-0.92,0];

disp('--------')

[dist,pu,CUTOFF,pu_unchecked,dist_unchecked] = distance_point_to_line_seg(p13,p23,p03)
