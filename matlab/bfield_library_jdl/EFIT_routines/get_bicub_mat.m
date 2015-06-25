function bicub_mat = get_bicub_mat()

bicub_mat = zeros(16,16);
% ;;;;;; function values at corners
bicub_mat(1,1) = 1.0;
bicub_mat(2,[1,2,3,4]) = [1.,1.,1.,1.];
bicub_mat(3,[1,5,9,13]) = [1.,1.,1.,1.];
bicub_mat(4,:) = 1.;

% ;;;;;;; 1st derivatives at corners: x direction
bicub_mat(5,2)=1.0;
bicub_mat(6,[2,3,4]) = [1.0,2.0,3.0];
bicub_mat(7,[2,6,10,14]) = [1.0,1.0,1.0,1.0];
bicub_mat(8,[2,3,4]) = [1.,2.,3.];
bicub_mat(8,[6,7,8]) = [1.,2.,3.];
bicub_mat(8,[10,11,12]) = [1.,2.,3.];
bicub_mat(8,[14,15,16]) = [1.,2.,3.];

% ;;;;;;; 1st derivatives at corners: y direction
bicub_mat(9,5) = 1.0;
bicub_mat(10,5:8) = [1.,1.,1.,1.];
bicub_mat(11,[5,9,13]) = [1.,2.,3.];
bicub_mat(12,[5,9,13]) = [1.,2.,3.];
bicub_mat(12,[6,10,14]) = [1.,2.,3.];
bicub_mat(12,[7,11,15]) = [1.,2.,3.];
bicub_mat(12,[8,12,16]) = [1.,2.,3.];

% ;;;;;;; cross derivatives at corners
bicub_mat(13,6) = 1.;
bicub_mat(14,[6,7,8]) = [1.,2.,3.];
bicub_mat(15,[6,10,14]) = [1.,2.,3.];
bicub_mat(16,[6,10,14]) = [1.,2.,3.];
bicub_mat(16,[7,11,15]) = [2.,4.,6.];
bicub_mat(16,[8,12,16]) =[3.,6.,9.];

% bicub_mat = bicub_mat.';