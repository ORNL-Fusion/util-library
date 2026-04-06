function fil = define_MPEX_coil_filaments
% nturns = number of turns per coil layer in each coil (see note below)
% nlayers = number of layers in each coil 
% rr1 = inner radius of each coil
% rr2 = outer radius of each coil
% cl = axial length of each coil
% z0 = (axial) starting position of each coil relative to z=0. The coil extends from        
% this location in the direction of larger z
%
% layers in radius, turns in axial
% coils centered around r=0
%
% Windings are constructed downstream in build_circular_coil*. The
% numbering used there is shown schematically below
%  r2 -------------------
%     |     |     |     |
%     |  2  |  4  |  6  |        % two layers
%     |     |     |     |
%     -------------------
%     |     |     |     |
%     |  1  |  3  |  5  |        % three turns
%     |     |     |     |
%  r1 -------------------
%    z1                  z1+dz
%
%
% Latest update from CoilSetup_MPEX.xlsx, dated 2021_04_01 ("Obtained from Earl from the MPEX DAC")
% Table is pasted below and implemented in the code
% JD Lore

data = mpex_coil_table;

fil.ncoils = size(data,2);
fil.nturns  = data(7,:);
fil.nlayers = data(8,:);
fil.nwind = fil.nturns.*fil.nlayers;

fil.z0  = data(3,:)./1000;
fil.cl  = data(4,:)./1000;
fil.rr1 = data(5,:)./1000;
fil.rr2 = data(6,:)./1000;

thick = fil.rr2 - fil.rr1;

fil.area = fil.cl.*thick;


function data = mpex_coil_table
data = [ ...
1	5	-1638.3	76.2	411	484.8	64	65
2	5	-1219.2	76.2	411	484.8	64	65
3	5	-380	64	241.3	358.6	11	6
4	5	-140	64	241.3	358.6	11	6
5	5	140	64	241.3	358.6	11	6
6	5	380	64	241.3	358.6	11	6
7	5	1104.9	76.2	411	484.8	64	65
8	5	1524	76.2	411	484.8	64	65
9	5	2108.2	76.2	411	484.8	64	65
10	5	2527.3	76.2	411	484.8	64	65
11	5	3111.5	76.2	411	484.8	64	65
12	5	3530.6	76.2	411	484.8	64	65
13	5	3810	76.2	411	484.8	64	65
14	5	4165.6	76.2	411	484.8	64	65
15	5	4521.2	76.2	411	484.8	64	65
16	5	4876.8	76.2	411	484.8	64	65
17	5	5283.2	76.2	411	484.8	64	65
18	5	5800.1	348	850.9	885.5	30	289
19	5	6224.3	348	850.9	885.5	30	289
20	5	6648.5	348	850.9	885.5	30	289
21	5	7402.8	348	850.9	885.5	30	289
22	5	7827	348	850.9	885.5	30	289
23	5	8251.2	348	850.9	885.5	30	289
].';


%% Raw table
% coil_index	ps	datum	z	dz	r_inner	r_outer	layers_z	layers_r	setup	dateEffective	comment
% 1	C1	5	-1638.3	76.2	411	484.8	64	65	1	2021_04_01	Obtained from Earl from the MPEX DAC
% 2	C2	5	-1219.2	76.2	411	484.8	64	65	1		
% 3	C3	5	-380	64	241.3	358.6	11	6			
% 4	C4	5	-140	64	241.3	358.6	11	6			
% 5	C5	5	140	64	241.3	358.6	11	6			
% 6	C6	5	380	64	241.3	358.6	11	6			
% 7	C7	5	1104.9	76.2	411	484.8	64	65			
% 8	C8	5	1524	76.2	411	484.8	64	65			
% 9	C9	5	2108.2	76.2	411	484.8	64	65			
% 10	C10	5	2527.3	76.2	411	484.8	64	65			
% 11	C11 	5	3111.5	76.2	411	484.8	64	65			
% 12	C12	5	3530.6	76.2	411	484.8	64	65			
% 13	C13	5	3810	76.2	411	484.8	64	65			
% 14	C14	5	4165.6	76.2	411	484.8	64	65			
% 15	C15	5	4521.2	76.2	411	484.8	64	65			
% 16	C16	5	4876.8	76.2	411	484.8	64	65			
% 17	C17	5	5283.2	76.2	411	484.8	64	65			
% 18	C18	5	5800.1	348	850.9	885.5	30	289			
% 19	C19	5	6224.3	348	850.9	885.5	30	289			
% 20	C20	5	6648.5	348	850.9	885.5	30	289			
% 21	C21	5	7402.8	348	850.9	885.5	30	289			
% 22	C22	5	7827	348	850.9	885.5	30	289			
% 23	C23	5	8251.2	348	850.9	885.5	30	289			

