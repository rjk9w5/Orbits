addpath '../include/tridorbits/'
addpath '../include/orbiteq'
addpath '../Constants/'

load earth.mat

% Given

dv_v = [.1 .1 .1]; % km/s 
r0_v = [-3864 751 -5342]; %km 
tof = 6500; % s 
riss_v = [-2965 3443 -5045]; % km

% from black box
a_tran = 7794.51; % km
e_tran = .151;
v1_tran = [-1.7353 -8.1141 .3894]; % km/s 

% velocity of parking orbit (prior to ingnition of the residual fuel)
v0_v = v1_tran - dv_v;
[a,e,i,RAAN,arg_periap,M] = cart2clas([r0_v v0_v],earth.mu)

v2 = vatr(earth.mu, a_tran, norm(riss_v,2));

h_vec = cross(r0_v,v1_tran);
A = [0 -riss_v(3) riss_v(2);riss_v(3) 0 riss_v(1);-riss_v(2) riss_v(1) 0];
B = h_vec';
v2_tran = inv(A)*B
norm(v2_tran)
v2

rmpath '../include/tridorbits'
rmpath '../Constants/'