addpath '../include/tridorbits/'
addpath '../include/orbiteq'
addpath '../Constants/'

load earth.mat

% Part 1
% Given
display('Part 1')
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
[a,e,i,RAAN,arg_periap,~] = cart2clas([r0_v v0_v],earth.mu,1);
a
e
i 

% part d
v2 = vatr(earth.mu, a_tran, norm(riss_v,2));
h_vec = cross(r0_v,v1_tran);
h = norm(h_vec,2);
riss = norm(riss_v,2);

gam = acos(h_tran/(v2*riss));
if(gam<0)
	gam = -gam
end
v2_polar = v2*[sin(gam) cos(gam) 0];
[a,e,i,RAAN,arg_p,E] = cart2clas([r0_v v1_tran],earth.mu,0);
nu = eccen2true(E,e);
v2_cart = dcm(v2_polar',i, RAAN, arg_p,nu)

% Part 2
% Given
display('Part 2');
a = 4*earth.radius_km;
e = .2;
i = 20*pi/180;
RAAN = 200*pi/180;
arg_p = 300*pi/180;
nu = 100*pi/180;
desc_nu = pi - arg_p;
r = ratnu(desc_nu, a, e)
r_cart = dcm([r;0;0],i,RAAN,arg_p,desc_nu)

v_desc = vatr(earth.mu,a,r)
h = sqrt(earth.mu*a*(1-e^2));
gam = acos(h/(v_desc*r));
gam = -gam; % quad check
v0_v = v_desc*[sin(gam) cos(gam) 0];
v0_cart = dcm(v0_v',i,RAAN,arg_p,desc_nu)

v_circ = vatr(earth.mu,r,r)
vn_cart = v_circ*[cos(RAAN) sin(RAAN) 0]

dv_v = vn_cart - v0_cart'


rmpath '../include/tridorbits'
rmpath '../Constants/'