% ----------------------- %
% Author: Ryan Krattiger  %
% Assignment: Homework #9 %
% Due Date: 4 Dec 2014    %
% ----------------------- %
clc
clear all
close all
% Load constants from table
load ~/Documents/Orbits/Constants/earth.mat
load ~/Documents/Orbits/Constants/sun.mat
load ~/Documents/Orbits/Constants/jupiter.mat

AU = 1.496e8;
N = 100; % defining the resolution of the graphs (points b/w 0 & 2pi)

% Units:
%	angle: radians
%	speed: km/s
%	distance: km
%	vector % [r_hat, theta_hat, z_hat]

% ---------------- %
% Pre-Calculations %
% ---------------- %

% Earth/Jupiter circular orbit) wrt sun
% letting the orbit radius equal the semi major axis from the constants table
v_earth = sqrt(sun.mu/earth.semi_major_axis_km)*[0,1,0]; 
v_jupiter = sqrt(sun.mu/jupiter.semi_major_axis_km)*[0,1,0];

r_EwrtS = [earth.semi_major_axis_km,0,0];
r_JwrtS = [jupiter.semi_major_axis_km,0,0];

% part 1
% vehicle burnout
v_burnout = 14.5; % km/s
alt_burnout = 150; % km
rad_burnout = earth.radius_km + alt_burnout; % km

a_hyp = -earth.mu/2*(v_burnout^2/2 - earth.mu/rad_burnout)^-1;
e_hyp = 1 - rad_burnout/a_hyp;
delta_ = 2*asin(1/e);
v_infin = sqrt(-earth.mu/a_hyp);

%part 2


% initializing vectors
true_anom = linspace(0,2*pi,N);
rp_vec = zeros(1,N);
ra_vec = zeros(1,N);

ang_adj = (pi + delta_)/2;

% ---------- %
% Problem #1 %
% ---------- %

% iterating between angles
for i = 1:1:N
	% calculate adjusted angle
	beta_ = true_anom(i) + ang_adj;

	% calculate velocity of s/c relative to the sun
	v_sc = v_infin*[sin(beta_),cos(beta_),0]; % [r_hat, theta_hat, z_hat] wrt sun
	v_net = v_sc + v_earth;
	v_mag = norm(v_net,2);

	% calculate semi major axis of s/c orbit wrt sun
	e_vec = ((v_mag^2-sun.mu/r_EwrtS(1))*r_EwrtS - (dot(r_EwrtS,v_net))*v_net)/sun.mu;
	e = norm(e_vec,2);
	a = (-sun.mu/2*(v_mag^2/2-sun.mu/earth.semi_major_axis_km)^-1);
	
	% solve for ra, rp, and e for the new orbit of the s/c wrt sun
	ra_vec(i) = a*(1+e);
	rp_vec(i) = a*(1-e);

	% grab the orbits that have aphelion at or beyond jupiter
	if(ra_vec(i)>=r_JwrtS(1))
		% ---------- %
		% Problem #2 %
		% ---------- %

		% find velocity of s/c at jupiter wrt sun
		v_sc = sqrt(2*sun.mu*(1/r_JwrtS(1)-1/(2*a)));
		% find angular momentum 
		h = sqrt(a*sun.mu*(1-e^2));
		% find gamma of s/c at jupiter
		% quad check --assume rendezvous as s/c is departing wrt sun
		v_sc_gamma = acos(h/(r_JwrtS(1)*v_sc)); 
		% get v_in as a vector wrt sun
		v_sc = v_sc*[sin(v_sc_gamma),cos(v_sc_gamma),0];
		% find v_infinity coming into jupiter
		v_infin2 = norm(v_sc + v_jupiter,2);
		% find a,e,delta,beta' for fly by of jupiter
		a_temp = v_infin2^-2*jupiter.mu;
		e_temp = 1 - 1e6/a_temp;
		delta_2 = 2*asin(1/e_temp);
		beta_2 = delta_2+v_sc_gamma;
		v_net = v_infin2*[sin(beta_2),cos(beta_2),0] + v_jupiter;
		v_mag = norm(v_net,2);

		e_vec = ((v_mag^2-sun.mu/r_JwrtS(1))*r_JwrtS - (dot(r_JwrtS,v_net))*v_net)/sun.mu;
	    e = norm(e_vec,2)
	    a = (-sun.mu/2*(v_mag^2/2-sun.mu/r_JwrtS(1))^-1);
	    rp = a*(1-e);
	    if e>1
	    	display(num2str(true_anom(i)));
	    end
	end
end

% plot results for problem #1
figure
hold on
title('Part 1')
xlabel('true anomaly (deg)');
ylabel('distance (km)');
plot(true_anom*180/pi, ra_vec, "*r;aphelion;");
plot(true_anom*180/pi, rp_vec, "-k;perihelion;");