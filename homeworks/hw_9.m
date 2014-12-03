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

% part 1
% vehicle burnout
v_burnout = 14.5; % km/s
alt_burnout = 150; % km
r_burnout = earth.radius_km + alt_burnout; % km

a_hyp = -earth.mu/2*(v_burnout^2/2 - earth.mu/r_burnout)^-1;
e_hyp = 1 - r_burnout/a_hyp;
delta_ = 2*asin(1/e);
v_infin = sqrt(-earth.mu/a_hyp);

%part 2


% initializing vectors
nu_ = linspace(0,2*pi,N);
rp_vec = zeros(1,N);
ra_vec = zeros(1,N);

delta_pi = (pi + delta_)/2;

% ---------- %
% Problem #1 %
% ---------- %

% iterating between angles
for i = 1:1:N
	% calculate adjusted angle
	beta_ = nu_(i) + delta_pi;

	% calculate velocity of s/c relative to the sun
	v_sc = v_infin*[sin(beta_),cos(beta_)]; % [r_hat, theta_hat] wrt sun
	v_net = v_sc + v_earth;
	v_mag = norm(v_net,2);

	% calculate semi major axis of s/c orbit wrt sun
	e_vec = ((v_mag^2-sun.mu/earth.semi_major_axis_km)*earth.semi_major_axis_km*[1 0] - (dot(earth.semi_major_axis_km*[1 0],v_net))*v_net)/sun.mu;
	e = norm(e_vec,2);
	a = (-sun.mu/2*(v_mag^2/2-sun.mu/earth.semi_major_axis_km)^-1);
	
	% solve for ra, rp, and e for the new orbit of the s/c wrt sun
	ra_vec(i) = a*(1+e);
	rp_vec(i) = a*(1-e);

	% grab the orbits that have aphelion at or beyond jupiter
	if(ra_vec(i)>=jupiter.semi_major_axis_km)
		% ---------- %
		% Problem #2 %
		% ---------- %


	end
end

% plot results for problem #1
figure
hold on
title('Part 1')
xlabel('true anomaly (deg)');
ylabel('distance (km)');
plot(nu_*180/pi, ra_vec, "*r;aphelion;");
plot(nu_*180/pi, rp_vec, "-k;perihelion;");