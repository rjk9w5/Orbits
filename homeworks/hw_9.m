% ----------------------- %
% Author: Ryan Krattiger  %
% Assignment: Homework #9 %
% Due Date: 4 Dec 2014    %
% ----------------------- %
clc
clear all
close all
% include custom libraries
addpath ../include/orbiteq

% Load constants from table
load ../Constants/earth.mat
load ../Constants/sun.mat
load ../Constants/jupiter.mat


AU = 1.496e8;
N = 10000; % defining the resolution of the graphs (points b/w 0 & 2pi)
flag1 = false;
flag2 = false;

% Units:
%	angle: radians
%	speed: km/s
%	distance: km
%	vector % [r_hat, theta_hat]

% ---------------- %
% Pre-Calculations %
% ---------------- %

% Earth/Jupiter circular orbit) wrt sun
% letting the orbit radius equal the semi major axis from the constants table
v_earth = vatr(sun.mu,earth.semi_major_axis_km,earth.semi_major_axis_km)*[0 1]; 
v_jupiter = vatr(sun.mu,earth.semi_major_axis_km,earth.semi_major_axis_km)*[0 1];

% part 1
% vehicle burnout
v_burnout = 14.5; % km/s
alt_burnout = 150; % km
r_burnout = earth.radius_km + alt_burnout; % km

a1 = -earth.mu/2*(v_burnout^2/2 - earth.mu/r_burnout)^-1;
e1 = 1 - r_burnout/a1;
delta_ = 2*asin(1/e1);
v_inf_1 = sqrt(-earth.mu/a1);
delta_pi = (pi + delta_)/2;

%part 2
esc_count=0;
jup_count=0;

% initializing vectors
nu_ = linspace(0,2*pi,N);
rp_vec = zeros(1,N);
ra_vec = zeros(1,N);

% ---------- %
% Problem #1 %
% ---------- %

% iterating between angles
for i = 1:1:N
	% calculate adjusted angle
	beta_ = nu_(i) + delta_pi;

	% calculate velocity of s/c relative to the sun
	v_sc = v_inf_1*[sin(beta_),cos(beta_)]; % [r_hat, theta_hat] wrt sun
	v_dp_ = v_sc + v_earth;
	v_dp = norm(v_dp_,2);

	% calculate a & e of s/c orbit wrt sun
	e_vec = eccenvec(sun.mu, v_dp_, earth.semi_major_axis_km*[1 0]);
	e = norm(e_vec,2);
	a = (-sun.mu/2*(v_dp^2/2-sun.mu/earth.semi_major_axis_km)^-1);
	
	% solve for ra, rp, and e for the new orbit of the s/c wrt sun
	ra_vec(i) = a*(1+e);
	rp_vec(i) = a*(1-e);

	% grab the orbits that have aphelion at or beyond jupiter
	if(ra_vec(i)>=jupiter.semi_major_axis_km)
		% ---------- %
		% Problem #2 %
		% ---------- %
		jup_count = jup_count + 1;
		if ~flag1
			flag1 = ~flag1;
			nu_jup_s = i;
		end

		% calculate gamma of v_ar at jupiter
		v_ar = vatr(sun.mu, a, jupiter.semi_major_axis_km);
		h = spmomentum(sun.mu,a,e);
		gamma_ = acos(h/(v_ar*jupiter.semi_major_axis_km));

		% calculate v_inf_2
		v_ar_ = v_ar*[sin(gamma_) cos(gamma_)];
		v_inf_2 = v_ar_ - v_jupiter;
		v_inf_2 = norm(v_inf_2,2);

		% calculate beta_p
		beta_p = (acos((v_ar*cos(gamma_) - norm(v_jupiter,2)) / (v_inf_2)));

		% calculate delta_
		a2 = -jupiter.mu/v_inf_2^2;
		e2 = 1 - 1e6/a2;
		delta_ = 2*asin(1/e2);

		% calculate alpha_
		alpha_ = beta_p - delta_;

		% calculate v_dp_ as a vector
		v_inf_2 = v_inf_2*[sin(alpha_) cos(alpha_)];
		v_dp_ = v_jupiter + v_inf_2;
		v_dp = norm(v_dp_,2);

		% calculate a & e after the flyby
		e_vec = eccenvec(sun.mu, v_dp_, jupiter.semi_major_axis_km*[1 0]);
		e = norm(e_vec,2);
		a = (-sun.mu/2*(v_dp^2/2-sun.mu/jupiter.semi_major_axis_km)^-1);

		% solve aphelion after jupiter flyby
		ra_jup(jup_count) = a*(1+e);
		if e>1
			esc_count = esc_count + 1;
			%display(nu_(i));
			%display('To the final frontier');
			gone(esc_count) = ra_jup(jup_count);
			if ~flag2
				flag2 = ~flag2;
				nu_esc_s = i;
			end
			nu_esc_f = i;
		end % check escape
		
		nu_jup_f = i;
		%display(nu_(i)*180/pi)
	end % part 2
end % main loop

% plot results for problem #1
figure 1
hold on
title('Part 1')
xlabel('true anomaly (deg)');
ylabel('distance (km)');
plot(nu_, ra_vec, ':b*');
plot(nu_, rp_vec, ':r*');

figure 2
hold on
title('Part 2')
xlabel('true anomaly (deg)');
ylabel('distance (km)');
plot(nu_(nu_jup_s:nu_jup_f), ra_jup, ':go')

figure 3
hold on
title('Part 2')
xlabel('true anomaly (deg)');
ylabel('distance (km)');
plot(nu_(nu_esc_s:nu_esc_f), gone, ':r*')


% remove paths
rmpath ../include/orbiteq