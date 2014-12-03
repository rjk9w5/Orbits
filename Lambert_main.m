clc
close all
clear all

% Position Vector to First Point (km)
X1 = [7000; 9000; -5000];

% Position Vector to Second Point (km)
X2 = [-2000; 8500; 0];

% Time of Flight (sec)
TOF = 86400;

% Gravitational Paramter (km^3/s^2)
mu = 398600;

% JJ flag (0 = parabolic/hyperbolic) (1 = ellipse)
JJ = 1;

% Change in iterated value small enough to consider convergence
tol = 1e-14;

% Maximum Interations
kmax = 100;

% Call Lambert Solver. Note: if conv = 0, convergence DID NOT HAPPEN.
% Increase kmax or decrease tol if applicable.

[A,P,V1,V2,conv] = Lambert(X1,X2,TOF,mu,JJ,tol,kmax)