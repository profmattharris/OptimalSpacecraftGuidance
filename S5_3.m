% Optimal Spacecraft Guidance
% Matthew W. Harris and M. Benjamin Rose
%
% This code is provided for academic/educational use only. 
% Use at your own risk.
%
% If you use the code, please cite the book.
% @book{OptimalSpacecraftGuidance,
%   author    = {Matthew W. Harris and M. Benjamin Rose},
%   title     = {Optimal Spacecraft Guidance},
%   publisher = {Utah State University},
%   address   = {Logan, UT},
%   year      = {2023}
% }
%
% Section 5.3
% MATLAB Implementation

clear all; close all; clc

% Data
h  = 50;    % m
L  = 100;   % m
ve = 25;    % m/s
a  = 10;    % m/s2

% Solve for the final time
tf = roots([-1/4*a^2,0,ve^2,2*L*ve,h^2+L^2]);
tf = tf(1);

% Solve for the thrust angle
theta = atan( h/(L+ve*tf) );

% Final position of the evader
xe = L+ve*tf;
ye = h;

% Final position of the pursuer
xp = 1/2*a*tf^2*cos(theta);
yp = 1/2*a*tf^2*sin(theta);

% Position error
error = norm([xe-xp,ye-yp])