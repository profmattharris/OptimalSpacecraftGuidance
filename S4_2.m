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
% Section 4.2
% MATLAB Implementation of Scalar Linear Quadratic Control

clear all; close all; clc

% Data
a  = 1.1;
b  = 1;
x0 = 1;
xN = 0;
N  = 10;

% Optimal Control Fixed Final State
L = b^2*( 1-a^(2*N) ) / (1-a^2);
for k = 0:N-1
    u(k+1) = -b/L * (a^N*x0-xN)*a^(N-k-1);
end

% Simulation
x(1) = x0;
for k = 1:N
    x(k+1) = a*x(k) + b*u(k);
end
figure, plot(0:N,x), grid on
xlabel('k'), ylabel('x')