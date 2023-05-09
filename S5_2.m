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
% Section 5.2
% MATLAB Implementation

clear all; close all; clc

% Data
a  = 1.1;
b  = 1;
x0 = 1;
xf = 0;
tf = 1;

% Compute Lambda
L = b^2/(2*a)*( exp(2*a*tf)-1 );

% Simulation
[t,x] = ode45(@ode,[0,tf],x0,[],a,b,x0,xf,tf,L);
figure, plot(t,x), grid on
xlabel('t'), ylabel('x')

% Function
function xdot = ode(t,x,a,b,x0,xf,tf,L)
u = -b/L*exp( a*(tf-t) )* (exp( a*(tf) )*x0 - xf);
xdot = a*x+b*u;
end