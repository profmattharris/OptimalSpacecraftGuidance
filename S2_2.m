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
% Section 2.2
% MATLAB Implementation of Polynomial Guidance

clear all; close all; clc

% Data
r0 = [1000; 1000];      % m
v0 = [-25; 0];          % m/s
rf = [0;0];             % m
vf = [0;0];             % m/s
tf = 100;               % s
t  = linspace(0,tf,1e3)';

% Guidance Simulation
[~,x] = ode45(@ode,t,[r0;v0;0;0],[],rf,vf,tf);

figure, plot(x(:,1),x(:,2)), grid on
xlabel('Range (m)'), ylabel('Altitude (m)')

% Control plots
for i = 1:length(t)
    [~,u(:,i),u_cmd(:,i)] = ode(t(i),x(i,:)',rf,vf,tf);
end
figure, plot(t,u.',t,u_cmd.','--'), grid on
xlabel('t (s)'), ylabel('Control (m/s2)')

% ode function
function [xdot,u,u_cmd] = ode(t,x,rf,vf,tf)
r = x(1:2);
v = x(3:4);
u = x(5:6);
g = [0;-9.81];

% Navigation Block
rhat = r + .02*r*sin(2*t);
vhat = v + .02*v*cos(3*t);

% Guidance Block
I = eye(2); tgo = tf-t;
if tgo >= .1
    C = [tgo^2*I, tgo^3*I; 2*tgo*I, 3*tgo^2*I] \ ...
        [rf-rhat-vhat*tgo; vf-vhat];
    c2 = C(1:2); c3 = C(3:4);
    u_cmd = 2*c2 + 6*c3*(t-t)-g;
else
    u_cmd = [0;0];
end

% Control Block
sigma = 10;
udot = sigma*u_cmd - sigma*u;


% Dynamics Block
rdot = v;
vdot = u + g - 1/2*1*1.5*10*norm(v)*v/1000;
xdot = [rdot; vdot; udot];
end