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
%   publisher = {Utah State University Libraries},
%   address   = {Logan, UT},
%   year      = {2023}
% }
%
% Section 1.2
% MATLAB Implementation of Q-guidance

clear all; close all; clc

% Data
r0 = [0;0];             % m
v0 = [0;0];             % m/s
rf = [90000; 45000];    % m
umax = 50;              % m/s2
tf = 100;               %s
t  = linspace(0,tf,1e3)';

% Guidance Simulation
[~,x] = ode45(@ode,t,[r0;v0;35;35],[],tf,rf,umax);
figure, plot(x(:,1),x(:,2)), grid on
xlabel('Range (m)')
ylabel('Altitude (m)')

% Control plots
for i = 1:length(t)
    [~,u(:,i),u_cmd(:,i)] = ode(t(i),x(i,:)',tf,rf,umax);
end
figure, subplot(3,1,1), plot(t,u(1,:),t,u_cmd(1,:)), grid on, ylabel('u1 (m/s2)')
subplot(3,1,2), plot(t,u(2,:),t,u_cmd(2,:)), grid on, ylabel('u2 (m/s2)')
subplot(3,1,3), plot(t,sqrt(u(1,:).^2+u(2,:).^2)), grid on, ylabel('unorm (m/s2)')
xlabel('t (s)')

% ode function
function [xdot,u,u_cmd] = ode(t,x,tf,rf,umax)
r = x(1:2);
v = x(3:4);
u = x(5:6);

% Navigation Block
rhat = r+0.05*r*cos(3*t);
vhat = v+0.05*v*sin(2*t);

% Guidance Block
u_cmd = [0;0];
tgo   = tf-t;
if tgo >= 1
    vr = 1/tgo * ( rf - rhat - 1/2*tgo^2*[0;-9.81] );
    vg = vr - vhat;
    if norm(vg) >= 1
        u_cmd = vg/norm(vg) * umax;
    end
end

% Control Block
sigma = 10;
udot  = sigma*u_cmd-sigma*u;

% Dynamics Block
rdot = v;
vdot = u + 0.75*[0;-9.81] + 0.25*sin(5*t);
xdot = [rdot; vdot; udot];
end