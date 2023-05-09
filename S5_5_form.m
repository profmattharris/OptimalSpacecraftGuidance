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
% Section 5.5
% MATLAB Implementation of Control Function Fitting

clear all; close all; clc

% Unit scalings
MU = 1 / 1;
DU = 1 / 149.6e9;
TU = 1 / 86400;

% Data
T    = 4        * MU*DU/TU^2;
m0   = 4500     * MU;
mdot = 7e-5     * MU/TU;
mu   = 1.327e20 * DU^3/TU^2;
tf   = 193;

% Initial conditions
r0 = 1;
u0 = 0;
v0 = sqrt(mu/r0);
x0 = [r0;u0;v0];

% Control guess
N = 50;
tspan = linspace(0,tf,N);
coeff = [150;100;30;170];

% Optimize
ops = optimoptions('fmincon','Display','iter',...
                   'EnableFeasibilityMode',true,...
                   'SubproblemAlgorithm','cg',...
                   'MaxFunEvals',5e3);
                   
LB = [100;50;5;100];
UB = [200;200;50;200];
pobj = @(coeff) obj(coeff,x0,m0,mdot,T,mu,tspan);
pcon = @(coeff) con(coeff,x0,m0,mdot,T,mu,tspan);
coeff = fmincon(pobj,coeff,[],[],[],[],LB,UB,pcon,ops);
theta = coeff(1)*tanh( (tspan-coeff(2))/coeff(3) )+coeff(4);
theta = theta*pi/180;

% Integrate the solution
[~,x] = ode45(@ode,tspan,x0,[],m0,mdot,T,mu,tspan,coeff);

% Plot the solution
figure, plot(tspan,x(:,1)), grid on
xlabel('t (days)'), ylabel('r (AU)')

figure, plot(tspan,theta*180/pi), grid on
xlabel('t (days)'), ylabel('theta (deg)')

% State differential equation
function xdot = ode(t,x,m0,mdot,T,mu,tvec,coeff)
r = x(1); u = x(2); v = x(3);
theta = coeff(1)*tanh( (t-coeff(2))/coeff(3) )+coeff(4);
theta = theta*pi/180;
rdot = u;
udot = v^2/r - mu/r^2 + T*sin(theta)/(m0-mdot*t);
vdot = -u*v/r + T*cos(theta)/(m0-mdot*t);
xdot = [rdot;udot;vdot];
end

% Objective function
function J = obj(coeff,x0,m0,mdot,T,mu,tspan)
[~,x] = ode45(@ode,tspan,x0,[],m0,mdot,T,mu,tspan,coeff);
J = -x(end,1);
end

% Constraint function
function [cin,ceq] = con(coeff,x0,m0,mdot,T,mu,tspan)
[~,x] = ode45(@ode,tspan,x0,[],m0,mdot,T,mu,tspan,coeff);
rf = x(end,1);
uf = x(end,2);
vf = x(end,3);
cin = [];
ceq(1,1) = uf - 0;
ceq(2,1) = vf - sqrt(mu/rf);
end