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
% MATLAB Implementation of Indirect Shooting

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

% Continuation procedure
% s = 1e-5; unk = 1e-3*[1;1;1];
% s = 1e-1; unk = [2.1454783539651; 61.7820129652734; 136.3627108440081];
s = 1; unk = [2.1455405424997; 61.7825876741792; 136.3528303294001];

% Shooting method
ops = optimoptions('fsolve','Display','iter-detailed','MaxFunEvals',500);
sol = fsolve(@shooting,unk,ops,x0,m0,mdot,T,mu,tf,s);

% Simulation and plots
[t,x] = ode45(@ode,[0,tf],[x0;sol],[],m0,mdot,T,mu);
figure, plot(t,x(:,1)), grid on
xlabel('t (days)'), ylabel('r (AU)')

theta = atan2(x(:,5),x(:,6));
theta = 180/pi*wrapTo2Pi(theta);
figure, plot(t,theta), grid on
xlabel('t (days)'), ylabel('theta (deg)')


function xdot = ode(t,x,m0,mdot,T,mu)

% States and costates
r  = x(1); u  = x(2); v  = x(3);
lr = x(4); lu = x(5); lv = x(6);

% Optimal thrust angle
theta = atan2(lu,lv);

% State equations
rdot = u;
udot = v^2/r - mu/r^2 + T*sin(theta)/(m0 - mdot*t);
vdot = -u*v/r + T*cos(theta)/(m0 - mdot*t);

% Costate equations
lrdot = -lu*(-v^2/r + 2*mu/r^3) - lv*(u*v/r^2);
ludot = -lr + lv*v/r;
lvdot = -lu*2*v/r + lv*u/r;

xdot = [rdot; udot; vdot; lrdot; ludot; lvdot];
end

function F = shooting(unk,x0,m0,mdot,T,mu,tf,s)

lambda0 = 1;

% Initial costates
lr = unk(1);
lu = unk(2);
lv = unk(3);

% Integrate state/costate system
[~,x] = ode45(@ode,[0,tf],[x0;lr;lu;lv],[],m0,mdot,T,mu);

% Terminal states and costates
rf  = x(end,1); uf  = x(end,2); vf  = x(end,3);
lrf = x(end,4); luf = x(end,5); lvf = x(end,6);

% Terminal constraints
F(1) = uf;
F(2) = vf - sqrt(mu/rf);
F(3) = ( lrf - lambda0 - lvf*sqrt(mu)/(2*rf^(3/2)) )*s;
end