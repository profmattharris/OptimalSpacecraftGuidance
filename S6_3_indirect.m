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
% Section 6.3
% MATLAB Implementation of Indirect Shooting

clear all; close all; clc

% Data
r0 = [2000; 100; 1500]; % m
v0 = [-100; -10; 0];    % m/s
m0 = 2100;              % kg
x0 = [r0; v0; m0];  
g  = [0; 0; -3.71];     % m/s2
Isp  = 225;             % s
rho1 = 13e3;            % N
rho2 = 20e3;            % N
alpha = 1/(9.81*Isp);

% Shooting method
ops = optimoptions('fsolve','Display','iter-detailed');
Lguess = [-0.0167; 0.0002; -0.0503; -0.0537; 0.0149; -0.3644; 0.9011];
tguess = 35;
guess  = [Lguess; tguess];
sol = fsolve(@shooting,guess,ops,x0,g,rho1,rho2,alpha);

% Solution Extraction
L = sol(1:7); tf = sol(8);
[t,x] = ode45(@ode,[0,tf],[x0;L],[],g,rho1,rho2,alpha);
r  = x(:,1:3);  v  = x(:,4:6);   m  = x(:,7);
Lr = x(:,8:10); Lv = x(:,11:13); Lm = x(:,14);
for i = 1:length(t)
    [T(i,:),Tmag(i)] = getControl(Lv(i,:),Lm(i,:),m(i),rho1,rho2,alpha);
end

% Plots
figure, plot(t,r,'linewidth',2), grid on
xlabel('Time (s)'), ylabel('Positions (m)')

figure, plot(t,Tmag/1000,'linewidth',2), grid on
xlabel('Time (s)'), ylabel('Thrust magnitude (kN)')



function xdot = ode(~,x,g,rho1,rho2,alpha)

% State/costate extraction
r  = x(1:3);  v  = x(4:6);   m  = x(7);
Lr = x(8:10); Lv = x(11:13); Lm = x(14);

% Optimal control
T = getControl(Lv,Lm,m,rho1,rho2,alpha);

% State equations
rdot = v;
vdot = g + T/m;
mdot = -alpha*norm(T);

% Costate equations
Lrdot = zeros(3,1);
Lvdot = -Lr;
Lmdot = Lv.'*T/m^2;

xdot = [rdot; vdot; mdot; Lrdot; Lvdot; Lmdot];

end

function [T,Tmag] = getControl(Lv,Lm,m,rho1,rho2,alpha)
That = Lv/norm(Lv);
f = norm(Lv) - m*Lm*alpha;
if f >= 0
    Tmag = rho2;
else
    Tmag = rho1;
end
T = That*Tmag;
end

function F = shooting(unk,x0,g,rho1,rho2,alpha)

L = unk(1:7); tf = unk(8);
[~,x] = ode45(@ode,[0,tf],[x0;L],[],g,rho1,rho2,alpha);

rf  = x(end,1:3).';  vf  = x(end,4:6).';   mf  = x(end,7);
Lrf = x(end,8:10).'; Lvf = x(end,11:13).'; Lmf = x(end,14);
Tf  = getControl(Lvf,Lmf,mf,rho1,rho2,alpha);
Hf  = Lrf.'*vf + Lvf.'*(g+Tf/mf) - Lmf*alpha*norm(Tf);
F   = [rf; vf; Lmf-1; Hf];
end