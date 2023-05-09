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
% Section 4.8
% MATLAB Implementation of Direct Transcription

clear all; close all; clc

% Unit scalings
MU = 1 / 1;
DU = 1 / 149.6e9;
TU = 1 / 86400;

% Data
F    = 4        * MU*DU/TU^2;
m0   = 4500     * MU;
mdot = 7e-5     * MU/TU;
mu   = 1.327e20 * DU^3/TU^2;
tf   = 193; % For book

% Initial conditions
r0 = 1;
u0 = 0;
v0 = sqrt(mu/r0);
x0 = [r0;u0;v0];

% Node time and differentiation matrix
N = 50;
k = 0:N;
t = -cos(pi*k/N)';

% Calculate T
for k = 1:N+1
    T(:,k)  = chebyshevT(k-1,t);
end

% Calculate Tdot
Tdot(:,1) = zeros(N+1,1);
for k = 1:N
    Tdot(:,k+1) = k*chebyshevU(k-1,t);
end

% Calculate the differentiation matrix
D = Tdot*inv(T);

% Scale the nodes and differentiation matrix
t = tf/2*(t+1);
D = 2/tf*D;

% Initial guess
rguess  = linspace(1,2,N+1)';
uguess  = linspace(0,0,N+1)';
vguess  = linspace(sqrt(mu/r0),sqrt(mu/2),N+1)';
thguess = linspace(30,300,N+1)' * pi/180;
xguess  = [rguess; uguess; vguess; thguess];

% Optimize
ops = optimoptions('fmincon','Display','iter',...
                   'EnableFeasibilityMode',true,...
                   'SubproblemAlgorithm','cg',...
                   'MaxFunEvals',5e4);
pobj = @(x) obj(x,N);
pcon = @(x) con(x,r0,u0,v0,m0,mdot,F,mu,t,D,N);
x = fmincon(pobj,xguess,[],[],[],[],[],[],pcon,ops);

% Extract the solution
x = reshape(x,N+1,4);
r = x(:,1); u = x(:,2); v = x(:,3); theta = x(:,4);

% Plot the solution
figure, plot(t,r), grid on
xlabel('t (days)'), ylabel('r (AU)')

figure, plot(t,theta*180/pi), grid on
xlabel('t (days)'), ylabel('theta (deg)')

% Objective function
function J = obj(x,N)
x = reshape(x,N+1,4);
J = -x(N+1,1);
end

% Constraint function
function [cin,ceq] = con(x,r0,u0,v0,m0,mdot,F,mu,t,D,N)
x = reshape(x,N+1,4);
r = x(:,1); u = x(:,2); v = x(:,3); theta = x(:,4);
rdot = u;
udot = v.^2./r - mu./r.^2 + F*sin(theta) ./ (m0-mdot*t);
vdot = -u.*v./r + F*cos(theta) ./ (m0-mdot*t);
cin = [];
ceq = [r(1)-r0; u(1)-u0; v(1)-v0];
ceq = [ceq; D*r-rdot; D*u-udot; D*v-vdot];
ceq = [ceq; u(N+1)-0; v(N+1)-sqrt(mu/r(N+1))];
end