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
% Section 3.6
% MATLAB Implementations including Direct Transcription

clear all; close all; clc

% Continuous-time CW System
w = 4 / 3600;
A = [0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    3*w^2, 0, 0, 0, 2*w, 0;
    0, 0, 0, -2*w, 0, 0;
    0, 0, -w^2, 0, 0, 0];
B = [zeros(3,3); eye(3,3)];

% Discretization
dt = 1;
sysc = ss(A,B,[],[]);
sysd = c2d(sysc,dt);
A = sysd.A;
B = sysd.B;

% Analytical Implementation
x0 = [0;-1;0;-.1;.1;0];
xF = [0;0;0;0;0;0];
N  = 5*60;
X  = xF - A^N*x0;
C = [];
for i = 0:N-1
    C = [C, A^(N-1-i)*B];
end
U1 = C.'*inv(C*C.')*X; 

% Simulation
u1 = reshape(U1,3,N);
x1(:,1) = x0;
for i = 1:N
    x1(:,i+1) = A*x1(:,i) + B*u1(:,i);
end

figure, plot(x1(2,:),x1(1,:)), grid on
xlabel('Vertical Position (km)'), ylabel('Horizontal Position (km)')

% Pseudoinverse
U2 = pinv(C)*X;

% Quadratic Programming
options = optimoptions('quadprog','Algorithm','trust-region-reflective');
U3 = quadprog(eye(3*N),zeros(3*N,1),[],[],C,X,[],[],[]);

% YALMIP
opts        = sdpsettings;
opts.solver = 'gurobi';

yalmip('clear')
x    = sdpvar(6,N+1);
u    = sdpvar(3,N);

con = [x(:,1) == x0];
obj = 0;
for k = 1:N
    con = [con, x(:,k+1) == A*x(:,k) + B*u(:,k)];
    obj = obj + u(:,k).'*u(:,k);
end
con = [con, x(:,N+1) == xF];
con = [con, -3.5 <= x(1,:) <= 3.5];
con = [con, -3.5 <= x(2,:) <= 3.5];
sol = solvesdp(con,obj,opts);
    
x = double(x);
u = double(u);
U4 = reshape(u,size(U3));

tx = 0:dt:N*dt;
tu = 0:dt:(N-1)*dt;
figure, plot(x(2,:),x(1,:)), grid on
xlabel('Vertical Position (km)'), ylabel('Horizontal Position (km)')