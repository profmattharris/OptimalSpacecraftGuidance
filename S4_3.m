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
% Section 4.3
%   MATLAB Implementation of Relative Orbit Control

clear all; close all; clc

% Continuous-time CW matrices
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
t  = 0:dt:600;
N  = length(t)-1;
sysc = ss(A,B,[],[]);
sysd = c2d(sysc,dt);
A = sysd.A;
B = sysd.B;

% Compute Lambda and simulate...
xN = [0;0;0;0;0;0];
R = eye(3);
L = 0;
for i = 0:N-1
    L = L + A^(N-i-1)*B*inv(R)*B.'*A.'^(N-i-1);
end

x0 = [100;0;0;0;5;0];
x(:,1) = x0;
for k = 1:N
    u(:,k) = -inv(R)*B.'*A.'^(N-k)*inv(L)*(A^N*x0-xN);
    x(:,k+1) = A*x(:,k) + B*u(:,k);
end

figure, plot(x(2,:),x(1,:)), grid on
xlabel('y (m)'), ylabel('x (m)')

figure, plot(t(1:N),u.'), grid on
xlabel('t (s)'), ylabel('u (m/s2)')