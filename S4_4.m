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
% Section 4.4
%   MATLAB Implementation of Relative Orbit Regulation

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

% Compute the feedback gains and simulate...
SN = eye(6);
Q  = 1e-1*eye(6);
R  = 1e+6*eye(3);

S = SN;
for k = N:-1:1
    K(:,:,k) = inv(R+B'*S*B)*B'*S*A;
    S = Q + A'*S*inv(eye(6)+B*inv(R)*B'*S)*A; 
end

x0 = [100;0;0;0;5;0];
x(:,1) = x0;
for k = 1:N
    u(:,k) = -K(:,:,k)*x(:,k);
    x(:,k+1) = A*x(:,k) + B*u(:,k);
end

figure, plot(x(2,:)',x(1,:)'), grid on
xlabel('y (m)'), ylabel('x (m)')

figure, plot(t(1:N),u.'), grid on
xlabel('t (s)'), ylabel('u (m/s2)')