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
% Section 4.6
% MATLAB Implementation of Descent Trajectory Following

clear all; close all; clc

% Nominal boundary conditions and time
r0 = [1000; 1000];
v0 = [-25; 0];
r1 = [0;0];
v1 = [0;0];

t1 = 100;
N  = 1e3;
t  = linspace(0,t1,N);

% Other constants
g  = [0;-9.81];
I  = eye(2);
Z  = zeros(2);

% Nominal trajectory
C = [t1^2*I, t1^3*I; 2*t1*I, 3*t1^2*I] \ [r1-r0-v0*t1; v1-v0];
c2 = C(1:2); c3 = C(3:4);
u_nom = 2*c2 + 6*c3*t-g;
v_nom = v0 + 2*c2*t + 3*c3*t.^2;
r_nom = r0 + v0*t + c2*t.^2 + c3*t.^3;
x_nom = [r_nom; v_nom];

% Discretization
A = [Z, I; Z Z];
B = [Z; I];
sysc = ss(A,B,[],[]);
sysd = c2d(sysc,t(2)-t(1));
A = sysd.A;
B = sysd.B;

% Compute the feedback gain
SN = eye(4);
Q  = eye(4);
R  = 1e3*eye(2);

S = SN;
for k = N:-1:1
    K(:,:,k) = inv(R+B'*S*B)*B'*S*A;
    S = Q + A'*S*inv(eye(4)+B*inv(R)*B'*S)*A; 
end

% Simulation
dx0 = [50; 75; -20; -5];
x(:,1) = x_nom(:,1)+dx0;
for k = 1:N-1
    dx(:,k)  = x(:,k) - x_nom(:,k);
    du(:,k)  = -K(:,:,1)*dx(:,k);
    u(:,k)   = u_nom(:,k)+du(:,k);
    x(:,k+1) = A*x(:,k) + B*u(:,k) + B*g;
end
figure, plot(x(1,:),x(2,:),x_nom(1,:),x_nom(2,:)), grid on
xlabel('Range (m)'), ylabel('Altitude (m)')
figure, plot(t(1:N-1)',u.',t',u_nom.'), grid on
xlabel('Time (s)'), ylabel('Control (m/s2)')
figure, plot(t(1:N-1)',dx.'), grid on
xlabel('Time (s)'), ylabel('State Deviation (m)')
figure, plot(t(1:N-1)',du.'), grid on
xlabel('Time (s)'), ylabel('Control Deviation (m)')