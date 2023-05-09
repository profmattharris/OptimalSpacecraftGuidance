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
% Section 4.5
%   MATLAB Implementation of Relative Orbit Tracking

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

% Create a reference trajectory to follow...
w  = 2*pi/300;
rx = 100*cos(w*t);
ry = 100*sin(w*t);
r = [rx; ry; 0*t];

% Define the data for the LQT problem...
C = [eye(3) zeros(3)];
P = eye(3);
Q = 1e2*eye(3);
R = eye(3);

% Compute the feedback and feedforward gains...
S = C'*P*C;
V = C'*P*r(:,N+1);

Ss(:,:,N+1) = S;
Vs(:,:,N+1) = V;
for k = N:-1:1
    K(:,:,k)  = inv(R+B'*S*B)*B'*S*A;
    Kv(:,:,k) = inv(R+B'*S*B)*B';
    V = C'*Q*r(:,k) + A'*V - A'*S*inv(eye(6)+B*inv(R)*B'*S)*B*inv(R)*B'*V;
    S = C'*Q*C + A'*S*inv(eye(6)+B*inv(R)*B'*S)*A;
    Ss(:,:,k) = S; 
    Vs(:,:,k) = V;
end

% Implement the control and simulate the system...
x0 = [100;0;0;0;5;0];
x(:,1) = x0;
for k = 1:N
    u(:,k)   = -K(:,:,k)*x(:,k) + Kv(:,:,k)*Vs(:,:,k+1);
    x(:,k+1) = A*x(:,k) + B*u(:,k);
end
figure, plot(x(2,:),x(1,:)), grid on
xlabel('y (m)'), ylabel('x (m)')

figure, plot(t(1:N),u.'), grid on
xlabel('t (s)'), ylabel('u (m/s2)')


% Implement the constant feedback control and simulate the system...
% K  = inv(B'*Ss(:,:,1)*B+R)*B'*Ss(:,:,1)*A;
% Kv = inv(B'*Ss(:,:,1)*B+R)*B';
% 
% y(:,1) = x0;
% for k = 1:N
%     v(:,k)   = -K*y(:,k) + Kv*Vs(:,:,k+1);
%     y(:,k+1) = A*y(:,k) + B*v(:,k);
% end
% plot(y(1,:),y(2,:),'--','linewidth',2)
% legend('Reference','Optimal','Suboptimal')

% Create animated plots for fun...
colorMap = [0,0,1; 1 0 0];
xylim = [-100 100 -100 100];
t_now = linspace(t(1),t(N+1),100);
figure
for i = 1:length(t_now)
    rx_now = interp1(t,rx,t_now(i));
    ry_now = interp1(t,ry,t_now(i));
    x_now  = interp1(t,x(1,:),t_now(i));
    y_now  = interp1(t,x(2,:),t_now(i));
    scatter([ry_now; y_now], [rx_now; x_now], [200;100], colorMap, 'filled'), grid on
    axis(xylim)
    title('Blue = Reference, Red = Tracker')
    xlabel('y (m)'), ylabel('x (m)')
    pause(.1)
end