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
% Section 2.3
% MATLAB Implementation of LQR Guidance

clear all; close all; clc


%------%
% Data %
%------%
x0 = [0;-10;0;0;0;0];


%-------------------------%
% Loop with 'pilot' input %
%-------------------------%
X = [];
T = [];
x(1,:) = x0.';
fprintf('\nHuman Pilot Scenario \n\n')
figure
for i = 1:10
    fprintf('Current position is [%4.2f,%4.2f] m. ',x(end,1),x(end,2));
    fprintf('Time is %4.2f s. ',i-1); 
    u = input('u = '); % Input as column [ux;uy;uz]
    [t,x] = ode45(@ode,[i-1,i],x(end,:),[],u);
    X = [X; x];
    T = [T; t];
    for j = 1:length(t)
        plot(x(j,2),x(j,1),'ko','MarkerSize',5), hold on, grid on, axis([-20 20 -20 20])
        pause(0.01)
    end
end


keyboard
%-----------------------%
% Loop with 'autopilot' %
%-----------------------%
clear x
X = [];
T = [];
x(1,:) = x0.';
fprintf('\n\nAuto Pilot Scenario \n\n')
figure
for i = 1:10
    u = getControl( x(end,:).' );
    fprintf('Current position is [%4.2f,%4.2f] m. ',x(end,1),x(end,2));
    fprintf('Time is %4.2f s. ',i-1); 
    fprintf('u = [%4.2f,%4.2f,%4.2f].\n',u(1),u(2),u(3)); pause(2)
    [t,x] = ode45(@ode,[i-1,i],x(end,:),[],u);
    X = [X; x];
    T = [T; t];
    for j = 1:length(t)
        plot(x(j,2),x(j,1),'ko','MarkerSize',5), hold on, grid on, axis([-20 20 -20 20])
        pause(0.01)
    end
end




function xdot = ode(t,x,u)

w = 4; % rad/hr
A = [0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    3*w^2, 0, 0, 0, 2*w, 0;
    0, 0, 0, -2*w, 0, 0;
    0, 0, -w^2, 0, 0, 0];
B = [zeros(3,3); eye(3,3)];
xdot = A*x+B*u;

end


function u = getControl(x)
w = 4; % rad/hr
A = [0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    3*w^2, 0, 0, 0, 2*w, 0;
    0, 0, 0, -2*w, 0, 0;
    0, 0, -w^2, 0, 0, 0];
B = [zeros(3,3); eye(3,3)];
Q = diag([10;1;1;10;1;1]);
K = lqrd(A,B,Q,eye(3),1);
u = -K*x;
end