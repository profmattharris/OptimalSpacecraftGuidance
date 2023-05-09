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
% Section 4.7
% MATLAB Implementation of Chebyshev Discretization

clear all; close all; clc

% N and node times
N = 10;
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

% Function evaluation and derivative approximation
f  = sin(t)+.5*t.^2;
df = D*f;

% Error computation
fp = cos(t)+1.0*t;
error = fp-df;

figure, plot(t,error), grid on
xlabel('t'), ylabel('error')