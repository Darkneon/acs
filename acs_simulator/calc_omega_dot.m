function [omega_dot] = omega_dot(I,omega,tau)

% function [omega_dot] = omega_dot(I,omega,tau)
%
% Calculated the angular acceleration due to torques on the body
% reaction wheels not included!
%
% Inputs: I = Inertia Matrix, 3x3
%	  omega = angular velocity vector, 3x1
%	  tau = torque vector, 3x1
%
% Outputs: omega_dot, vector of angular acceleration
%
% Author: Scott Gleason, 2012
% License: GPLv3
%

% total angular momentum of the spacecraft: L = I*w
L = I*omega;

%  omega_dot = tau/I
omega_dot = tau*inv(I);




