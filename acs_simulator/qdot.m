function [qdot] = qdot(omega,q)

% function [qdot] = qdot(omega,q)
%
% calculated change in quaternian for given angular rates
%
% Inputs: omega = attitude rate
%	  q = current quaternian 
%
% Outputs: quaternian rates, 4 elements
%
% Author: Scott Gleason, 2012
% License: GPLv3
%
% Ref: Shuster, M., "Survey of Attitude Representations," Journal of the Astronautical Sciences,
% Vol. 41, No. 4, Oct.-Dec. 1993. pp. 439-517.14

qdot(1) = 0.5*(omega(1)*q(4) - omega(2)*q(3) + omega(3)*q(2));
qdot(2) = 0.5*(omega(1)*q(3) + omega(2)*q(4) - omega(3)*q(1));
qdot(3) = 0.5*(-1*omega(1)*q(2) + omega(2)*q(1) + omega(3)*q(4));
qdot(4) = 0.5*(-1*omega(1)*q(1) - omega(2)*q(2) - omega(3)*q(3));


