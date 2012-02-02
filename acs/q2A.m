function [A] = q2A(q)

% function [A] = q2A(q)
%
% Transforms Attitude matrix into quaternian
%
% Inputs: A = 3x3 attitude matrix
%
% Outputs: 4 element quaternian
%
% Author: Scott Gleason, 2012
% License: GPLv3
%
% Ref: Shuster, M., "Survey of Attitude Representations," Journal of the Astronautical Sciences,
% Vol. 41, No. 4, Oct.-Dec. 1993. pp. 439-517.14

A(1,1) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
A(2,2) = -1*q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2;
A(3,3) = -1*q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2;

A(1,2) = 2*(q(1)*q(2) + q(4)*q(3));
A(1,3) = 2*(q(1)*q(3) - q(4)*q(2));
A(2,1) = 2*(q(1)*q(2) - q(4)*q(3));
A(2,3) = 2*(q(2)*q(3) + q(4)*q(1));
A(3,1) = 2*(q(3)*q(1) + q(4)*q(2));
A(3,2) = 2*(q(3)*q(2) - q(4)*q(1));


