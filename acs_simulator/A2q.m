function [q] = A2q(A)

% function [q] = A2q(A)
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

trC = A(1,1) + A(2,2) + A(3,3);
q(4) = (0.5)*sqrt(1 + trC);
q(1) = 1/(4*q(4))*(A(2,3) - A(3,2));
q(2) = 1/(4*q(4))*(A(3,1) - A(1,3));
q(3) = 1/(4*q(4))*(A(1,2) - A(2,1));
