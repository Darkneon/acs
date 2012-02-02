function [roll,pitch,yaw] = A2rpy_321(A)

% function [roll,pitch,yaw] = A2rpy_321(A)
%
% Transforms Attitude matrix into roll pitch yaw angles
% using a 313 Euler rotation
%
% Inputs: A = 3x3 attitude matrix
%
% Outputs: roll angle (radians)
%	   pitch angle (radians)
%	   yaw angle (radians)
%
% Author: Scott Gleason, 2012
% License: GPLv3
%
%  Ref: Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors
%  By James Diebel, Stanford University, 2006
      
% EULER321 rotation

roll = atan2(-1*A(3,2),A(3,3));
pitch = asin(A(3,1));
yaw = atan2(-1*A(2,1),A(1,1));


