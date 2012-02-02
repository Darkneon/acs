function [A] = orbitref2body_321(roll,pitch,yaw)

% function [A] = orbitref2body_321(roll,pitch,yaw)
%
% Transforms roll pitch and yaw angles to an attitude matrix
%
% Inputs: roll angle (radians)
%	   pitch angle (radians)
%	   yaw angle (radians)
%
% Output: A = 3x3 attitude matrix
%
% Author: Scott Gleason, 2012
% License: GPLv3
%
%  Ref: Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors
%  By James Diebel, Stanford University, 2006
      
% EULER321 rotation

A(1,1) =  cos(roll)*cos(pitch);
A(2,1) =  -1*sin(roll)*cos(pitch);
A(3,1) =  sin(pitch);
        
A(1,2) = cos(roll)*sin(pitch)*sin(yaw) + sin(roll)*cos(yaw);
A(2,2) = -1*sin(roll)*sin(pitch)*sin(yaw) + cos(roll)*cos(yaw);
A(3,2) = -1*cos(pitch)*sin(yaw);
        
A(1,3) = -1*cos(roll)*sin(pitch)*cos(yaw) + sin(roll)*sin(yaw);
A(2,3) = sin(roll)*sin(pitch)*cos(yaw) + cos(roll)*sin(yaw);
A(3,3) = cos(pitch)*cos(yaw);        


