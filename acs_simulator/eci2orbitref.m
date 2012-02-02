function [a_ref] = eci2orbitref(r_eci,v_eci,a_eci);

% function [r_ref,v_ref] = eci2orbitref(r_eci,v_eci);
%
% Transforms vector in ECI to orbit reference frame
%
% Inputs: r_eci = satellite position vector in ECI
% 	  v_eci = satellite velocity vector in ECI
% 	  a_eci = vector to be transformed in ECI
%
% Outputs: a_ref = vector a in the orbit reference frame
%
% Orbit reference frame defined as:
% X axis = Z x Y 
% Y axis = anti angular momentum vector (-satpos x satvel)
% Z axis = opposite spacecraft position 
%
% Author: Scott Gleason, 2012
% License: GPLv3
%

% only roughly tested. No error checking.

% Build transformation matrix from eci to orbit reference frame
r_mag = norm(r_eci);

% 3rd column of trans matrix, opposite position unit vector (towards Earth)
trans_eci2ref(:,3) =  (-1.*r_eci)./r_mag;

% 2nd column of trans matrix, h = anti-angular momentum unit vector
anti_h = cross(v_eci,r_eci);
anti_h_mag = norm(anti_h);
trans_eci2ref(:,2) = anti_h./anti_h_mag; 

% 1st column of trans matrix, r_unit x h_unit
r_cross_h = cross(trans_eci2ref(:,3),trans_eci2ref(:,2));
r_cross_h_mag = norm(r_cross_h);
trans_eci2ref(:,1) = r_cross_h./r_cross_h_mag; 

% Transform a vector to orbit ref frame
a_ref = a_eci*trans_eci2ref;

trans_eci2ref;
