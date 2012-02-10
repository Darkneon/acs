function [AziEl range] = XYZ2AziEl(XYZ)
% xyz2rae - Transforms cartesian coordinates (XYZ = East,North,Up) to radar coords R,Az,El
%   input  format #1: [...] = xyz2rae(x, y, z)
%   input  format #2: [...] = xyz2rae(xyz)
%   output format #1: [range, azimuth, elevation] = xyz2rae(...)
%   output format #2: rae = xyz2rae(...)
%
%   Inputs:
%     Format #1:
%       x - see definition in help for cart2sph; positive x = East
%       y - see definition in help for cart2sph; positive y = North
%       z - see definition in help for cart2sph; positive z = Up
%
%     Format #2:
%       xyz - same as format #1, with all 3 values/vector bunched together in a vector/matrix
%
%   Outputs:
%     Format #1:
%       range     - sqrt(x^2 + y^2 + z^2)
%       azimuth   - radian angle clockwize from north (= positive y axis)
%       elevation - radian angle from xy plane to positive z axis
%
%     Format #2:
%       rae - same as format #1, with all 3 values/vector bunched together in a vector/matrix
%
%   Example:
%     [range,az,el] = xyz2rae(1,1,1) => range=1.732, az=-0.785, el=0.615
%     rae = xyz2rae([1,1,1])         => rae = [1.732, -0.785, 0.615]
%
%   Notes:
%     Note the different definitions of azimuth here vs. Malab's cart2sph.
%     Also note the different format of input and output args: The input
%     coordinates here may be either singular values or a vector of
%     coordinate points.
%
%     Use the corresponding rae2xyz function for the reverse transformation.
%
%     xyz2rae does NOT take into account earth curvature, Ionosphere beam
%     curving etc. - this simple function uses a simple flat-earth free-space
%     model.
%
%   See also: rae2xyz, cart2sph, cart2pol

% Programmed by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.4 $  $Date: 2007/08/24 10:00:39 $

  %try
      % Process input args
      
      disp('size of XYZ')
      disp(size(XYZ));
      
      if (nargin == 1)
          
          for t=1:length(XYZ),
            x = XYZ(1,t);
            y = XYZ(2,t);
            z = XYZ(3,t);
            
            % Convert using Matlab's generic cart2sph
            [AziEl(1,t),AziEl(2,t),range(t)] = cart2sph(x,y,z);
            
          end
      else
          disp('invalid arguments');
      end
      
      disp('size of x')
      disp(size(z));
      disp('size of y')
      disp(size(y));
      disp('size of z')
      disp(size(z));

  end