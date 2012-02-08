function XYZ = AziEl2XYZ(AziEl, r)
% rae2xyz - Transforms radar coordinates R,Az,El to cartesian (XYZ = East,North,Up)
%   input  format #1: [...] = rae2xyz(range, azimuth, elevation)
%   input  format #2: [...] = rae2xyz(rae)
%   output format #1: [east, north, up] = rae2xyz(...)
%   output format #2: xyz = rae2xyz(...)
%
%   Inputs:
%     Format #1:
%       range     - sqrt(x^2 + y^2 + z^2)
%       azimuth   - radian angle clockwize from north (= positive y axis)
%       elevation - radian angle from xy plane to positive z axis
%
%     Format #2:
%       rae - same as format #1, with all 3 values/vector bunched together in a vector/matrix
%
%   Outputs:
%     Format #1:
%       x - see definition in help for cart2sph; positive x = East
%       y - see definition in help for cart2sph; positive y = North
%       z - see definition in help for cart2sph; positive z = Up
%
%     Format #2:
%       xyz - same as format #1, with all 3 values/vector bunched together in a vector/matrix
%
%   Example:
%     [east,north,up] = rae2xyz(1,1,1) => east=0.455, north=0.292, up=0.841
%     xyz = rae2xyz([1,1,1])           => xyz = [0.455, 0.292, 0.841]
%
%   Notes:
%     Note the different definitions of azimuth here vs. Malab's sph2cart.
%     Also note the different format of input and output args: The input
%     coordinates here may be either singular values or a vector of
%     coordinate points.
%
%     Use the corresponding xyz2rae function for the reverse transformation.
%
%     rae2xyz does NOT take into account earth curvature, Ionosphere beam
%     curving etc. - this simple function uses a simple flat-earth free-space
%     model.
%
%   See also: xyz2rae, sph2cart, pol2cart

% Programmed by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.4 $  $Date: 2007/08/24 10:00:39 $

      % Process input args
      if nargin == 2
          azimuth = AziEl(1,:)
          elevation = AziEl(2,:);
          range = r;
      else
          disp('invalid arguments');
      end

      % Transform the azimuth and convert using Matlab's generic sph2cart
      % Note: use -a instead of pi/2-a if azimuth 0 is Eastward, not Northward
      for t=1:length(azimuth),
        XYZ(t) = sph2cart(pi/2-azimuth(t), elevation(t), range(t));
      end
      
  end