            function p_e = LL2XYZ(mat,alt)
% %===========================================================%
% %           function p_e = wgslla2xyz(lat,lon,alt)          %
% %                                                           %
% %   This function returns the position vector p_e given in  %
% %   the WGS-84 Earth Centered Earth Fixed (ECEF) coordinate %
% %   for a user located at the goedetic coordinates lat,     %
% %   lon and alt.  The units of the output position vector,  %
% %   p_e, are meters while the inputs have the following     %
% %   units: latitude and longitude in degrees and            %
% %   altitude in meters.                                     %
% %                                                           %
% %   Programmer:     Demoz Gebre-Egziabher                   %
% %   Created:        July 2, 1998                            %
% %   Last Modified:  March 26, 2009                          %
% %   License:  BSD  see bsd.txt                              %
% %                                                           %
% %===========================================================%     

%   Load ellipsoid constants

wgs_84_parameters;
deg2rad = pi/180;
%   Compute East-West Radius of curvature at current position
[r,c]=size(mat);
for k=1:c
    for p=1:length(mat)
    R_E = R_0/(sqrt(1 - (e*sin(deg2rad*mat(2,k))^2)));

    %   Compute ECEF coordinates
if(exist('alt','var') == 0)
    alt = zeros(1000);
end
    
    p_e(1,k) = (R_E + alt(p))*cos(deg2rad*mat(2,k)*cos(deg2rad*mat(1,k)));
    p_e(2,k) = (R_E + alt(p))*cos(deg2rad*mat(2,k))*sin(deg2rad*mat(1,k));
    p_e(3,k) = ((1 - e^2)*R_E + alt(p))*sin(deg2rad*mat(2,k));
    end
end
%===========================================================%