%  getRange	Determine S/C range from ground station
%		using G/S elevation and S/C radius
%
%	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 29 Dec 02
%
%  [Range, ANGL] = SCrange(Elevation, Rsc)
%	Elevation : Elevation of object S/C  [rad]
%	Rsc :	distance array of object S/C from earth center [km]
%	Range :	distance array of object S/C from G/S [km]
%	Theta : orbit portion for G/S contact [rad]
%
%  See also: trackSC, identSC, getDoppler

function [Range, Theta] = getRange(Elevation, Rsc)

global Re
if isempty(Re)
   Re = 6378.13619;		%[km]
end

delta = pi/2 + Elevation;				%[rad] opposite angle from Rsc
Range = zeros(length(Elevation), length(Rsc));		%[km]
ANGL = zeros(length(Elevation), length(Rsc));		%[rad]
for i=1:length(Elevation)
   Range(i,:) = Re*cos(delta(i)) + sqrt((Re*cos(delta(i)))^2 + (Rsc.^2 - Re^2));    %[km]
   lambda = asin( sin(delta(i))*Re./Rsc );		%[rad] opposite angle from Re
   Theta(i) = 2*(pi-delta(i)-lambda);			%[rad] orbit portion for G/S contact
end

