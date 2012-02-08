%  getShadow -	Compute shadow area and solar ray incidence
%		on orbital plane
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 29 Dec 02
%
%  [beta, theta] = getShadow(LLsun, LLorbit, Rsc)
%	LLsun :	[Longitude; Latitude] column-vector array [rad]
%	LLorbit : [Longitude; Latitude] column-vector array [rad]
%	Rsc :	distance array of S/C from earth center
%	beta :	Solar ray incidence angle array on orbit plane
%	theta :	Shadow portion of orbit (array)
%
%  See also: getEclipse

function [beta, theta] = getShadow(LLsun, LLorbit, Rsc)

global Re
if isempty(Re)
   Re = 6378.13619;		%[km]
end

global EclipseCone
if isempty(EclipseCone)
   EclipseCone = 5*pi/180;		%[rad]
end

IJKsun = LL2XYZ(LLsun);
IJKorbit = LL2XYZ(LLorbit);

% Incidence of Solar Radiation on orbital plane
beta = pi/2 - acos(dot(IJKsun, IJKorbit)) - EclipseCone;	%[rad]

% Orbit portion in Shadow
theta = 2*asin(sqrt( (Re./Rsc)^2 - (sin(beta)).^2 ));	%[rad]
for i=1:1:length(theta)
   if ~isreal(theta(i))					% filter for complex elements
      theta(i) = 0;					% (out of shadow area)
   end
end

