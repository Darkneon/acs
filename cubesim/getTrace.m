%  getTrace -	Compute ground trace subpoint vector array
%		(LON, LAT) of an orbiting object
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 29 Dec 02
%
%  LL = getTrace(incl, raan, u, tVE)
%	incl :	Orbit Inclination (scalar)  [rad]
%	raan :	Right Ascention of Ascending Node (array)  [rad]
%	u :	Argument of Perigee + Orbit Motion (array)  [rad]
%	tVE :	Time (row-vector) referenced to Vernal Equinox [s]
%
%  See also: getIntersect

function LL = getTrace2(incl, raan, u, tVE)

global we
if isempty(we)
   we = 0.00417807*pi/180;			%[rad/s] 
end

LAT = asin(sin(incl).*sin(u));			%[rad]
LONu = atan(cos(incl).*tan(u));			%[rad]
for (i=1:1:length(LAT));
   if ((LAT(i)*sign(incl) >= 0) & (LONu(i) < 0))
      LONu(i) = LONu(i)+pi;
   end
   if ((LAT(i)*sign(incl) < 0) & (LONu(i) > 0))
      LONu(i) = LONu(i)-pi;
   end
end   

LON = raan + LONu - we.*tVE;			%[rad]

LL = [LON; LAT];				%[rad]
