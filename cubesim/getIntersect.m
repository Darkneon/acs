%  getIntersect	Compute intersect subpoint vector array
%		(LON, LAT) of an orbiting object pointing
%		to the spherical Earth
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 16 July 03
%
%  [LLtg, Rtg] = getTarget(LLsc, Rsc, Vtg)
%	LLsc  :	[LON; LAT] S/C location array [rad]
%	Rsc   :	distance array of object S/C from earth center [km]
%	Vtg   :	[X; Y; Z] column-vector array
%		of attitude pointing resolved in geo coordinates
%	LLtg  :	[LON; LAT] S/C Location [rad]
%	Rtg   :	distance array of ground Target from S/C [km]
%
%  See also: getTrace, getMagField, getPtError

function [LLtg, Rtg] = getIntersect(LLsc, Rsc, Vtg)

global Re
if isempty(Re)
   Re = 6378.13619;		%[km]
end

Nvtg = LEN(Vtg);
Vtg = Vtg./ ([1;1;1]*Nvtg);	%normalize pointing vector

IJKsc = LL2XYZ(LLsc);		%

gamma = acos(dot(-IJKsc, Vtg));		%[rad] angle between position vector and pointing vector

% by cos law:  Re^2 = Rsc^2 + Rtg^2 - 2*Rsc*Rtg*cos(gamma)
% the solved as = (-b -sqrt(b^2-4ac)) /2a
Rtg = Rsc.*cos(gamma) - sqrt( (Rsc.*cos(gamma)).^2 - (Rsc.^2-Re^2) ); %[km] distance of ground Target from S/C

% trim for out of earth surface intersect
for (i=1:1:length(Rtg))
   if (~isreal(Rtg(i)) | Rtg(i) < 0)	%Target intersects out of Earth surface
      Rtg(i) = NaN;			%=> exclude value
   end
end   

XYZtg = [Rsc;Rsc;Rsc].*IJKsc + [Rtg;Rtg;Rtg].*Vtg;
LLtg = XYZ2LL(XYZtg);
