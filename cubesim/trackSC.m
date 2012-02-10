%  trackSC	Identify ground station tracking parameters
%		from orbiting object location inputs
%
%	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 29 Dec 02
%
%  [AziEl, Range, XYZgs2sc_earth] = trackSC(LLgs, LLsc, Rsc, filter)
%	LLgs : 	[LON; LAT] G/S location [rad]
%	LLsc :	[LON; LAT] S/C location array [rad]
%	Rsc :	distance array of object S/C from earth center [km]
%	filter : scalar parameter to trim data when out of LOS
%	AziEl :	[Azimut; Elevation] column-vector array
%		of object (S/C) to identify [rad]
%	Range :	distance array of object S/C from G/S [km]
%	XYZgs2sc_earth : [X; Y; Z] column-vector array
%			of object S/C to track
%			in G/S referential [km]
%
%  See also: identSC, getRange

function [AziEl, Range, XYZgs2sc_earth] = trackSC(LLgs, LLsc, Rsc, filter)

global Re
if isempty(Re)
   Re = 6378.13619;		%[km]
end

XYZgs_earth = LL2XYZ(LLgs) * Re;
XYZsc_earth = LL2XYZ(LLsc, Rsc);

Xgs2sc_earth =  XYZsc_earth(1,:) - XYZgs_earth(1,1);
Ygs2sc_earth =  XYZsc_earth(2,:) - XYZgs_earth(2,1);
Zgs2sc_earth =  XYZsc_earth(3,:) - XYZgs_earth(3,1);
XYZgs2sc_earth = [Xgs2sc_earth; Ygs2sc_earth; Zgs2sc_earth];	%Earth referential

Range = LEN(XYZgs2sc_earth);
disp('size of the range');
disp(size(Range));

ROTgs_earth = Rx(pi)*Ry(-LLgs(2,1)-pi/2)*Rz(LLgs(1,1));
XYZgs2sc_gs = ROTgs_earth * XYZgs2sc_earth;        		%G/S referential

AziEl = XYZ2AziEl(XYZgs2sc_gs);        				%G/S referential

% exclude out-of-view data if 'filter' is requested
if (filter)				
   for (i = 1:length(Range) )
      if ( (AziEl(2,i) < 0) | ~isreal(AziEl(2,i)))
         AziEl(:,i) = NaN;		
         Range(i) = NaN;
      end
   end
end
