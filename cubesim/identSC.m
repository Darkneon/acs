%  identSC	Identify location of orbiting object
%		from ground station tracking inputs
%
%	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 29 Dec 02
%
%  [LLsc, Rsc] = identSC(LLgs, AziEl, Range)
%	LLgs : 	[LON; LAT] G/S location [rad]
%	AziEl :	[Azimut; Elevation] column-vector array
%		of object (S/C) to identify [rad]
%	Range :	distance array of object S/C from G/S [km]
%	LLsc :	[LON; LAT] S/C location array [rad]
%	Rsc :	distance array of object S/C from earth center [km]
%	XYZsc :	[X; Y; Z] column-vector array
%		of object (S/C) to identify [km]
%
%  See also: trackSC, getRange

function [LLsc, Rsc, XYZsc_earth] = identSC(LLgs, AziEl, Range)

global Re
if isempty(Re)
   Re = 6378.13619;		%[km]
end

XYZgs_earth = LL2XYZ(LLgs) * Re;

XYZgs2sc_gs = AziEl2XYZ(AziEl, Range);
ROTgs_earth = Rx(pi)*Ry(-LLgs(2,1)-pi/2)*Rz(LLgs(1,1));

disp('ROTgs_earth^(-1):'); %TEST LINE
disp(ROTgs_earth^(-1));

%TEST START

disp('XYZgs2sc_gs :');
disp(XYZgs2sc_gs);
disp('END TEST');
%TEST END

XYZgs2sc_earth = (ROTgs_earth^(-1)) * XYZgs2sc_gs'; %ERROR ON THIS LINE

Xsc_earth =  XYZgs2sc_earth(1,:) + XYZgs_earth(1,1);
Ysc_earth =  XYZgs2sc_earth(2,:) + XYZgs_earth(2,1);
Zsc_earth =  XYZgs2sc_earth(3,:) + XYZgs_earth(3,1);
XYZsc = [Xsc_earth; Ysc_earth; Zsc_earth];

[LLsc, Rsc] = XYZ2LL(XYZsc);
