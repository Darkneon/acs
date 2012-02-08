%  getPtError -	Compute pointing angle and range error
%		of a S/C with attitude stabilization
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 29 Dec 02
%
%  [AGLtg2pt, Rtg2pt] = getPtError(LLtg, LLpt, Rpt, Rsc)
%	LLtg : 	[LON; LAT] Target location [rad]
%	LLpt :	[LON; LAT] Point location array [rad]
%	Rpt2sc : distance array of Point from S/C [km]
%	Rtg2sc : distance array of S/C from Target [km]
%	AGLtg2pt : S/C Pointing angle error array [rad]
%	Range :	Pointing error distance array [km]
%
%  See also: getMegVec, getIntersect

function [AGLtg2pt, Rtg2pt] = getPtError(LLtg, LLpt, Rpt2sc, Rtg2sc)

global Re
if isempty(Re)
   Re = 6378.13619;		%[km]
end

XYZpt = LL2XYZ(LLpt) .* Re;
XYZtg = LL2XYZ(LLtg) .* Re;

Xtg2pt =  XYZpt(1,:) - XYZtg(1,1);
Ytg2pt =  XYZpt(2,:) - XYZtg(2,1);
Ztg2pt =  XYZpt(3,:) - XYZtg(3,1);
XYZtg2pt_geo = [Xtg2pt; Ytg2pt; Ztg2pt];	%Earth referential

Rtg2pt = LEN(XYZtg2pt_geo);			%[km] spatial distance 

AGLtg2pt = acos ( (Rtg2pt.^2 - Rtg2sc.^2 - Rpt2sc.^2) ./ (-2.*Rtg2sc.*Rpt2sc)); %[rad]
for i = 1:length(Rtg2pt);
   if ( ~isreal(AGLtg2pt(i)) )
      AGLtg2pt(i) = NaN;		
   end
   if ( ~isreal(Rtg2pt(i)) )
      Rtg2pt(i) = NaN;
   end
end

%XYZgsArray = XYZgs*ones(1,length(XYZtg(1,:)));
%theta = acos( dot(XYZgsArray, XYZtg)./Re^2 );	%
%Arc = Re*theta;				%[km] ground distance (arc of earth)
