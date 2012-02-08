%  getMagField - Compute Earth's magnetic field line vector at location
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 15 July 03
%
%  XYZfield_geo = getMagField(LLsc, Rsc)
%  XYZfield_geo = getMagField(XYZsc_geo)
%	LLsc :	[Longitude, Latitude] column-vector array
%		of spacecraft location (Geo referential)
%       Rsc :   Radial distance from Earth's center
%	XYZsc_geo : column-vector array
%	 	    of spacecraft location (Geo referential)
%	XYZfield_geo : [X; Y; Z] Magnetic Field vector array
%		        at S/C location (Geo referential)
%
%	The calculations are based on the SMAD L-Shell model.
%	The magnetic line angle is independant from the distance
%	of the spacecraft for a given magnetic latitude 
%	--------------------
%	 lambda = LLsc_mag(2,:) : magnetic latitude
%	 R = L.cos²(lambda)	: from SMAD
%	 u = R.cos(lambda)
%	 v = R.sin(lambda)
%	 slope = Dz/Du
%	 phi = atan(Dv/Du)
%	 u = x^2 + y^2
%        x = cos(phi) .cos(LON)
%        y = cos(phi) .sin(LON)
%        z = sin(phi)
%        B(R,lambda) = Bo*Re^3/R^3* (1+sin(lambda)^2)^.5
%	--------------------
%
%  See also: getIntersect

function XYZfield_geo = getMagField(arg1, arg2)

global LON_MAG_SOUTH LAT_MAG_SOUTH B_EARTH Re
if isempty(LON_MAG_SOUTH)
   LON_MAG_SOUTH = -104.0*pi/180;	%[rad]
end
if isempty(LAT_MAG_SOUTH)
   LAT_MAG_SOUTH = 78.3*pi/180;		%[rad]
end
if isempty(B_EARTH)
   B_EARTH = 3e-5;				%[Tesla] Earth's magnetic induction (field strength)
end
if isempty(Re)	
   Re = 6378.13649;			%[km] Radius at equator
end

ROTgeo_mag = Ry(pi/2-LAT_MAG_SOUTH)*Rz(LON_MAG_SOUTH);

if nargin == 2
   LLsc = arg1;
   Rsc = arg2;
   IJKsc_geo = LL2XYZ(LLsc);
elseif nargin == 1
   Rsc = LEN(arg1);
   IJKsc_geo = [arg1(1,:)./Rsc ; arg1(2,:)./Rsc ; arg1(3,:)./Rsc];
end

%Magnetic latitude calculation
IJKsc_mag = ROTgeo_mag * IJKsc_geo;
LLsc_mag = XYZ2LL(IJKsc_mag);
lambda = LLsc_mag(2,:);					%[rad] Magnetic Latitude

% Field direction calculation
du_dlambda = -3*(cos(lambda)).^2.*sin(lambda);
dv_dlambda = (cos(lambda)).^3 - 2*(sin(lambda)).^2.*cos(lambda);

phi = atan2(dv_dlambda, du_dlambda);

x = cos(phi) .*cos(LLsc_mag(1,:));
y = cos(phi) .*sin(LLsc_mag(1,:));
z = sin(phi);

IJKfield_mag = [x ; y ; z];

IJKfield_geo = ROTgeo_mag^(-1) * (IJKfield_mag);	%Magnetic field vector in geo Referential

% Field Strength calculation
B = B_EARTH*Re^3./Rsc.^3.* (1 + sin(lambda).^2).^(1/2);

XYZfield_geo = [IJKfield_geo(1,:).*B ; IJKfield_geo(2,:).*B ; IJKfield_geo(3,:).*B];

