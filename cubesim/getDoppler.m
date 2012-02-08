%  getDoppler	Determine Doppler shift and Rate
%
%	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 29 Dec 02
%
%  [Doppler, DRate] = getDoppler(XYZgs2sc, fs, f0, filter)
%	XYZgs2sc :	[X; Y; Z] collumn-vector array
%			of object S/C from G/S [km]
%	fs :	sampling frequency (scalar)
%	f0 : 	communication link frequency (scalar)
%	Doppler:	frequency shift array [Hz]
%	DRate :	frequency shift rate array [Hz/s]
%
%  See also: getRange

function [Doppler, DRate] = getDoppler(XYZgs2sc, fs, f0)

global Re
if isempty(Re)
   Re = 6378.13619;	%[km]
end

global c0
if isempty(c0)
   c0 = 3e+8;		%[m/s]
end

Len = length(XYZgs2sc(1,:));
Vsc = zeros(3,Len);
Vd = zeros(1,Len);
DRate = zeros(1,Len);

for i = 1:Len-1
    Vsc(:,i) = (XYZgs2sc(:,i+1) - XYZgs2sc(:,i)) .*fs.*1000;		%[m/s] absolute velocity vector
    Vsc2gs(i) = dot(Vsc(:,i), XYZgs2sc(:,i)) /LEN(XYZgs2sc(:,i));	%[m/s] projected velocity vector
end
Doppler = [-Vsc2gs * (f0/c0), 0];			%[Hz] Doppler shift
	
for i = 1:Len-2
    DRate(i) = (Doppler(i+1) - Doppler(i)) .*fs;	%[Hz/s] Doppler shift rate
end

