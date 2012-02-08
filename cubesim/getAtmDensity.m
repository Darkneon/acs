%   getAtmDensity - 	Compute atmospheric density in altitude
%			using exponential intrapolation
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 18 July 03
%
%  [Density, ScaleHt, DensityMax, DensityMin] = getAtmDensity(alt)
%	alt 		: Altitude [km]
%	Density 	: Mean density [kg/m³]
%	ScaleHt 	: Atmospheric Scale Height [km]
%	DensityMax 	: Max density [kg/m³]
%	DensityMin 	: Min Density [kg/m³]

function [Density, ScaleHt, DensityMax, DensityMin] = getAtmDensity(alt)

%atmospheric density (SMAD, WERTZ)
%Altitude   Min         Mean        Max
%(km)       (kg/m³)     (kg/m³)     (kg/m³)
%0          1.225       1.225       1.225
%50         1.77e-2	1.77e-2     1.77e-2
%100        4.61e-7     4.79e-7     5.10e-7
%150        1.65e-9     1.81e-9     2.04e-9
%200        1.78e-10    2.53e-10    3.52e-10
%250        3.35e-11    6.24e-11    1.06e-11
%300        8.19e-12    1.95e-11    3.96e-11
%350        2.34e-12    6.98e-12    1.66e-11
%400        7.32e-13    2.72e-12    7.55e-12
%450        2.47e-13    1.13e-12    3.61e-12
%500        8.98e-14    4.89e-13    1.80e-12
%550        3.63e-14    2.21e-13    9.25e-13
%600        1.68e-14    1.04e-13    4.89e-13
%650        9.14e-15    5.15e-14    2.64e-13
%700        5.74e-15    2.72e-14    1.47e-13
%750        3.99e-15    1.55e-14    8.37e-13
%800        2.96e-15    9.63e-15    4.39e-14
%850        2.28e-15    6.47e-15    3.00e-14
%900        1.80e-15    4.66e-15    1.91e-14
%950        1.44e-15    3.54e-15    1.27e-14
%1000       1.17e-15    2.79e-15    8.84e-15
%1250       4.67e-16    1.11e-15    2.59e-15
%1500       2.30e-16    5.21e-16    1.22e-15

AltDensity = [
1.225       1.225       1.225
1.77e-2     1.77e-2     1.77e-2
4.61e-7     4.79e-7     5.10e-7
1.65e-9     1.81e-9     2.04e-9
1.78e-10    2.53e-10    3.52e-10
3.35e-11    6.24e-11    1.06e-11
8.19e-12    1.95e-11    3.96e-11
2.34e-12    6.98e-12    1.66e-11
7.32e-13    2.72e-12    7.55e-12
2.47e-13    1.13e-12    3.61e-12
8.98e-14    4.89e-13    1.80e-12
3.63e-14    2.21e-13    9.25e-13
1.68e-14    1.04e-13    4.89e-13
9.14e-15    5.15e-14    2.64e-13
5.74e-15    2.72e-14    1.47e-13
3.99e-15    1.55e-14    8.37e-13
2.96e-15    9.63e-15    4.39e-14
2.28e-15    6.47e-15    3.00e-14
1.80e-15    4.66e-15    1.91e-14
1.44e-15    3.54e-15    1.27e-14
1.17e-15    2.79e-15    8.84e-15
4.67e-16    1.11e-15    2.59e-15
2.30e-16    5.21e-16    1.22e-15]; %[kg/m^3]

%Atmospheric Scale Height by 50km steps 
%from ground level to 1000km (SMAD)
AtmScale =[
   8.44
   7.95
   5.9
   25.5
   37.5
   44.8
   50.3
   54.8
   58.2
   61.3
   64.5
   68.7
   74.8
   84.4
   99.3
   121
   151
   188
   226
   263
   296
   408
   516];
   
%h=linspace(0,1500,200);		%for graph plot below

% Linear interpolation
%if alt>1000
%   Density = 0;
%   ScaleHt = 1000;
%else
%   DensityHi = AltDensity(ceil(alt/50)+1);
%   DensityLo = AltDensity(floor(alt/50)+1);
%   Density = (DensityHi - DensityLo)* (mod(alt,50))/50 + DensityLo; %linear interpolation
%   
%   ScaleHtHi = AtmScale(ceil(alt/50)+1);
%   ScaleHtLo = AtmScale(floor(alt/50)+1);
%   ScaleHt = (ScaleHtHi - ScaleHtLo)* (mod(alt,50))/50 + ScaleHtLo; %linear interpolation
%end

% Exponential interpolation/extrapolation
if alt>1000 & alt<1500
   DensityHi = AltDensity(ceil((alt-1000)/250+0.0001)+21, 2);
   DensityLo = AltDensity(floor((alt-1000)/250)+21, 2);
   k = 1/250*log(DensityLo/DensityHi);
   A = DensityLo/exp(-k*250*floor(alt/250));
   Density = A*exp(-k*alt); 		%exponential interpolation
   
   DensityMaxHi = AltDensity(ceil((alt-1000)/250+0.0001)+21, 3);
   DensityMaxLo = AltDensity(floor((alt-1000)/250)+21, 3);
   k = 1/250*log(DensityMaxLo/DensityMaxHi);
   A = DensityMaxLo/exp(-k*250*floor(alt/250));
   DensityMax = A*exp(-k*alt); 		%exponential interpolation

   DensityMinHi = AltDensity(ceil((alt-1000)/250+0.0001)+21, 1);
   DensityMinLo = AltDensity(floor((alt-1000)/250)+21, 1);
   k = 1/250*log(DensityMinLo/DensityMinHi);
   A = DensityMinLo/exp(-k*250*floor(alt/250));
   DensityMin = A*exp(-k*alt); 		%exponential interpolation
   
   ScaleHtHi = AtmScale(ceil((alt-1000)/250+0.0001)+21);
   ScaleHtLo = AtmScale(floor((alt-1000)/250)+21);
   k = 1/250*log(ScaleHtLo/ScaleHtHi);
   A = ScaleHtLo/exp(-k*250*floor(alt/250));
   ScaleHt = A*exp(-k*alt); 		%exponential interpolation
elseif alt >=1500
   DensityHi = AltDensity(length(AltDensity), 2);
   DensityLo = AltDensity(length(AltDensity)-1, 2);
   k = 1/250*log(DensityLo/DensityHi);
   A = DensityHi/exp(-k*1500);
   Density = A*exp(-k*alt); 		%exponential extrapolation
   
   DensityMaxHi = AltDensity(length(AltDensity), 3);
   DensityMaxLo = AltDensity(length(AltDensity)-1, 3);
   k = 1/250*log(DensityMaxLo/DensityMaxHi);
   A = DensityMaxHi/exp(-k*1500);
   DensityMax = A*exp(-k*alt); 		%exponential extrapolation
   
   DensityMinHi = AltDensity(length(AltDensity(:,1)), 1);
   DensityMinLo = AltDensity(length(AltDensity(:,1))-1, 1);
   k = 1/250*log(DensityMinLo/DensityMinHi);
   A = DensityMinHi/exp(-k*1500);
   DensityMin = A*exp(-k*alt); 		%exponential extrapolation
   
   ScaleHtHi = AtmScale(length(AltDensity));
   ScaleHtLo = AtmScale(length(AltDensity)-1);
   k = 1/250*log(ScaleHtLo/ScaleHtHi);
   A = ScaleHtHi/exp(-k*1500);
   ScaleHt = A*exp(-k*alt); 		%exponential extrapolation
else
   DensityHi = AltDensity(ceil(alt/50+0.0001)+1, 2);
   DensityLo = AltDensity(floor(alt/50)+1, 2);
   k = 1/50*log(DensityLo/DensityHi);
   A = DensityLo/exp(-k*50*floor(alt/50));
   Density = A*exp(-k*alt); 		%exponential interpolation
   
   DensityMaxHi = AltDensity(ceil(alt/50+0.0001)+1, 3);
   DensityMaxLo = AltDensity(floor(alt/50)+1, 3);
   k = 1/50*log(DensityMaxLo/DensityMaxHi);
   A = DensityMaxLo/exp(-k*50*floor(alt/50));
   DensityMax = A*exp(-k*alt); 		%exponential interpolation
   
   DensityMinHi = AltDensity(ceil(alt/50+0.0001)+1, 1);
   DensityMinLo = AltDensity(floor(alt/50)+1, 1);
   k = 1/50*log(DensityMinLo/DensityMinHi);
   A = DensityMinLo/exp(-k*50*floor(alt/50));
   DensityMin = A*exp(-k*alt); 		%exponential interpolation
   
   ScaleHtHi = AtmScale(ceil(alt/50+0.0001)+1);
   ScaleHtLo = AtmScale(floor(alt/50)+1);
   k = 1/50*log(ScaleHtLo/ScaleHtHi);
   A = ScaleHtLo/exp(-k*50*floor(alt/50));
   ScaleHt = A*exp(-k*alt); 		%exponential interpolation
end



%F1=figure;
%set(F1, 'Position', [0 225 1300 700])
%title('Air density in altitude')
%xlabel('Density [kg/m³]')
%ylabel('Altitude [km]')
%semilogy(alt,AltDensity)
%grid;
