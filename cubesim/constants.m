%  constants -	CubeSat Orbital Analysis Constants Definition
%
% 	Jean-Francois Levesque, MS
%  	jflev@yahoo.ca
%  	Last Update: 21 July 03

%============================================================
% Constants at Epoch 2000
%============================================================
global JD2000 JD1900;
JD2000 = 2451543.5;		%[Julian Day] Epoch time for Jan 0.0 2000.0
JD1900 = 2415020.0;		%[Julian Day] Epoch time for Jan 0.5 1900.0
GST1900 = JD2000;		%[Julian Day] Current Epoch time for Jan 2000.0
T = (JD2000-2451545.0)/36525;	%[Julian Century] Curent Epoch time (SMAD p.898)

%---------------------------------
% Astronomical Constants
%---------------------------------
global c0 AU Te
c0 = 2.99792458e+8;		%[m/s] light speed in vacuum
AU = 1.49597870691e+8;		%[km] astronomical unit
Te = 86400;			%[s/day] one day

%---------------------------------
% Sun Constants
%---------------------------------
GMsun = 1.327124e+11;		%[km3/s2] Gravity parameter (mu)
msun = 1.9891e+30;		%[kg] Mass
Rsun = 6.95508e+5;		%[km] Mean Radius of Photosphere (+-0.00026e+5)

%---------------------------------
% Earth Constants
%---------------------------------
global GMe me incle
GMe = 398600.4418;		%[km3/s2] Gravity parameter (mu)
me = 5.9737e+24; 		%[kg] Mass
incle = (23.439281083-46.815/3600*T)*pi/180;	%[rad] Obliquity of Ecliptic

global Re we Pe
Re = 6378.13649;		%[km] Radius at equator
Re_polar = 6356.7517;		%[km]
Re_mean = 6371.0003;		%[km]
we = 7.2921158553e-5;		%[rad/s] Rotation velocity
Pe = 2*pi/we;			%[s/rev] Rotation period

global B_EARTH uo SGalnico SGhymu J2
B_EARTH = 3e-5;		%[Tesla] Earth's magnetic induction (field strength)
uo = 4*pi*1e-7;	 	%[1/A²] permittivity of vacuum
SGalnico = 7.010;
SGhymu = 8.740;

ECCsun = 0.016751;		%[] Earth orbital eccentricity
S = 1.000001057*AU;		%[km] Mean distance to sun
Rpsun = S * (1-ECCsun);		%[km] Perihelion
Rasun = S * (1+ECCsun);		%[km] Aphelion
Psun = 365.256*Te;		%[s/rev] Revolution period
wsun = 2*pi/Psun;		%[rad/s] Revolution velocity

J2 = 0.00108263;		%[] Second Geopotential Parameter


%---------------------------------
% Time Constants (Equinoxes, Solstices, Peri/Aphelions) 
%---------------------------------
%2000  3/20  7:30  6/21  1:42  9/22 17:16  12/21 13:28   1/04  0:02   7/04 14:57
%2001  3/20 13:19  6/21  7:30  9/22 23:04  12/21 19:17   1/03  6:16   7/04 21:10
%2002  3/20 19:08  6/21 13:18  9/23  4:53  12/22  1:07   1/03 12:30   7/05  3:24
%2003  3/21  0:58  6/21 19:06  9/23 10:42  12/22  6:57   1/03 18:44   7/05  9:38
%2004  3/20  6:47  6/21  0:55  9/22 16:30  12/21 12:47   1/04  0:58   7/04 15:52
%2005  3/20 12:36  6/21  6:43  9/22 22:19  12/21 18:36   1/03  7:12   7/04 22:06
%2006  3/20 18:25  6/21 12:31  9/23  4:08  12/22  0:26   1/03 13:26   7/05  4:20
%2007  3/21  0:14  6/21 18:19  9/23  9:56  12/22  6:16   1/03 19:40   7/05 10:34
%2008  3/20  6:04  6/21  0:07  9/22 15:45  12/21 12:05   1/04  1:54   7/04 16:48
%2009  3/20 11:53  6/21  5:55  9/22 21:34  12/21 17:55   1/03  8:08   7/04 23:02
%2010  3/20 17:42  6/21 11:43  9/23  3:23  12/21 23:45   1/03 14:22   7/05  5:16
global VERNAL_EQX APHELION
VERNAL_EQX = [(20+28+31)+(13+19/60)/24,...	%2001 (indice = 1 = year-2000)
      (20+28+31)+(19+08/60)/24,...		%2002 (indice = 2)
      (21+28+31)+(0+58/60)/24,...
      (20+29+31)+(6+47/60)/24,...
      (20+28+31)+(12+36/60)/24,...
      (20+28+31)+(18+25/60)/24,...
      (21+28+31)+(0+14/60)/24,...
      (20+29+31)+(6+04/60)/24,...
      (20+28+31)+(11+53/60)/24,...
      (20+28+31)+(17+42/60)/24,...
         ];	                        %[day] approx from January 1, 0:00

APHELION = [(4+30+31+30+31+28+31)+(21+10/60)/24,...	%2001 (indice = 1 = year-2000)
      (5+30+31+30+31+28+31)+(3+24/60)/24,...		%2002 (indice = 2)
      (5+30+31+30+31+28+31)+(9+38/60)/24,...
      (4+30+31+30+31+29+31)+(15+52/60)/24,...
      (4+30+31+30+31+28+31)+(22+06/60)/24,...
      (5+30+31+30+31+28+31)+(4+20/60)/24,...
      (5+30+31+30+31+28+31)+(10+34/60)/24,...
      (4+30+31+30+31+29+31)+(16+48/60)/24,...
      (4+30+31+30+31+28+31)+(23+02/60)/24,...
      (5+30+31+30+31+28+31)+(5+16/60)/24,...
         ];	                        %[day] approx from  January 1, 0:00

SummerSol = (21+31+30+31+28+31)+(19+6/60)/24;	                %[day] approx
AutomnalEq = (23+31+31+30+31+30+31+28+31)+(10+42/60)/24;	%[day] approx
WinterSol = (22+30+31+30+31+31+30+31+30+31+28+31)+(6+57/60)/24; %[day] approx
Aphelion = (5+30+31+30+31+28+31)+(9+38/60)/24;	                %[day] approx
Perihelion = (3)+(18+44/60)/24;	                                %[day] approx

%---------------------------------
% Atmospheric constants
%---------------------------------
Po = 101325;     %[Pa] standard pressure at sea level
To = 288.15;     %[K] standard temperature at sea level
g = 9.80665;     %[m/s²] gravitational constant
L = 6.5;         %[K/km] temperature lapse rate
R = 8.31432;     %[J/mol.K] gas constant
M = 28.9644;     %[g/mol] gas molecular weight

%---------------------------------
% Radiation Constants
%---------------------------------
global Gsolar Albedo GIR Boltzmann
Gsolar = 1367;		%[W/m2] Sun 
Gsolar_min = 1326;	%[W/m2] Sun 
Gsolar_max = 1418;	%[W/m2] Sun
SOLAR_YEAR = 2000.7;
global Albedo;
Albedo = 0.3;		%[] Sun Albedo Ratio
global GIR;
GIR = 237;		%[W/m2] Earth IR
global Boltzmann;
Boltzmann = 5.670e-8;   	%[W/m^2-K^4] Boltzmann constant

% --- Radiation Coefficients ---
emi_alu = 0.0346;
abs_alu = 0.379;
emi_bk = 0.874;
abs_bk = 0.975;
emi_cell = 0.825;
abs_cell = 0.805;
emi_gold = 0.023;
abs_gold = 0.299;

%---------------------------------
% Magnetic Poles
%---------------------------------
global LON_MAG_NORTH LAT_MAG_NORTH LON_MAG_SOUTH LAT_MAG_SOUTH
LON_MAG_SOUTH = -104.0*pi/180;		%[rad]
LAT_MAG_SOUTH = 78.3*pi/180;		%[rad]
LON_MAG_NORTH = LON_MAG_SOUTH +pi;	%[rad]
LAT_MAG_NORTH = -LAT_MAG_SOUTH;		%[rad]

%http://aom.giss.nasa.gov/srvernal.html
%
%Orbital Events    Tropical Year = 365.2425 (days)    Greenwich Mean Time is used
%      Vernal      Summer      Autumnal     Winter    
%Year     Equinox    Solstace     Equinox     Solstace   Perihelion    Aphelion 
%----  ----------  ----------  ----------   ----------   ----------   ----------
%2000  3/20  7:30  6/21  1:42  9/22 17:16  12/21 13:28   1/04  0:02   7/04 14:57
%2001  3/20 13:19  6/21  7:30  9/22 23:04  12/21 19:17   1/03  6:16   7/04 21:10
%2002  3/20 19:08  6/21 13:18  9/23  4:53  12/22  1:07   1/03 12:30   7/05  3:24
%2003  3/21  0:58  6/21 19:06  9/23 10:42  12/22  6:57   1/03 18:44   7/05  9:38
%2004  3/20  6:47  6/21  0:55  9/22 16:30  12/21 12:47   1/04  0:58   7/04 15:52
%2005  3/20 12:36  6/21  6:43  9/22 22:19  12/21 18:36   1/03  7:12   7/04 22:06
%2006  3/20 18:25  6/21 12:31  9/23  4:08  12/22  0:26   1/03 13:26   7/05  4:20
%2007  3/21  0:14  6/21 18:19  9/23  9:56  12/22  6:16   1/03 19:40   7/05 10:34
%2008  3/20  6:04  6/21  0:07  9/22 15:45  12/21 12:05   1/04  1:54   7/04 16:48
%2009  3/20 11:53  6/21  5:55  9/22 21:34  12/21 17:55   1/03  8:08   7/04 23:02
%2010  3/20 17:42  6/21 11:43  9/23  3:23  12/21 23:45   1/03 14:22   7/05  5:16


%atmospheric density (SMAD)
%Altitude   Min         Mean        Max
%(km)       (kg/m³)     (kg/m³)     (kg/m³)
%0          1.2         1.2         1.2
%50         
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

%---------------------------------------
% Ground Station Parameters
%---------------------------------------
% Weingarten DE
%LON_WG = 9+37/60;			%[deg]
%LAT_WG = 48+46/60;			%[deg]
%MEL_WG = 10;				%[deg] minimum horizon inclinaison
%LL_WG = [LON_WG; LAT_WG]*pi/180;	%[rad]
%S_WG = SCrange(MEL_WG, rp);   %[km]

% Wurzburg DE
%LON_WU = 9+55/60;			%[deg]
%LAT_WU = 49+47/60;			%[deg]
%MEL_WU = 10;				%[deg] minimum horizon inclinaison
%LL_WU = [LON_WU; LAT_WU]*pi/180;	%[rad]
%S_WU = SCrange(MEL_WU, rp);   %[km]

% New York, USA
%LON_NY = -74;				%[deg]
%LAT_NY = 40.7;				%[deg]
%MEL_NY = 10;				%[deg] minimum horizon inclinaison
%LL_NY = [LON_NY; LAT_NY]*pi/180;	%[rad]
%S_NY = SCrange(MEL_NY, rp);   %[km]

% Palo Alto, USA
%LON_PA = -122.12;			%[deg]
%LAT_PA = 37.41;				%[deg]
%MEL_PA = 10;				%[deg] minimum horizon inclinaison
%LL_PA = [LON_PA; LAT_PA]*pi/180;	%[rad]
%S_PA = SCrange(MEL_PA, rp);   %[km]



%============================================================
