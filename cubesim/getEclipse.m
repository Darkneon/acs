%  getEclipse -	Compute S/C eclipses occurence and durations 
%
%  Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 17 July 03
%
%  [State, VF, Eclipse_time, Eclipse_dur] = getEclipse(XYZsun, XYZsc, Ts)
%	XYZsun 		: [X; Y; Z] column-vector array [km]
%	XYZsc 		: [X; Y; Z] column-vector array [km]
%	Ts 		: sampling period [s] (scalar)
%	State 		: Eclipse state array 
%			  0 = No eclipse
%			  1 = Partial eclipse
%			  2 = Total eclipse
%			  3 = Anular eclipse
%	VF 		: View factor of Sun for intersity calculation
%			  between 0 and 1 (0 = total eclipse, 1 = complete Sun)
%	Eclipse_time 	: Eclipse occurence times [s] (array)
%	Eclipse_dur 	: Eclipse occurence duration [s] (array)
%
%  See also: getShadow

function [State, VF, Eclipse_time, Eclipse_dur] = getEclipse(XYZsun, XYZsc, Ts)

global Re Rsun
if isempty(Re)
   Re = 6378.13619;			%[km]
end
if isempty(Rsun)
   Rsun = 6.95508e+5;			%[km]
end

% Astral bodies relative to S/C (in earth coordinates though)
XYZearth_sc = -XYZsc;
XYZsun_sc = XYZearth_sc + XYZsun;
Dearth_sc = LEN(XYZearth_sc);
Dsun_sc = LEN(XYZsun_sc);
Dsun_earth = LEN(XYZsun);
C = Re.*Dsun_earth ./(Rsun - Re);

N=length(Dsun_sc);

% Angular radius of bodies (WERTZ p.76)
rhos = asin(Rsun./Dsun_sc);
rhoe = asin(Re./Dearth_sc);
theta = acos(dot(XYZearth_sc./([1;1;1]*Dearth_sc), XYZsun_sc./([1;1;1]*Dsun_sc) ));

% eclipse computation, 4 states: (0)No, (1)Partial, (2)Total, (3)Annular (WERTZ p.76)
State = zeros(1,N);  % set initialy  all value as no eclipse
VF = ones(1,N);	% Sun Area View Factor (also intersity factor) initially 1 (no eclipse)
for i=1:N
   if (Dsun_sc(i) > Dsun_earth(i)) & (theta < rhos(i) + rhoe(i)) & (theta(i) > abs(rhoe(i) - rhos(i)))
      State(i) = 1;
      a = (cos(rhoe(i)) - cos(rhos(i))*cos(theta(i))) /sin(rhos(i))/sin(theta(i));
      b = (cos(rhos(i)) - cos(rhoe(i))*cos(theta(i))) /sin(rhoe(i))/sin(theta(i));
      c = (cos(theta(i)) - cos(rhos(i))*cos(rhoe(i))) /sin(rhos(i))/sin(rhoe(i));
      VF(i) = 1 - ((pi - cos(rhos(i))*acos(a) - cos(rhoe(i))*acos(b) - acos(c)) /pi/(1-cos(rhos(i))));
   elseif (Dsun_sc(i) > Dsun_earth(i)) & (Dsun_sc(i) < Dsun_earth(i)+ C(i)) & (theta(i) < rhoe(i) - rhos(i))
      State(i) = 2;
      VF(i) = 0;
   elseif (Dsun_sc(i) > Dsun_earth(i) + C(i)) & (theta(i) < rhoe(i) - rhos(i))
      State(i) = 3; 
      VF(i) = 1 - ((1-cos(rhoe(i))) / (1-cos(rhos(i))) );
   end
end

%Eclipse times cumulation, assume any form of eclipse as full eclipse
Eclipse_dur = zeros(10,1);			%initial matrices
Eclipse_time = zeros(10,1);
j=0;
previous = 'sun';
for i=1:N
   if State(i)					%1,2,3 = eclipse, 0 = no eclipse
      if (previous == 'sun')
         j = j+1;
         if j>10				%enlarge initial matrices if more than 10 occurences
            Eclipse_dur = [Eclipse_dur; 0];
            Eclipse_time = [Eclipse_time; 0];            
         end
         Eclipse_time(j) = i*Ts - Ts/2;			%time when eclipse begins (half sample time before)
	 Eclipse_dur(j) = Eclipse_dur(j) + Ts;	%2xhalf-sample durations added for begining and end compensation
         previous = 'ecl';
      end
      Eclipse_dur(j) = Eclipse_dur(j) + Ts;	%sample time is cumulated until exit of eclipse
   else
      previous = 'sun';
   end
end

for i=1:length(Eclipse_dur)
   if (Eclipse_dur(i) == 0)
%      Eclipse_dur(i) = NaN;
      Eclipse_time(i) = NaN;		%remove unused matrices cells for graph plot purpose
   end
end
