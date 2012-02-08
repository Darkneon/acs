%  getD1900 -	Compute number of julian days since Jan 0.5 1900 
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 17 July 03
%
%  function [D1900, C1900] = getD1900(YEAR, DAY)
%
%  See also:

function [D1900, C1900] = getD1900(YEAR, DAY)

global JD1900 JD2000
if isempty(JD2000)
   JD2000 = 2451543.5;		%[Julian Day] Epoch time for Jan 0.0 2000.0
end
if isempty(JD1900)
   JD1900 = 2415020.0;		%[Julian Day] Epoch time for Jan 0.5 1900.0
end

JD = JD2000 + mod(YEAR,2000)*365 + floor(mod(YEAR,2000)/4)+1 + DAY;
D1900 = JD - JD1900;	%[day] from Jan 0.5 1900
C1900 = D1900./36525;	%[Julian century] from Jan 0.5 1900
