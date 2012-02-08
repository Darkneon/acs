%  plotTrace	Plot ground trace from a data set 
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 03 Jan 03
%
%  H = plotTrace(LON, DATA, type)
%	LON : Longitude array [deg]
%	DATA : Data array
%	type :	graphic plot type
%	H :	graphic handle
%  
%  See also: getShadow, trimTrace

function H = plotTrace(LON, DATA, type)

N = length(LON);

if LON(N-1) < LON(N)
   Dir = 1;
else
   Dir = -1;
end

%trim data set for ±180°  LON
for ( j=1:1:(N-1) )
   if isnan(LON(j))
   else
      if  ( ( (LON(j+1)<LON(j)-90) & (Dir>0) ) | ( (LON(j+1)>LON(j)+90) & (Dir<0) ) )
         LON = [LON(1:j), NaN, LON((j+1):N)];
         DATA = [DATA(1:j), NaN, DATA((j+1):N)];
         N = N+1;
      end
      
      if LON(j+1) > LON(j)
         Dir = 1;
      else
         Dir = -1;
      end
   end
end
H = plot(LON, DATA, type);


