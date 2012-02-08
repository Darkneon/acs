%  getLink	Compute G/S contacts occurence and durations 
%
%	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 12 August 03
%
%  [tlink, Pass_dur, Link_state] = SClink(Elevation, Min_El, Ts)
%  [tlink, Pass_dur, Link_state, Link_dur] = SClink(Elevation, Min_El, Ts, AntPt, Beam)
%	IJKray : [I; J; K] unit column-vector array [rad]
%	IJKsc : [I; J; K] unit column-vector array [rad]
%	Rsc :	distance array of S/C from earth center [km]
%	Ts :	sampling period [s] (scalar)
%	AntPt :	angle array of S/C antenna pointing to G/S loc [rad/s]
%	Beam :	Beam width [rad/s]
%	tlink :	G/S contact occurences times [s] (array)
%	Pass_dur : G/S pass occurences duration [s] (array)
%	Link_state : G/S contact state array (1 = contact, 0 = no-contact)
%	Link_dur : G/S link occurences duration [s] (array)
%
%  See also: trackSC

function [tlink, Pass_dur, Link_state, Link_dur] = SClink(Elevation, Min_El, Ts, AntPt, Beam)

Pass_dur = zeros(1,10);
tlink = zeros(1,10);
Link_state = zeros(1,length(Elevation));
Link_dur = zeros(1,10);

if nargin > 3
   j=0;
   loop = 0;
   L=length(Elevation);
   for i=1:L
      if (Elevation(i) >= Min_El)
         Link_state(i) = 1;		% 0 = no link, 1 = link
         if (loop == 0)
            j = j+1;
            if j>10
               Pass_dur = [Pass_dur, 0];
               Link_dur = [Link_dur, 0];
               tlink = [tlink, 0];            
            end
            tlink(j) = i*Ts - Ts/2;		%[s]
            Pass_dur(j) = Pass_dur(j) + Ts;	%2 half-samples are added
            Link_dur(j) = Link_dur(j) + Ts;	%2 half-samples are added
            loop = 1;
         end
         Pass_dur(j) = Pass_dur(j) + Ts;
         if AntPt(i) <= Beam/2
            Link_dur(j) = Link_dur(j) + Ts;
         else
            Link_state(i) = 0;		% override
         end
      else
         loop = 0;
      end
   end
   
else
   j=0;
   loop = 0;
   L=length(Elevation);
   for i=1:L
      if (Elevation(i) >= Min_El)
         Link_state(i) = 1;		% 0 = no link, 1 = link
         if (loop == 0)
            j = j+1;
            if j>10
               Pass_dur = [Pass_dur, 0];
               tlink = [tlink, 0];            
            end
            tlink(j) = i*Ts - Ts/2;		%[s]
            Pass_dur(j) = Pass_dur(j) + Ts;	%2 half-samples are added
            loop = 1;
         end
         Pass_dur(j) = Pass_dur(j) + Ts;
      else
         loop = 0;
      end
   end
   Link_dur = Pass_dur;
end

for i=1:length(Pass_dur)
    if (Pass_dur(i) == 0)
%       Pass_dur(i) = NaN;
       tlink(i) = NaN;
    end
end
 
