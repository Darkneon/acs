%  getElements - Retrieve orbit elements from AMSAT string format
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 14 July 03
%
%  ElementStruct = getElements(ElementStr)
%
%***********************************************
%*  AMSAT Orbit Element Format
%*  -----------------------------------
%*  Satellite: AO-07 
%*  Catalog number: 07530
%*  Epoch time:  02303.68051080
%*  Element set:  0148
%*  Inclination:  101.7801 deg
%*  RA of node:   347.4291 deg
%*  Eccentricity:  0.0012182
%*  Arg of perigee:   292.1694 deg
%*  Mean anomaly: 067.8083 deg
%*  Mean motion:  12.53560994 rev/day
%*  Decay rate:   -2.9e-07 rev/day^2
%*  Epoch rev:   27929
%*  Checksum:  294
%* 
%   String example:
%   ---------------
%ElementStr = ['Satellite: AO-07',10,'Catalog number: 07530',10,
%'Epoch time:  02303.68051080',10,'Element set:  0148',10,
%'Inclination:  101.7801 deg',10,'RA of node:   347.4291 deg',10,
%'Eccentricity:  0.0012182',10,'Arg of perigee:   292.1694 deg',10,
%'Mean anomaly: 067.8083 deg',10,'Mean motion:  12.53560994 rev/day',
%10,'Decay rate:   -2.9e-07 rev/day^2',10,
%'Epoch rev:   27929',10,'Checksum:  294'];
%***********************************************

function ElementStruct = getElements(ElementStr)

CHAR_LF = 10;
ElementStr2 = [CHAR_LF, ElementStr, CHAR_LF];	%Add LineFeed at begining and end of string
LF = findstr(ElementStr2, CHAR_LF);
Colon = findstr(ElementStr2, ':');

ElementArray = cell(length(LF)-1,2);
for i= 1:(length(LF)-1)
   Field = ElementStr2(LF(i)+1:Colon(i)-1);
   Space = findstr(Field, ' ');
   for j= 1:length(Space)
      Field(Space(j)) = ('_');			% replace all spaces by underscore in Fields
   end
   ElementArray{i,1} = Field;
   
   Value = ElementStr2(Colon(i)+1:LF(i+1)-1);
   while Value(1)==(' ')
      Value = Value(2:length(Value));	% remove all spaces at begining of Values
   end
   Space = findstr(Value, ' ');
   if ~isempty(Space)
      Value = Value(1:Space(1)-1); 	% remove all characters after Value (ex. units)
   end
   ElementArray{i,2} = Value;
end

ElementStruct = struct(ElementArray{1,1},ElementArray{1,2},...
   ElementArray{2,1},ElementArray{2,2},...
   ElementArray{3,1},ElementArray{3,2},...
   ElementArray{4,1},ElementArray{4,2},...
   ElementArray{5,1},ElementArray{5,2},...
   ElementArray{6,1},ElementArray{6,2},...
   ElementArray{7,1},ElementArray{7,2},...
   ElementArray{8,1},ElementArray{8,2},...
   ElementArray{9,1},ElementArray{9,2},...
   ElementArray{10,1},ElementArray{10,2},...
   ElementArray{11,1},ElementArray{11,2},...
   ElementArray{12,1},ElementArray{12,2},...
   ElementArray{13,1},ElementArray{13,2});

