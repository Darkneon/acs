%  image2vector	Converts binary image array to XY vector array
%
%	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 29 Dec 02
%
%  vector = image2vector(array)
%	array :	binary image array [mxn]
%	vector:	XY vector array
%
%  See also: imread

function vector = image2vector(array)

[m,n] = size(array);

vector = zeros(2, m*n);		% allocate memory space

k=1;
for i=1:m
    for j=1:n
        if array(i,j)==0
            vector(1,k)=i;	%
            vector(2,k)=j;	%
            k = k+1;	
        end
    end
end
vector = vector( : , 1:(k-1) );	%trim vector to last value stored
    
