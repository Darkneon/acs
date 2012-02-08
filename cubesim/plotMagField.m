%  plotMagField	Plot Earth's magnetic L-shell curves
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 13 July 02
%
%  [] = plotMagField(incl, a, e)
%	incl : 	orbit inclisation [rad](array) 
%	a :	semi-major axis [km] (scalar)
%	e :	eccentricity (scalar)
%
%  See also: getMagVec

function [] = plotMagField(incl, a, e)

% load constants
global Re

rp = a*(1-e);				%[km] periapsis
ra = a*(1+e);				%[km] apoapsis
n = 15;

%---------------------------------------
% Magnetic Field
%---------------------------------------
o = linspace(-pi, pi, 200);

%L-Shell
lp = linspace(-pi, pi, 200);
LR = (linspace(Re, 5*a, n))';
R = LR*(cos(lp)).^2;
xR = LR*(cos(lp)).^3;
yR = LR*((cos(lp)).^2.*sin(lp));

%Van Allen Radiation Belts
LVA1 = (linspace(Re+1000, 2.3*Re, n))';
LVA2 = (linspace(4*Re, 6*Re, n))';
RLVA1 = LVA1*(cos(lp)).^2;
xRVA1 = LVA1*(cos(lp)).^3;
yRVA1 = LVA1*((cos(lp)).^2.*sin(lp));
RLVA2 = LVA2*(cos(lp)).^2;
xRVA2 = LVA2*(cos(lp)).^3;
yRVA2 = LVA2*((cos(lp)).^2.*sin(lp));
for i=1:n
   for j=1:length(lp)
      if RLVA1(i,j) < Re+1000
         xRVA1(i,j) = NaN;
         yRVA1(i,j) = NaN;
      end
      if RLVA2(i,j) < Re+1000
         xRVA2(i,j) = NaN;
         yRVA2(i,j) = NaN;
      end
   end
end

% String definition
str_incl = sprintf('i = %0.4g°',incl*180/pi);
str_rp = sprintf('rp = %0.5gkm',rp);
str_ra = sprintf('ra = %0.5gkm',ra);


F0 = figure('Position', [350 450 900 500], 'Tag', 'PlotMagField');
A0 = axes('Parent', F0, 'Tag', 'PlotMagField');
title('Magnetic L-Shell Field Lines - Earth Cross-Section at Magnetic Poles', 'Parent', A0, 'Tag', 'PlotMagField')
xlabel('Distance [km]', 'Parent', A0, 'Tag', 'PlotMagField')
ylabel('Distance [km]', 'Parent', A0, 'Tag', 'PlotMagField')
axis('equal')
axis([-3*a 3*a -2*a 2*a  ])
grid
hold on

plot(Re.*cos(o), Re.*sin(o), 'k-', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotMagField');	%Earth surface
plot([-ra*cos(incl); ra*cos(incl)],[-ra*sin(incl); ra*sin(incl)],'b-',...
   'linewidth',2, 'Parent', A0, 'Tag', 'PlotMagField') %Orbit
plot(xR(1,:),yR(1,:),'g-', 'Parent', A0, 'Tag', 'PlotMagField');			% L-Shell
plot(xRVA1(1,:),yRVA1(1,:),'r-', 'Parent', A0, 'Tag', 'PlotMagField'); 			%Van Allen Radiation Belts 

H=legend('Earth Surface', ['Orbit ',str_incl,', ',str_ra,', ',str_rp],...
   'L-Shell  R = L cos²(i)', 'Van Allen Radiation Belts');
set(H, 'Tag', 'PlotMagField');
%set(get(H, 'Children'), 'Tag', 'PlotMagField'); !bad: overwrite the labels in the legend!

for (i=1:1:n)
   plot(xR(i,:),yR(i,:),'g-', 'Parent', A0, 'Tag', 'PlotMagField');				% L-Shell
   plot(xRVA1(i,:),yRVA1(i,:),'r-', 'Parent', A0, 'Tag', 'PlotMagField');
   plot(xRVA2(i,:),yRVA2(i,:),'r-', 'Parent', A0, 'Tag', 'PlotMagField');
end


plot(Re.*cos(o), Re.*sin(o), 'k-', 'linewidth',3, 'Parent', A0, 'Tag', 'PlotMagField');
for i=1:length(o)
    if (o(i) > incl) & (o(i) < pi-incl);
        o(i) = NaN;
    elseif (o(i) > incl-pi) & (o(i) < -incl);
        o(i) = NaN;
    end
end
plot([-ra*cos(incl); ra*cos(incl)],[-ra*sin(incl); ra*sin(incl)],'b-',...
   'linewidth',3, 'Parent', A0, 'Tag', 'PlotMagField')
plot(ra.*cos(o), ra.*sin(o), 'b-', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotMagField');
plot(rp.*cos(o), rp.*sin(o), 'b-', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotMagField');
plot([ra*cos(incl); -ra*cos(incl)],[-ra*sin(incl); ra*sin(incl)],'b-',...
   'linewidth',3, 'Parent', A0, 'Tag', 'PlotMagField')

plot([0,0],[0,1e+4], 'b-', 'linewidth',4, 'Parent', A0, 'Tag', 'PlotMagField')
plot([0,0],[0,-1e+4], 'r-', 'linewidth',4, 'Parent', A0, 'Tag', 'PlotMagField')
text(-500, 1.1e+4, 'S', 'color', [0 0 1], 'fontweight', 'bold', 'Parent', A0, 'Tag', 'PlotMagField')
text(-500, -1.1e+4, 'N', 'color', [1 0 0], 'fontweight', 'bold', 'Parent', A0, 'Tag', 'PlotMagField')
