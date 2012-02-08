%   plotHysteresis - Sketches 3 magnetic hysteresis models
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 23 Sept 03
%
%  plotHysteresis(Br, Bs, Hc)
%	Br	: Remanence [A/m^2]
%	Bs 	: Saturation induction [A/m^2]
%	Hc 	: Coercivity [Tesla]

function plotHysteresis(Br, Bs, Hc)

[a1 b1 c1]=solve('c=-Bs','a*exp(0)+c=-Br','a*exp(Hc/b)+c=0','a','b','c');
[a2 b2 c2]=solve('c=-Bs','a*exp(0)+c=-Br*1.2','a*exp(Hc/b)+c=0','a','b','c');
[a3 b3 c3]=solve('c=Bs','a*exp(0)+c=Br','a*exp(-Hc/b)+c=0','a','b','c');
[a4 b4 c4]=solve('c=Bs','a*exp(0)+c=Br*1.2','a*exp(-Hc/b)+c=0','a','b','c');
[a5 b5 c5]=solve('c=Bs','a*exp(0)+c=Br*1.1','a*exp(-Hc/b)+c=0','a','b','c');

a1 = eval(a1);b1 = eval(b1);c1 = eval(c1);
a2 = eval(a2);b2 = eval(b2);c2 = eval(c2);
a3 = eval(a3);b3 = eval(b3);c3 = eval(c3);
a4 = eval(a4);b4 = eval(b4);c4 = eval(c4);
a5 = eval(a5);b5 = eval(b5);c5 = eval(c5);

x=linspace(-4*Hc,4*Hc,200);
y1=a1*exp((x)./b1)+c1;
y2=a2*exp((x+2*Hc)./b2)+c2;
y3=a3*exp((x)./b3)+c3;
y4=a4*exp((x-2*Hc)./b4)+c4;
y5=a5*exp((x-Hc)./b5)+c5;

for i = 1:length(x)
   if y1(i) > 0
      y1(i) = NaN;
   end
   if y2(i) > 0
      y2(i) = NaN;
   end
   if y3(i) < 0
      y3(i) = NaN;
   end
   if y4(i) < 0
      y4(i) = NaN;
   end
   if y5(i) < 0
      y5(i) = NaN;
   end
end

F0 = figure('Tag', 'PlotHyst');
A0 = axes('Parent', F0, 'Tag', 'PlotHyst');
axis([min(x) max(x) -Bs Bs]*1.1)
hold on
xlabel('Field Strength H [Tesla]', 'Parent', A0, 'Tag', 'PlotHyst')
ylabel('Magnetic Induction B [A/m^2]', 'Parent', A0, 'Tag', 'PlotHyst')
title('Magnetic Hysteresis Models', 'Parent', A0, 'Tag', 'PlotHyst')
plot([x x],[y1 y4],'r', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([-Hc-(Hc/Br*Bs) -Hc+(Hc/Br*Bs)], [-Bs Bs], 'g', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([-3*Hc Hc], [-Bs Bs], 'b', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([-Hc max(x)], [Bs Bs], 'color', [0.7 0.7 0.7], 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')

H=legend('Hysteresis curve','Model 1','Model 2','Model 3', 4);
set(H, 'Tag', 'PlotHyst');

plot([0 0], 1.1*[-Bs Bs], 'k', 'linewidth',1, 'Parent', A0, 'Tag', 'PlotHyst')
plot(1.2*[min(x) max(x)], [0 0], 'k', 'linewidth',1, 'Parent', A0, 'Tag', 'PlotHyst')
plot([-Hc -Hc], 1.05*[-Bs Bs], 'k--', 'linewidth',1, 'Parent', A0, 'Tag', 'PlotHyst')
plot([Hc Hc], 1.05*[-Bs Bs], 'k--', 'linewidth',1, 'Parent', A0, 'Tag', 'PlotHyst')

plot([x x],[y2 y3],'r', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([x],[y5],'r--', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')

plot([Hc-(Hc/Br*Bs) Hc+(Hc/Br*Bs)], [-Bs Bs], 'g', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([-Hc+(Hc/Br*Bs) max(x)], [Bs Bs], 'g', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([min(x) -Hc-(Hc/Br*Bs)], [-Bs -Bs], 'g', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')

plot([-Hc 3*Hc], [-Bs Bs], 'b', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([Hc max(x)], [Bs Bs], 'b', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([min(x) -Hc], [-Bs -Bs], 'b', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')

plot([-Hc -Hc], [-Bs Bs], 'color', [0.7 0.7 0.7], 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([Hc Hc], [-Bs Bs], 'color', [0.7 0.7 0.7], 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([min(x) Hc], [-Bs -Bs], 'color', [0.7 0.7 0.7], 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')

text(Hc*1.1, -Br/10,'Hc', 'Parent', A0, 'Tag', 'PlotHyst')
text(-Hc*0.95, -Br/10,'-Hc', 'Parent', A0, 'Tag', 'PlotHyst')
text(-Hc/3, Br,'Br', 'Parent', A0, 'Tag', 'PlotHyst')
text(-Hc/2.5, -Br,'-Br', 'Parent', A0, 'Tag', 'PlotHyst')
text(-Hc/3, Bs*0.95,'Bs', 'Parent', A0, 'Tag', 'PlotHyst')
text(-Hc/2.5, -Bs*0.95,'-Bs', 'Parent', A0, 'Tag', 'PlotHyst')
plot([-Hc Hc]/15,[Br Br],'k', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([-Hc Hc]/15,-[Br Br],'k', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([-Hc Hc]/15,[Bs Bs],'k', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')
plot([-Hc Hc]/15,-[Bs Bs],'k', 'linewidth',2, 'Parent', A0, 'Tag', 'PlotHyst')

