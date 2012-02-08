%  getEphemeris - Compute S/C ephemeris from Orbit elements
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 18 July 03
%
%  [XYZsc_geo, S] = getEphemeris(Arg1, Cb, t)
%	Arg1 :	AMSAT Orbit Elements in string
%		or structure format
%	Cb :	Ballistic coefficient = m/(Cd.A) [kg/m^2]
%	t :	Time row-vector [s]
%	XYZsc_geo : S/C position in Earth's coordinates [km]
%	S : 	Structure containing orbit properties
%		using SI units (km, kg, sec, rad)
%
%³  See also: getElements

function [XYZsc_geo, S] = getEphemeris(Arg1, Cb, t)

global GMe Re Te VERNAL_EQX J2

N = length(t);
tf = t(N);
Ts = t(2)-t(1);

if isstr(Arg1)|isstruct(Arg1)
   if isstr(Arg1)
      ElementStruct = getElements(Arg1);
   elseif isstruct(Arg1)
      ElementStruct = Arg1;
   end
   
   %=== Orbit Elements ===
   %Satellite = getfield(ElementStruct, 'Satellite');             	%[]
   %CalatogNo = getfield(ElementStruct, 'Catalog_number');         	%[]
   EpochStr = getfield(ElementStruct, 'Epoch_time');              	%[yearday]
   EpochYear = sscanf(EpochStr(1:2),'%f');                   		%[year] from 2000
   EpochDay = sscanf(EpochStr(3:length(EpochStr)),'%f');          	%[day]
   t_VE = (EpochDay-VERNAL_EQX(EpochYear)) * Te + t;			%[sec] time from vernal equinox
   %ElementSet = sscanf(getfield(ElementStruct, 'Element_set'),'%f');	%
   incl = sscanf(getfield(ElementStruct, 'Inclination'),'%f')*pi/180;	%[rad]
   RA = sscanf(getfield(ElementStruct, 'RA_of_node'),'%f')*pi/180;   	%[rad]
   e = sscanf(getfield(ElementStruct, 'Eccentricity'),'%f');      	%
   ArgP = sscanf(getfield(ElementStruct, 'Arg_of_perigee'),'%f')*pi/180;%[rad]
   M = sscanf(getfield(ElementStruct, 'Mean_anomaly'),'%f')*pi/180;     %[rad]
   n = sscanf(getfield(ElementStruct, 'Mean_motion'),'%f');      	%[rev/day]
   %DecayRate = sscanf(getfield(ElementStruct, 'Decay_rate'),'%f');	%[rev/day^2]
   %EpochRev = sscanf(getfield(ElementStruct, 'Epoch_rev'),'%f');	%
   
   %=== Orbit Properties ===
   a = (sqrt(GMe)/(2*pi)*Te/n)^(2/3);	%[km] semi-major axis
   rp = a*(1-e);                 	%[km] perigee
   ra = a*(1+e);                 	%[km] apogee
   c = a - rp;                      	%[km] orbit center distance from focii
   b = sqrt(a^2-c^2);               	%[km] semi-minor axis
   hp = rp - Re;                     	%[km] altitude at perigee

   Psc = 2*pi*sqrt(a^3/GMe);         	%[s/rev] orbital period
   wsc = sqrt(GMe/a^3);             	%[rad/s] mean angular velocity

   Va = sqrt(GMe*(2/rp-1/a));        	%[km/s] velocity at apogee
   Vp = sqrt(GMe*(2/ra-1/a));        	%[km/s] velocity at perigee
   Vmean = sqrt(GMe/a);             	%[km/s] mean circular velocity
   Vesc =  sqrt(2*GMe/a);           	%[km/s] mean escape velocity
   
   
   %=== Secular rates of change from astral bodies and geopotential function (SMAD p.142) ===
   Draan_sun = -0.00338*cos(incl)/n*(pi/180)/Te;			%[rad/s]
   Draan_moon = -0.00154*cos(incl)/n*(pi/180)/Te;			%[rad/s]
   Draan_geo = -1.5*wsc*J2*(Re/a)^2*cos(incl)/(1-e^2)^2;		%[rad/s]
   Draan = Draan_sun + Draan_moon + Draan_geo;
   
   Dargp_sun = 0.00169*(4-5*(sin(incl))^2)/n*(pi/180)/Te;		%[rad/s]
   Dargp_moon = 0.00077*(4-5*(sin(incl))^2)/n*(pi/180)/Te;		%[rad/s]
   Dargp_geo = 0.75*wsc*J2*(Re/a)^2*(4-5*(sin(incl))^2)/(1-e^2)^2;	%[rad/s]
   Dargp = Dargp_sun + Dargp_moon + Dargp_geo;
   
   
   %=== Solar radiation perturbations ===
   %aR = -4.5e-6*(1+Rsc)*Asc/msc            	%[m/s2]
   
   
   %=== Atmospheric drag for low eccentricity orbits(SMAD p.145)===
   % calculations take into account matching of SI units
   % since astral dimensions are in [km] while density uses [km/m^3]
   [rhoa, ScaleHt] = getAtmDensity(hp);		%[kg/m³], [km]
   adrag = -0.5*rhoa/Cb*(Vp*1000)^2;		%[m/s^2] acceleration from drag
   DecayRate = adrag/(2*pi*a*1000)*(Te^2);  	%[rev/day^2] acceleration
   
   Da_rev = -2*pi/Cb*rhoa*(a*1000)^2/1000;	      	%[km/rev] variation of semi-major axis
   DPsc_rev = -6*pi^2/Cb*rhoa*(a*1000)^2/(Vmean*1000);	%[s/rev/rev] variation of orbit period
   DVmean_rev = pi/Cb*rhoa*(a*1000)*(Vmean*1000)/1000;	%[km/s/rev] variation of velocity
   De_rev = 0;                          		%[] variation of eccentricity
   
   Life = -ScaleHt/Da_rev;              	%[rev] estimated orbit lifetime
   
   S = struct('TimeVE',t_VE(1),'a',a,'rp',rp,'ra',ra,'c',c,'b',b,'hp',hp,'P',Psc,'w',wsc,...
      'Va',Va,'Vp',Vp,'Vmean',Vmean,'Vesc',Vesc,'rhoa',rhoa,'adrag',adrag,'Life',Life,...
      'Draan_sun',Draan_sun,'Draan_moon',Draan_moon,'Draan_geo',Draan_geo);
   
   
   %=== Time-variant orbit properties ===   
   a_t = a + Da_rev/Psc.*t;   		%[km] semi-major axis
   ra_t = a_t*(1+e);         		%[km] perigee
   rp_t = a_t*(1-e);         		%[km] apogee
   c_t = a_t - rp_t;         		%[km] orbit center distance from focii
   b_t = sqrt(a_t.^2-c_t.^2);		%[km] semi-minor axis
   hp_t = rp_t - Re;         		%[km] altitude at perigee
   
   w_t = (GMe./(a_t.^3)).^(1/2);	%[rad/s] mean angular velocity
   P_t = 2*pi./w_t;             	%[s/rev] orbital period
   n_t = Te./P_t;              		%[rev/day] mean motion
   
   Va_t = sqrt(GMe*(2./rp_t-1./a_t));	%[km/s] velocity at apogee
   Vp_t = sqrt(GMe*(2./ra_t-1./a_t));	%[km/s] velocity at perigee
   Vmean_t = sqrt(GMe./a_t);         	%[km/s] mean circular velocity
   Vesc_t =  sqrt(2*GMe./a_t);       	%[km/s] mean escape velocity
   
      
   %=== Time-variant secular rates of change from astral bodies and geopotential function ===
   Draan_sun = -0.00338*cos(incl)*(pi/180)/Te./n_t;			%[rad/s]
   Draan_moon = -0.00154*cos(incl)*(pi/180)/Te./n_t;			%[rad/s]
   Draan_geo = -1.5.*w_t.*J2.*(Re./a_t/n).^2.*cos(incl)/(1-e^2)^2;	%[rad/s]
   
   %Computation of RA taking into acound the variation of secular rates
   RA_t = zeros(1,N);
   RA_t(1) = RA;
   for i = 2:N
      RA_t(i) = RA_t(i-1) + Draan_sun(i)*Ts + Draan_moon(i)*Ts + Draan_geo(i)*Ts;	%[rad]
   end
     
   Dargp_sun = 0.00169*(4-5*(sin(incl))^2)*(pi/180)/Te./n_t;		%[rad/s]
   Dargp_moon = 0.00077*(4-5*(sin(incl))^2)*(pi/180)/Te./n_t;		%[rad/s]
   Dargp_geo = 0.75.*w_t.*J2.*(Re./a_t).^2.*(4-5*(sin(incl))^2)/(1-e^2)^2;	%[rad/s]
   
   %Computation of ArgP taking into acound the variation of secular rates
   ArgP_t = zeros(1,N);
   ArgP_t(1) = ArgP;
   for i = 2:N
      ArgP_t(i) = ArgP_t(i-1) + Dargp_sun(i)*Ts + Dargp_moon(i)*Ts + Dargp_geo(i)*Ts;	%[rad]
   end
   
   %===================================================
   % Orbit calculations (WERTZ p.134, SMAD p.140)
   %===================================================
   %nu_t = M_t + 2*e.*sin(M_t) + 1.25*e^2.*sin(2.*M_t);	%[rad] True anomaly (approximation for small eccentricity)
   %E_t = acos(e+cos(nu_t))./(1+e.*cos(nu_t));	%[rad] Eccentric anomaly
   %E_t = 2.*atan(tan(nu_t./2)./sqrt((1+e)./(1-e)));%[rad] Eccentric anomaly 
   %R_t = (rp_t.*(1+e))./(1 + e.*cos(nu_t));	%[km] Distance from focus 
      
   M_t = M + w_t.*t;					%[rad] Mean anomaly
   
   %=== Newton iterative algorythm to solve M = E - e.sin(E) (WERTZ p.134)===
   count = 0;
   E_t = zeros(1,N);					%[rad] Eccentric anomaly
   for i=1:N
      Ek = 0;
      Ek2= min(max(mod(M_t(i),2*pi), 0.2),2*pi-0.2);	%first approx trimmed to avoid singularity and
      while (abs(Ek2-Ek)> 1e-12)			%optimized for fast convergence over 0~2pi range
         Ek = Ek2;
         Ek2 = Ek + (M_t(i) + e*sin(Ek)- Ek)/(1 - e*cos(Ek));
         count = count+1;
      end
      E_t(i) = Ek2;
   end
   for i=1:N						%quicksolve the remaining singularities
      if isnan(E_t(i))
         E_t(i) = (E_t(i-1)+E_t(i+1))/2;
         E_t(i) = E_t(i) + (M_t(i) + e*sin(E_t(i))- E_t(i))/(1 - e*cos(E_t(i)));
         E_t(i) = E_t(i) + (M_t(i) + e*sin(E_t(i))- E_t(i))/(1 - e*cos(E_t(i)));
         E_t(i) = E_t(i) + (M_t(i) + e*sin(E_t(i))- E_t(i))/(1 - e*cos(E_t(i)));         
      end
   end
   %count   

   nu_t = 2.*atan(tan(E_t./2)./sqrt((1-e)./(1+e)));	%[rad] True anomaly (exact)
   R_t = a_t.*(1 - e.*cos(E_t));			%[km] Distance from focus (based on approx)
   V_t = sqrt(GMe.*(2./R_t-1./a_t));			%[km/s] Local velocity
   u_t = ArgP_t + nu_t;					%[rad] Angle to ascending node
   
   %=== Ephemeris ===
   LLsc = getTrace(incl, RA_t, u_t ,t_VE);		%[rad]
   LATsc = LLsc(2,:);                                  	%[rad]
   LONsc = LLsc(1,:);                                  	%[rad] 
   LONsc = mod(LONsc+pi,2*pi*ones(size(LONsc)))-pi;    	%[rad] normalisation @ ±Pi
   LLsc = [LONsc; LATsc];                             	%[rad] S/C position in geo referential
   
   XYZsc_geo = LL2XYZ(LLsc, R_t);			%[km] S/C position in geo referential
   
else
   fprintf('\nERROR: first argument must be AMSAT orbit element structure or string\n')
   XYZgeo = [];
end




