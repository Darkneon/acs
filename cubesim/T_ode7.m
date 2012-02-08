%  T_ode7	ODE file for Temperature simulation
%		using 7rd order cube model
%
%  	Jean-Francois Levesque, MS
%	jflev@yahoo.ca
%	Last Update: 18 July 03
%
%  T_dot = T_ode7(ti, T, option, t, Qin, Rsc, Eclipse_state,...
%   Qearth, Qsun, Qalbedo, VFearth, PVpacking, wsc, Asc, csc, wallsc, msc)
%	ti : 	instantaneous time [s]
%	T : 	insantaneous temperature array [K]
%	option : (not used)
%	t :	time vector [s]
%	Qin :	Heat generated [W]
%	Rsc :	S/C distance array [km] (array)
%	Eclipse_state : array (1 = eclipse, 0 = sun input)
%	Qearth : Radiation IR on faces [W] (array)
%	Qsun : Radiation Sun on faces [W] (array)
%	Qalbedo : Radiation Sun on faces [W] (array)
%	VFearth : Earth View factor of faces (array)
%	wsc : 	S/C Revolution [rad/s](scalar)
%	Asc : 	S/C face area  [m�] (scalar)
%	csc : 	S/C size [m] (scalar)
%	wallsc : S/C wall thickness [m] (scalar)
%	msc : 	S/C mass [kg] (scalar)
%	T_dot :	temperature rate array [K/s]
%
%   Notes:
%   Gin is almost negligible, producing a temperature change of << 0.01K/s
%   Thus, the internal equilibrium temperature is realy close to the faces temperature,
%   but with a slow dynamic
%   The conduction process is much faster than the radiaton process, such that the
%   temperature is almost even on every faces of the cube
%   This Model is not numerically accrate/stable since Qcond varies a lot
%
%   See Also: T_ode2


function T_dot = T_ode7(ti, T, option, t, Qin, Rsc, ...%Eclipse_state,...
   Qearth, Qsun, Qalbedo, VFearth, PVpacking, wsc, Asc, csc, wallsc, msc)

global Boltzmann

global Re
global Qtot7
global Qtot7b

T0 = 4;		% [K] Deep Space temperature

% --- Materials ---
k_alu = 167.7;	% [W/m-K] conduction constant
cp_alu = 920;	% [J/kg-K] specific heat
k_cu = 389;	% [W/m-K] conduction constant
cp_cu = 390;	% [J/kg-K] specific heat
k_si = 0.151;	% [W/m-K] conduction constant
cp_si = 712;	% [J/kg-K] specific heat
k_kapton = 0;	% [W/m-K] conduction constant
cp_kapton = 1006;	% [J/kg-K] specific heat
k_nylon = 0.0332;	% [W/m-K] conduction constant
cp_nylon = 1680;	% [J/kg-K] specific heat

% --- Radiation Coefficients ---
emi_alu = 0.0346;
abs_alu = 0.379;
emi_bk = 0.874;
abs_bk = 0.975;
emi_cell = 0.825;
abs_cell = 0.805;
emi_gold = 0.023;
abs_gold = 0.299;

% --- Computing Data ---
i = find(t>=ti);
if length(Qin) > 1
   Qin = Qin(i(1));
end
%Gs = Gsolar(i(1));
Rsc = Rsc(i(1));
%Eclipse_state = Eclipse_state(i(1));
Qsun = Qsun(:,i(1));
Qearth = Qearth(:,i(1));
Qalbedo = Qalbedo(:,i(1));
VFearth = VFearth(:,i(1));
% --- Temperature ---
Twall = T(1:6,1);
Tin = T(7);

% --- Calculations ---
rho_earth = asin(Re/Rsc); 				%[rad]

m_wall = 0.1*msc;
cp_wall = cp_alu;
k_wall = k_alu;
emi_wall = PVpacking * emi_cell + (1-PVpacking) * emi_alu;
abs_wall = PVpacking * abs_cell + (1-PVpacking) * abs_alu;

m_in = msc - 6*m_wall;
cp_in = cp_alu;
emi_in = emi_bk;

% --- Space Radiation ---
q_radspace = Boltzmann .* emi_wall .* (Twall.^4 - T0^4);
Qspace = q_radspace .* (1-VFearth) * Asc;	%[W]

% --- Internal Radiation (assume view factor of 1 for each face) ---
q_radin = Boltzmann * emi_in .* emi_wall .* (Tin^4 - Twall.^4);		%[W/m�]
Qradin = q_radin * Asc;					%[W]
%Qradin = [0;0;0;0;0;0];

% --- Conduction Between Faces ---
Qcond = zeros(6,6);
for j=1:6
   for k=1:6
      Qcond(j:k) = k_wall*0.085*csc/(0.085) * (Twall(j) - Twall(k));	%[W/m�]
   end
end

% Differential Equations
Twall_dot = ( Qsun + Qalbedo + Qearth + Qradin - Qspace...
   -(Qcond(:,1) + Qcond(:,2) + Qcond(:,3) + Qcond(:,4) + Qcond(:,5) + Qcond(:,6)) )...
   ./ (m_wall * cp_wall);
Tin_dot = ( Qin - (Qradin(1) + Qradin(2) + Qradin(3) + Qradin(4) + Qradin(5) + Qradin(6)) ) / (m_in* cp_in);

T_dot = [Twall_dot; Tin_dot];

Qtot7(:,i(1)) = [( Qsun + Qalbedo + Qearth + Qradin - Qspace -(Qcond(:,1) + Qcond(:,2) + Qcond(:,3) + Qcond(:,4) + Qcond(:,5) + Qcond(:,6)) )
                ( Qin - (Qradin(1) + Qradin(2) + Qradin(3) + Qradin(4) + Qradin(5) + Qradin(6)) )];
             
Qtot7b(:,i(1)) = [( Qsun + Qalbedo + Qearth + Qradin - Qspace )
                ( Qin - (Qradin(1) + Qradin(2) + Qradin(3) + Qradin(4) + Qradin(5) + Qradin(6)) )];
             
             
