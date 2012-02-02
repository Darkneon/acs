function omega = calc_omega(omega_old,tau,timestep,I); 

% function omega = calc_omega(tau_coil,tau_jets,tau_exttimestep,I); 
%
% Rotational equation, calculated at each RK time step
% No Wheels!
%
% Inputs: omega_old = starting angular velocity
%	  tau = total torque on S/C, 3x1
%	  timestep = duration of above torques 
%	  I = Inertia Matrix, 3x3
%
% Outputs: omega = vector of angular rate at end of timestep
%
% Author: Scott Gleason, 2012
% License: GPLv3
%

% total angular momentum of the spacecraft: L = I*w
L = I*omega_old';

%  omega_dot = tau/I
omega_dot = tau*inv(I);

% new omega
omega = omega_old + omega_dot*timestep;

