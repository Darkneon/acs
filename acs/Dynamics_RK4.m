function [roll,pitch,yaw,A_new,omega_new] = Dynamics_RK4(timestep,I,omega,tau_coil,tau_jets,tau_ext,A)

% function [roll,pitch,yaw,A_new,omega_new] = Dynamics_RK4(timestep,I,omega,tau_coil,tau_jets,tau_ext,A)
%
% Calculates roll, pitch and yaw and angular rates over time interval
% reaction wheels not included!
%
% Inputs: timestep = time to propagate over
%	  I = Inertia Matrix, 3x3
%	  omega = current angular velocity vector, 3x1
%	  tau_coils = torque on S/C due to coils, 3x1
%	  tau_jets = torque on S/C due to thruster jets, 3x1
%	  tau_ext = total external torque on S/C, 3x1
%         A = attitude matrix
%
% Outputs: roll = spacecraft roll angle (radians)
%	   pitch = spacecraft pitch angle (radians)
%	   yaw = spacecraft yaw angle (radians)
%	   A_new = new attitude matrix
%	   omega_new = angular velocity vector
%
% Author: Scott Gleason, 2012
% License: GPLv3
%
% Begin state update 
%
% torque applied over entire interval
%
%

% Times of RK stages
t(1) = 0.0;
t(2) = timestep*0.5;
t(3) = timestep*0.5;
t(4) = timestep;

% total torques over this interval
tau = tau_coil + tau_jets + tau_ext;

% convert A matrix to quaternian
[q1] = A2q(A);

% calculate first stage, S1
omega_S1 = calc_omega(omega,tau,t(1),I);
% calculate q rates
[qdot1] = qdot(omega,q1);  % rate at start
k1 = qdot1.*timestep;   % delta value using slope at start

% calculate second stage, S2
omega_S2 = calc_omega(omega,tau,t(2),I);   % rate at halfway point, with torques
% calculate q rates
[qdot2] = qdot(omega_S2,q1+k1/2);  % rate at halfway point, y+k1
k2 = qdot2.*timestep;   % delta value using slope at halfway point

% calculate second stage, S3
omega_S3 = calc_omega(omega,tau,t(3),I);  % rate at halfway point, with torques
% calculate q rates
[qdot3] = qdot(omega_S3,q1+k2/2);  % rate at halfway point, y+k2
k3 = qdot3.*timestep;   % delta value using slope at halfway point

% calculate second stage, S4
omega_S4 = calc_omega(omega,tau,t(4),I);   % rate at end, with torques
% calculate q rates
[qdot4] = qdot(omega_S4,q1+k3);  % rate at end, y+k3
k4 = qdot4.*timestep;   % delta value using slope at end
	
% Integrate to end of timestep

delta = (1.0/6.0).*k1 + (1.0/3.0).*k2 + (1.0/3.0).*k3 + (1.0/6.0).*k4;

q_final = q1 + delta;

omega_new = ((1.0/6.0).*omega_S1 + (1.0/3.0).*omega_S2 + ... 
	(1.0/3.0).*omega_S3 + (1.0/6.0).*omega_S4).*timestep;

% normalize final quaternian
q_final_norm = q_final./norm(q_final);

% convert q to Attitude matrix
A_new = q2A(q_final_norm);

% convery new Attitude matrix to roll pitch and yaw 
[roll,pitch,yaw] = A2rpy_321(A_new);


