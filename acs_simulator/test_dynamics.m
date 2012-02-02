% Test dynamics model
% STG, Jan 2012, GPLv3
close all
clear all

timestep = 1;  % seconds
%timestep = 10;  % seconds
%timestep = 60;  % seconds

% Inertia matrix
I = [1 0 0;0 1 0;0 0 1];
%I = [10 0 0;0 10 0;0 0 10];

% rates at t0
%omega = [0 0 0];
%omega = [0.001 0 0];
%omega = [-0.001 0 0];
%omega = [0 0.001 0];
%omega = [0 -0.001 0];
%omega = [0 0 0.001];
%omega = [0 0 -0.001];
omega = [0.002 0.0005 -0.001];

% initial attitude
A = [1 0 0;0 1 0;0 0 1];

% forces
tau_coil = [0 0 0];
tau_jets = [0 0 0];
tau_ext = [0 0 0];

% approx 1 orbit, 100 minutes
%duration = 60*100; % seconds 
duration = 500; % 5 minutes 

index = 1;
for t=0:timestep:duration


	if((t>0) && (t<60))	
		%tau_ext = [0 0 0];
		%tau_ext = [0.0001 0 0];
		%tau_ext = [-0.0001 0 0];
		%tau_ext = [0 0.0001 0];
		%tau_ext = [0 -0.0001 0];
		%tau_ext = [0 0 0.0001];
		%tau_ext = [0 0 -0.0001];
		tau_ext = [0.0001 0.00005 -0.0001];
	else
		tau_ext = [0 0 0];
	endif


[roll,pitch,yaw,A,omega] = ... 
	Dynamics_RK4(timestep,I,omega,tau_coil,tau_jets,tau_ext,A);

roll_deg(index) = rad2deg(roll);
pitch_deg(index) = rad2deg(pitch);
yaw_deg(index) = rad2deg(yaw);

omega_save(index,:) = rad2deg(omega);

xtime(index) = t;

index = index + 1;

end

figure(1);clf
hold on
plot(xtime,roll_deg,'k*')
plot(xtime,pitch_deg,'b*')
plot(xtime,yaw_deg,'g*')
legend('roll','pitch','yaw');
xlabel('time (seconds)')
ylabel('degrees')
title('attitude vs time')

figure(2);clf
hold on
plot(xtime,omega_save(:,1),'k*')
plot(xtime,omega_save(:,2),'b*')
plot(xtime,omega_save(:,3),'g*')
legend('roll rate','pitch rate','yaw rate');
xlabel('time (seconds)')
ylabel('degrees/sec')
title('attitude rates vs time')

