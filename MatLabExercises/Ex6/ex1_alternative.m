% This is a alternative solution, which uses dsolve to solve the differential equation 
clear; clc; close all 

%% Initializing constants
mu_0=4*pi*1e-7;
h_bar=1.05457*1e-34;
gamma_H=2.675e8;
gamma_N=-2.716e7;
B0=14.1;
r12=1.1e-10;
delta_iso=0;
delta_iso_H=0;
tau_c=1e-10;

%% elements of the relaxation matrix
omega_H = -gamma_H*B0;omega_N=-gamma_N*B0;
delta_HN = -2*mu_0*gamma_H*gamma_N*h_bar/4/pi/(r12^3);
% intitial conditions
S_z0 = -1; % initital condition for S spin after inversion
S_zeq = 1; % thermal equilibrium
I_z0 = S_zeq*omega_H/omega_N; % I spin in thermal equilibrium = I_z(eq)
I_zeq = I_z0;
% elements of the relaxation matrix
J_tau = @(omega) 2/5*tau_c/(1+(omega*tau_c).^2); % isotropic rotational tumbling 
G_Iz = (delta_HN/4)^2*(J_tau(omega_H-omega_N)+3*J_tau(omega_H)+ ... 
    6*J_tau(omega_H+omega_N))+3/4*(omega_H*delta_iso_H)^2*J_tau(omega_H); % matrix elements of relaxation matrix
G_Sz = (delta_HN/4)^2*(J_tau(omega_H-omega_N)+3*J_tau(omega_N)+ ...
    6*J_tau(omega_H+omega_N))+3/4*(omega_N*delta_iso)^2*J_tau(omega_N);
G_IzSz = (delta_HN/4)^2*(-J_tau(omega_H-omega_N)+6*J_tau(omega_H+omega_N));

%% Coupled differential equations
syms Iz(t) Sz(t);
% integrated solver for coupled differential equations in matlab
G = -[G_Iz G_IzSz; G_IzSz G_Sz]; % relaxation matrix
Y = [Iz; Sz]; % time dependence of spins
B = [I_zeq; S_zeq]; % equilibrium value
C = Y(0) == [I_z0; S_z0]; % initial condition
eqn = diff(Y) == G*(Y-B);
[IzSol(t), SzSol(t)] = dsolve(eqn, C);
% plot S_z
figure(1)
fplot(SzSol,[0 15])
hold on;
title('')
xlabel('t / s')
ylabel('<S_Z(t)>')
ylim([-1 1])

%% Monoexponential decay 
T1 = 1/G_Sz; % T1 for monoexponential decay of S_z
t=(0:0.02:15); % generate time axis
Sz_mono = 1-2*exp(-t/T1); % evolution of S_z magnetization during inversion recovery
figure(1)
plot(t,Sz_mono,'r')
legend('coupled system', 'T_1 decay')


