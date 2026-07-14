% =========================================================================
% NONLINEAR SIMULATION: CART-POLE SYSTEM WITH LQI + LUENBERGER OBSERVER
% =========================================================================
clearvars; close all; clc;

%% Physical Parameters & Hardware Constraints
p.M      = 1.0;     % [kg] Cart mass
p.m      = 0.5;     % [kg] Pendulum mass
p.l      = 0.5;     % [m] Length from pivot to Center of Mass
p.g      = 9.81;    % [m/s^2] Gravity acceleration

% Viscous Friction (Proportional to velocity)
p.beta_M = 0.5;     % [N*s/m] Cart viscous friction
p.beta_m = 0.005;   % [N*m*s/rad] Joint viscous friction

% Nonlinear Friction (Stribeck & Coulomb)
p.mu_c   = 0.025;    % [-] Kinetic (Coulomb) friction coefficient
p.mu_s   = 0.04;    % [-] Static friction coefficient (Stiction)
p.v_s    = 0.05;    % [m/s] Stribeck velocity threshold (transitions static->kinetic)
p.k      = 5;       % [1/m] Stribeck sharpness factor (higher = sharper transition)
p.ff_compensation = 0.90; 

% Hardware Limits
U_MAX       = 20;   % [N] Maximum force from the DC motor
TRACK_LIMIT = 0.6;  % [m] Maximum physical rail distance from center (+/- 0.6m)

%% Linearized Model Matrices (Evaluated at the upright equilibrium)
A = [
     0, 1, 0, 0;
     0, -p.beta_M/p.M, -p.g*p.m/p.M, p.beta_m/(p.M*p.l);
     0, 0, 0, 1;
     0, p.beta_M/(p.M*p.l), p.g*(p.M + p.m)/(p.M*p.l), p.beta_m*(-p.M - p.m)/(p.M*p.l^2*p.m)
     ];
B_full = [
     0, 0, 0;
     1/p.M, 1/p.M, 0;
     0, 0, 0;
     -1/(p.M*p.l), -1/(p.M*p.l), 1/(p.l^2*p.m)
     ];
C = [1, 0, 0, 0;   % Sensor 1: Cart position encoder
     0, 0, 1, 0];  % Sensor 2: Pendulum angle encoder
D = [
     0, 0, 0;
     0, 0, 0
     ];

B_u = B_full(:, 1);
B_d = B_full(:, 2:3);

%% Control & Estimation Design (LQI + Observer)
C_i = [1, 0, 0, 0];
A_aug = [A, zeros(4,1); C_i, 0];
B_aug = [B_u; 0];

Q_aug = diag([10, 1, 1000, 1, 50]);
R_aug = 1;

K = lqr(A_aug, B_aug, Q_aug, R_aug);

% LQE (Kalman Filter)
Q_n = diag([1e-2, 1e-1, 1e-2, 1e-1]);
R_n = diag([1e-3, 1e-3]);
G_noise = eye(4);
G = eye(4);

[L, P, E] = lqe(A, G, C, Q_n, R_n);

% Luenberger Observer (Kalman Filter is recommended for real hardware)
% observer_poles = [-15, -16, -17, -18];
% L = place(A', C', observer_poles)';