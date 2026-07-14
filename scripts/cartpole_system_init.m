% =========================================================================
% NONLINEAR SIMULATION: CART-POLE SYSTEM WITH LQR + LUENBERGER OBSERVER
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
p.mu_c   = 0.05;    % [-] Kinetic (Coulomb) friction coefficient
p.mu_s   = 0.08;    % [-] Static friction coefficient (Stiction)
p.v_s    = 0.05;    % [m/s] Stribeck velocity threshold (transitions static->kinetic)
p.k      = 100;     % [1/m] Stribeck sharpness factor (higher = sharper transition)
p.ff_compensation = 0.8;

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

%% Control & Estimation Design (LQR + Observer)
% LQR Controller
Q = diag([100, 1, 500, 10]);
R = 1;
K = lqr(A, B_u, Q, R);

% Luenberger Observer (Kalman Filter is recommended for real hardware)
observer_poles = [-15, -16, -17, -18];
L = place(A', C', observer_poles)';

% LQE (Kalman Filter)
% Q_n = diag([10, 100, 10, 100]);
% R_n = diag([1e-5, 1e-5]);

% G = eye(4);
% [L, P, E] = lqe(A, G, C, Q_n, R_n);