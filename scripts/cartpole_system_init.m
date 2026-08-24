% NONLINEAR SIMULATION: CART-POLE SYSTEM WITH LQI + CONSTANT-GAIN NONLINEAR OBSERVER

clearvars; close all; clc;

%% Physical Parameters & Hardware Constraints
p.M      = 0.8;     % [kg] Cart mass
p.m      = 0.15;     % [kg] Pendulum mass
p.l      = 0.3;     % [m] Length from pivot to Center of Mass
p.g      = 9.81;    % [m/s^2] Gravity acceleration

% Viscous Friction (Proportional to velocity)
p.beta_M = 2;     % [N*s/m] Cart viscous friction
p.beta_m = 0.002;   % [N*m*s/rad] Joint viscous friction

% Nonlinear Friction (Stribeck & Coulomb)
p.mu_c   = 0.025;   % [-] Kinetic (Coulomb) friction coefficient
p.mu_s   = 0.04;    % [-] Static friction coefficient (Stiction)
p.v_s    = 0.05;    % [m/s] Stribeck velocity threshold (transitions static->kinetic)
p.k      = 5;       % [1/m] Stribeck sharpness factor (higher = sharper transition)
p.ff_compensation = 0.85; % [-] Feedforward friction compensation factor

% Hardware Limits
p.U_MAX       = 20;   % [N] Maximum force from the DC motor
p.TRACK_LIMIT = 0.6;  % [m] Maximum physical rail distance from center (+/- 0.6m)

% Disturbance Configuration
p.dist.max_impulse_force = 10.0; % [N] Maximum force of a discrete hit
p.dist.max_impulse_torque = 2.0; % [Nm] Maximum torque hit
p.dist.continuous_force_std = 0.001; % [N] Baseline continuous friction noise
p.dist.continuous_torque_std = 0.0005; % [Nm] Baseline continuous aero noise

%% Linearized Model Matrices (Evaluated at the upright equilibrium)

A = [0 1 0 0; 0 -p.beta_M./p.M -p.g.*p.m./p.M p.beta_m./(p.M.*p.l); 0 0 0 1; 0 p.beta_M./(p.M.*p.l) p.g.*(p.M + p.m)./(p.M.*p.l) p.beta_m.*(-p.M - p.m)./(p.M.*p.l.^2.*p.m)];
B_full = [0 0 0; 1./p.M 1./p.M 0; 0 0 0; -1./(p.M.*p.l) -1./(p.M.*p.l) 1./(p.l.^2.*p.m)];
C = [1, 0, 0, 0;   % Sensor 1: Cart position encoder
     0, 0, 1, 0];  % Sensor 2: Pendulum angle encoder
D = [0 0 0; 0 0 0];

B_u = B_full(:, 1); % Control input: Horizontal force applied to the cart
% B_d = B_full(:, 2:3); % Disturbance inputs: [External horizontal force; External torque at the pivot]

%% Control & Estimation Design (LQI + Observer)
C_i = [1, 0, 0, 0];
A_aug = [A, zeros(4,1); C_i, 0];
B_aug = [B_u; 0];

% LQI Optimal Weights Synthesis (Bryson's Rule)
epsilon = 1e-6;

% Max tolerances: [x [m], dx [m/s], th [rad], dth [rad/s], xi [m*s]]
max_devs = max([p.TRACK_LIMIT * 0.9, 1.0, deg2rad(15), deg2rad(100), 0.5], epsilon);
Q_aug = diag(1 ./ (max_devs.^2));
R_aug = 1 / (max(p.U_MAX, epsilon)^2);

K_aug = lqr(A_aug, B_aug, Q_aug, R_aug);

% LQE (Kalman Filter)

% Sensor Measurement Noise Covariance (R_n)
% Assuming uniform quantization error distribution: variance = (resolution^2) / 12
res_x  = 1e-3;           % [m] Cart linear encoder resolution (e.g., 1 mm)
res_th = deg2rad(0.1);   % [rad] Pendulum rotary encoder resolution (e.g., 0.1 deg)

var_sens_x  = (res_x.^2) / 12;
var_sens_th = (res_th.^2) / 12;

R_n = diag([var_sens_x, var_sens_th]);

% Process Noise Covariance (Q_n) & Disturbance Mapping (G_noise)
% Process noise enters as physical forces/torques, affecting only accelerations.
var_force  = (p.dist.max_impulse_force / 3.0)^2;      % [N^2] Variance of unmodeled horizontal forces (e.g., friction errors)
var_torque = (p.dist.max_impulse_torque / 3.0)^2;     % % [(N*m)^2] Variance of unmodeled torques at the pivot joint  

Q_n = diag([var_force, var_torque]);

% Disturbance mapping matrix: maps [force_noise; torque_noise] to states [x, dx, th, dth]
G_noise = [0, 0;
    1, 0;
    0, 0;
    0, 1];

[L, P, E] = lqe(A, G_noise, C, Q_n, R_n);

% Continuous Observer Matrices (for ODE simulation)
% Pre-computing state-space matrices to avoid O(N) operations in the solver loop
A_obs = A - L * C;
B_obs = [B_u, L];
C_obs = eye(4);
D_obs = zeros(4, 3);

%% Swing-Up Controller Parameters (Energy-Based)
p.swing.k_E = 25.0;
p.swing.k_p = 40.0;
p.swing.k_d = 12.0;
p.swing.th_thresh = deg2rad(40);
p.swing.kick_amp  = 3.0; % [N] Forza asimmetrica per rompere l'equilibrio a theta = pi

%% Simulation Initialization
t_span = [0 180]; % [s] Total simulation time

angle_init = pi;
pos_init   = 0.3;

% Initial Physical State (Reality): Pendulum downwards (pi rad)
x0_phys = [pos_init; 0; angle_init; 0]; % [ position (m); velocity (m/s); angle (rad); angular velocity (rad/s) ]

% Initial Software State (Microcontroller)
x0_hat = [0; 0; 0; 0];

% Initial Integral State: Starts at zero  (no accumulated error)
x0_i = 0;

%% DISCRETE MODEL & CONTROL ADDITIONS (Kalman / DLQE)
p.Ts = 0.005; % [s] Microcontroller sampling time

% Exact ZOH Discretization of the Plant
sys_c = ss(A, B_u, C, D(:,1));
sys_d = c2d(sys_c, p.Ts, 'zoh');

Ad   = sys_d.A;
Bd_u = sys_d.B;
Cd   = sys_d.C;

% Discrete LQI Controller Design
Ad_aug = [Ad, zeros(4,1); p.Ts*Cd(1,:), 1];
Bd_aug = [Bd_u; 0];

Qd_aug = Q_aug;
Rd_aug = R_aug;

[Kd_aug, ~, ~] = dlqr(Ad_aug, Bd_aug, Qd_aug, Rd_aug);
Kd_x = Kd_aug(1:4);
Kd_i = Kd_aug(5);

% Discrete Kalman Filter (DLQE)
Q_n_discrete = Q_n * p.Ts;
R_n_discrete = R_n / p.Ts;
G_discrete   = G_noise;

[Ld, ~, ~] = dlqe(Ad, G_discrete, Cd, Q_n_discrete, R_n_discrete);

% Discrete Observer Matrices
Ad_obs = Ad - Ld * Cd;
Bd_obs = [Bd_u, Ld];
Cd_obs = eye(4);
Dd_obs = zeros(4, 3);

%% Disturbances & Reference Signals

% disturbances_func = @(t) [
%     5 * (t >= 15.0 & t <= 15.1);         % F_cart [N]
%     0.5 * (t >= 30.0 & t <= 30.05)       % M_ext [Nm]
% ];

ref_func = @(t) 0.3 * (t >= 40) - 0.4 * (t >= 70.0);

% Stochastic Disturbance Generation

% rng('shuffle')
rng(42);

% Continuous Band-Limited Noise (e.g., wind gusts, uneven friction)
t_noise = t_span(1):0.5:t_span(2);

F_noise_base = p.dist.continuous_force_std * randn(1, length(t_noise));
M_noise_base = p.dist.continuous_torque_std * randn(1, length(t_noise));

pp_F = spline(t_noise, F_noise_base);
pp_M = spline(t_noise, M_noise_base);

% Discrete Random Impulses (Simulating unexpected physical hits)
N_impulses = 12; % Number of random hits during the simulation
time_needed_to_swing_up = 10; % [s] Time needed to swing up

% Adjusted random window to prevent out-of-bounds generation and overlap
impulse_times = t_span(1) + time_needed_to_swing_up + (t_span(2) - t_span(1) - time_needed_to_swing_up - 0.1) * rand(1, N_impulses);
impulse_forces = (rand(1, N_impulses) * 2 - 1) * p.dist.max_impulse_force;

% Applied boolean masking to continuous noise during swing-up phase
disturbances_func = @(t) [
    (t >= time_needed_to_swing_up) * ppval(pp_F, t) + sum(impulse_forces .* (t >= impulse_times & t <= impulse_times + 0.1));
    (t >= time_needed_to_swing_up) * ppval(pp_M, t)
    ];