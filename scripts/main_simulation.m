cartpole_system_init;

%% Simulation Initialization
% Initial Physical State (Reality): Pendulum tilted by 15 degrees (0.26 rad)
x0_phys = [0; 0; 0.26; 0];

% Initial Software State (Microcontroller): Starts at zero (no knowledge)
x0_hat = [0; 0; 0.26; 0];

% Initial Integral State: Starts at zero (no accumulated error)
x0_i = 0;

% global state vector for ODE solver: [Physical States; Estimated States; Integral State]
z0 = [x0_phys; x0_hat; x0_i];

% stop simulation if the cart hits the track limits
ode_opts = odeset('Events', @(t,z) track_limit(t, z, TRACK_LIMIT), 'MaxStep', 0.01);

disturbances_func = @(t) [
    5 * (t >= 15.0 & t <= 15.1);         % F_cart [N]
    0.5 * (t >= 40.0 & t <= 40.05)         % M_ext [Nm]
];

%% Numerical Integration
t_span = [0 60];
disp('Starting nonlinear ODE integration...');
[t_out, Z] = ode15s(@(t, z) system_dynamics(t, z, p, K, L, A, B_u, C, U_MAX, disturbances_func), t_span, z0, ode_opts);

%% Data Extraction & Analysis
x_real     = Z(:, 1);
theta_real = Z(:, 3);
theta_hat  = Z(:, 7);

u_history = zeros(length(t_out), 1);
for i = 1:length(t_out)
    x_est = Z(i, 5:8)';
    
    x_i_est = Z(i, 9);
    
    th_hat_sim  = x_est(3);
    dth_hat_sim = x_est(4);
    N_est = (-p.beta_m*sin(th_hat_sim)*dth_hat_sim + p.l*(p.M*p.g + p.m*(p.g*cos(th_hat_sim) - p.l*dth_hat_sim^2)*cos(th_hat_sim)))/p.l;
    
    F_c_est = p.mu_c * N_est;
    u_ff = p.ff_compensation * F_c_est * tanh(p.k * x_est(2));
    
    u_lqi = -K(1:4) * x_est - K(5) * x_i_est;
    u_req = u_lqi + u_ff;
    
    u_history(i) = max(min(u_req, U_MAX), -U_MAX);
end

%% Data Visualization
figure('Name', 'Nonlinear Simulation Results', 'Position', [100 100 900 800]);

subplot(3,1,1);
plot(t_out, x_real, 'b', 'LineWidth', 1.5); hold on;
plot([0 t_out(end)], [TRACK_LIMIT TRACK_LIMIT], 'r--', 'LineWidth', 1.5);
plot([0 t_out(end)], [-TRACK_LIMIT -TRACK_LIMIT], 'r--', 'LineWidth', 1.5);
grid on; ylabel('Position [m]');
title('Cart Spatial Tracking & Collision Verification');
legend('Cart Position', 'Hardware Limits');

subplot(3,1,2);
plot(t_out, theta_real*180/pi, 'k', 'LineWidth', 1.5); hold on;
plot(t_out, theta_hat*180/pi, 'g--', 'LineWidth', 1.2);
grid on; ylabel('Angle [deg]');
title('Nonlinear Angular Stabilization');
legend('\theta True', '\theta Estimated');

subplot(3,1,3);
plot(t_out, u_history, 'm', 'LineWidth', 1.5); hold on;
plot([0 t_out(end)], [U_MAX U_MAX], 'r--');
plot([0 t_out(end)], [-U_MAX -U_MAX], 'r--');
grid on; ylabel('Force [N]'); xlabel('Time [s]');
title('Actuator Effort & Saturation Limits');

if max(abs(x_real)) >= TRACK_LIMIT - 1e-3
    warning('CRITICAL: The cart hit the track limits. System failed.');
else
    disp('SUCCESS: System stabilized without hardware collisions.');
end