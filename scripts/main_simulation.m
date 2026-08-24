cartpole_system_init;

%% Simulation Initialization

% global state vector for ODE solver: [Physical States; Estimated States; Integral State]
z0 = [x0_phys; x0_hat; x0_i];

ode_opts = odeset('Events', @(t,z) track_limit(t, z, p.TRACK_LIMIT), 'MaxStep', 0.01);

%% Numerical Integration
disp('Starting nonlinear ODE integration...');
[t_out, Z] = ode15s(@(t, z) system_dynamics(t, z, p, K_aug, L, p.U_MAX, disturbances_func, ref_func), t_span, z0, ode_opts);

%% Data Extraction & Analysis
x_real     = Z(:, 1);
x_hat      = Z(:, 5);
theta_real = Z(:, 3);
theta_hat  = Z(:, 7);

% Memory Pre-allocation (O(1) allocation overhead)
u_history    = zeros(length(t_out), 1);
u_ff_history = zeros(length(t_out), 1); % Extract feedforward for plotting

for i = 1:length(t_out)
    t_i      = t_out(i);
    x_est    = Z(i, 5:8)';
    x_i_est  = Z(i, 9);
    y_meas_i = [Z(i, 1); Z(i, 3)];
    
    [u_history(i), ~, ~, u_ff_history(i)] = controller(t_i, x_est, x_i_est, y_meas_i, p, K_aug, L, p.U_MAX, ref_func);
end

%% Data Visualization
figure('Name', 'Nonlinear Simulation Results', 'Position', [100 100 900 800]);

subplot(3,1,1);
plot(t_out, x_real, 'b', 'LineWidth', 1.5); hold on;
plot(t_out, x_hat, 'r--', 'LineWidth', 1.2);
plot([0 t_out(end)], [p.TRACK_LIMIT p.TRACK_LIMIT], 'k--', 'LineWidth', 1.2);
plot([0 t_out(end)], [-p.TRACK_LIMIT -p.TRACK_LIMIT], 'k--', 'LineWidth', 1.2);
grid on; ylabel('Position [m]');
title('Cart Spatial Tracking & Collision Verification');
legend('Cart Position (Real)', 'Cart Position (Estimated)', 'Hardware Limits', 'Location', 'Best');

subplot(3,1,2);
plot(t_out, theta_real*180/pi, 'k', 'LineWidth', 1.5); hold on;
plot(t_out, theta_hat*180/pi, 'g--', 'LineWidth', 1.2);
grid on; ylabel('Angle [deg]');
title('Nonlinear Angular Stabilization');
legend('\theta True', '\theta Estimated');

subplot(3,1,3);
plot(t_out, u_history, 'm', 'LineWidth', 1.5); hold on;
plot(t_out, u_ff_history, 'c--', 'LineWidth', 1.2); % Overlay feedforward signal
plot([0 t_out(end)], [p.U_MAX p.U_MAX], 'r--');
plot([0 t_out(end)], [-p.U_MAX -p.U_MAX], 'r--');
grid on; ylabel('Force [N]'); xlabel('Time [s]');
title('Actuator Effort & Saturation Limits');
legend('Total Effort (u)', 'Feedforward (u_{ff})', 'U_{MAX}');

if max(abs(x_real)) >= p.TRACK_LIMIT - 1e-3
    warning('CRITICAL: The cart hit the track limits. System failed.');
else
    disp('SUCCESS: System stabilized without hardware collisions.');
end