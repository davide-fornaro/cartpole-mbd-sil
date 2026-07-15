cartpole_system_init;

%% Simulation Initialization

% global state vector for ODE solver: [Physical States; Estimated States; Integral State]
z0 = [x0_phys; x0_hat; x0_i];

% stop simulation if the cart hits the track limits
ode_opts = odeset('Events', @(t,z) track_limit(t, z, p.TRACK_LIMIT), 'MaxStep', 0.01);

disturbances_func = @(t) [
    5 * (t >= 15.0 & t <= 15.1);         % F_cart [N]
    0.5 * (t >= 30.0 & t <= 30.05)         % M_ext [Nm]
];

ref_func = @(t) 0.3 * (t >= 40) - 0.4 * (t >= 70.0);

%% Numerical Integration
t_span = [0 90];
disp('Starting nonlinear ODE integration...');
[t_out, Z] = ode15s(@(t, z) system_dynamics(t, z, p, K, L, A, B_u, C, p.U_MAX, disturbances_func, ref_func), t_span, z0, ode_opts);

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
    
    u_history(i) = max(min(u_req, p.U_MAX), -p.U_MAX);
end

%% Data Visualization
figure('Name', 'Nonlinear Simulation Results', 'Position', [100 100 900 800]);

subplot(3,1,1);
plot(t_out, x_real, 'b', 'LineWidth', 1.5); hold on;
plot([0 t_out(end)], [p.TRACK_LIMIT p.TRACK_LIMIT], 'r--', 'LineWidth', 1.5);
plot([0 t_out(end)], [-p.TRACK_LIMIT -p.TRACK_LIMIT], 'r--', 'LineWidth', 1.5);
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
plot([0 t_out(end)], [p.U_MAX p.U_MAX], 'r--');
plot([0 t_out(end)], [-p.U_MAX -p.U_MAX], 'r--');
grid on; ylabel('Force [N]'); xlabel('Time [s]');
title('Actuator Effort & Saturation Limits');

if max(abs(x_real)) >= p.TRACK_LIMIT - 1e-3
    warning('CRITICAL: The cart hit the track limits. System failed.');
else
    disp('SUCCESS: System stabilized without hardware collisions.');
end