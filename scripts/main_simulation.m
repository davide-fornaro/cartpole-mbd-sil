cartpole_system_init;

%% Simulation Initialization

% global state vector for ODE solver: [Physical States; Estimated States; Integral State]
z0 = [x0_phys; x0_hat; x0_i];

ode_opts = odeset('Events', @(t,z) track_limit(t, z, p.TRACK_LIMIT), 'MaxStep', 0.01);

%% Numerical Integration
disp('Starting nonlinear ODE integration...');
[t_out, Z] = ode15s(@(t, z) system_dynamics(t, z, p, K_aug, L, disturbances_func, ref_func), t_span, z0, ode_opts);

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
    
    [u_history(i), ~, ~, u_ff_history(i)] = controller(t_i, x_est, x_i_est, y_meas_i, p, K_aug, L, ref_func);
end

%% Data Visualization
plot_cartpole_results(t_out, x_real, x_hat, theta_real, theta_hat, u_history, u_ff_history, p, false);