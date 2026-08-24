cartpole_system_init;

%% Simulation Initialization (Discrete Paradigm)

t_sim = t_span(1) : p.Ts : t_span(2);
N_steps = length(t_sim);

% O(1) Memory Pre-allocation to avoid dynamic resizing overhead
X_real_log = zeros(N_steps, 4);
X_hat_log  = zeros(N_steps, 4);
U_log      = zeros(N_steps, 1);
U_ff_log   = zeros(N_steps, 1); % Extract discrete feedforward 

% Initial states
x_phys = x0_phys;
x_hat  = x0_hat;
x_i    = x0_i;

disp('Starting Discrete integration...');

for k = 1:N_steps
    
    t_k = t_sim(k);
    
    % Sensor acquisition (Sampling)
    y_meas = [x_phys(1); x_phys(3)];
    ref    = ref_func(t_k);
    
    % Microcontroller execution (ZOH generation)
    [u_real, x_hat_next, x_i_next, u_ff_real] = controller_discrete_step(x_hat, x_i, y_meas, p, Kd_aug, Ld, p.U_MAX, ref);
    
    % Data logging
    X_real_log(k, :) = x_phys';
    X_hat_log(k, :)  = x_hat';
    U_log(k)         = u_real;
    U_ff_log(k)      = u_ff_real;
    
    if abs(x_phys(1)) >= p.TRACK_LIMIT
        disp(['Track limit hit at t = ', num2str(t_k), 's.']);
    end
    
    % Physical plant progression (Between k and k+1)
    if k < N_steps
        [~, Z_ode] = ode15s(@(t, z) cartpole_plant(t, z, u_real, p, disturbances_func), ...
                           [t_k, t_k + p.Ts], x_phys);
        
        x_phys = Z_ode(end, :)';
    end
    
    % State shifting
    x_hat = x_hat_next;
    x_i   = x_i_next;
end

%% Data Extraction
t_out      = t_sim(1:N_steps)';
x_real     = X_real_log(1:N_steps, 1);
x_hat      = X_hat_log(1:N_steps, 1);
theta_real = X_real_log(1:N_steps, 3);
theta_hat  = X_hat_log(1:N_steps, 3);
u_history  = U_log(1:N_steps);
u_ff_hist  = U_ff_log(1:N_steps);

%% Data Visualization
figure('Name', 'Discrete Simulation Results', 'Position', [150 150 900 800]);

subplot(3,1,1);
plot(t_out, x_real, 'b', 'LineWidth', 1.5); hold on;
plot(t_out, x_hat, 'r--', 'LineWidth', 1.2);
plot([0 t_out(end)], [p.TRACK_LIMIT p.TRACK_LIMIT], 'k--', 'LineWidth', 1.2);
plot([0 t_out(end)], [-p.TRACK_LIMIT -p.TRACK_LIMIT], 'k--', 'LineWidth', 1.2);
grid on; ylabel('Position [m]');
title('Cart Spatial Tracking (ZOH Control)');
legend('Cart Position (Real)', 'Cart Position (Estimated)', 'Hardware Limits', 'Location', 'Best');

subplot(3,1,2);
plot(t_out, theta_real*180/pi, 'k', 'LineWidth', 1.5); hold on;
plot(t_out, theta_hat*180/pi, 'g--', 'LineWidth', 1.2);
grid on; ylabel('Angle [deg]');
title('Discrete Angular Stabilization');
legend('\theta True', '\theta Estimated');

subplot(3,1,3);
stairs(t_out, u_history, 'm', 'LineWidth', 1.5); hold on; % 'stairs' reflects exact MCU ZOH behavior
stairs(t_out, u_ff_hist, 'c--', 'LineWidth', 1.2); % Discrete feedforward tracking
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