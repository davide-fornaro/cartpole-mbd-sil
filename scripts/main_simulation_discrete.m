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
        N_steps = k;
        break;
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
t_out        = t_sim(1:N_steps)';
x_real       = X_real_log(1:N_steps, 1);
x_hat        = X_hat_log(1:N_steps, 1);
theta_real   = X_real_log(1:N_steps, 3);
theta_hat    = X_hat_log(1:N_steps, 3);
u_history    = U_log(1:N_steps);
u_ff_history = U_ff_log(1:N_steps);

%% Data Visualization
plot_cartpole_results(t_out, x_real, x_hat, theta_real, theta_hat, u_history, u_ff_history, p, true);