function [u, x_hat_next, x_i_next, u_ff] = controller_discrete_step(x_hat, x_i, y_meas, p, Kd_aug, Ld, U_MAX, ref)
    
    x_hat_cart = x_hat(1);
    dx_hat     = x_hat(2);
    th_hat     = x_hat(3);
    dth_hat    = x_hat(4);

    % Maps angle to [-pi, pi] for continuous LQI feedback
    th_err_lqi = atan2(sin(th_hat), cos(th_hat));
    
    % LQI Error State Vector (Corrected Reference Tracking)
    x_lqi = [x_hat_cart - ref; dx_hat; th_err_lqi; dth_hat];

    % Estimated friction computation
    N_approx_est  = max(0, (-p.beta_m*sin(th_hat)*dth_hat + p.l*(p.M*p.g + p.m*(p.g*cos(th_hat) - p.l*dth_hat^2)*cos(th_hat)))/p.l);
    F_coulomb_est = p.mu_c * N_approx_est;
    u_ff = p.ff_compensation * F_coulomb_est * tanh(p.k * dx_hat);

    
    % LQI Stabilization
    u_lqi = -Kd_aug(1:4) * x_lqi - Kd_aug(5) * x_i;

    % Energy-Based Swing-Up
    E_kinetic   = 0.5 * p.m * (p.l * dth_hat)^2;
    E_potential = p.m * p.g * p.l * cos(th_hat); 
    E_current   = E_kinetic + E_potential;
    E_target    = p.m * p.g * p.l;
    E_error     = E_current - E_target;

    % Virtual Soft-Wall (Cubic Penalty for hardware limits)
    x_penalty = p.swing.k_p * ((x_hat_cart - ref) / p.TRACK_LIMIT)^3;

    % Symmetry Breaker (C-infinity Starter Kick at bottom dead center)
    is_down_and_still = exp(-0.5 * (dth_hat/0.2)^2) * exp(-0.5 * ((cos(th_hat) + 1)/0.05)^2);
    u_kick = p.swing.kick_amp * is_down_and_still;

    u_swing = p.swing.k_E * E_error * dth_hat * cos(th_hat) ...
              - x_penalty ...
              - p.swing.k_d * dx_hat ...
              + u_kick;

    % C-infinity Convex Blending
    sigma = p.swing.th_thresh / 1.5; 
    weight_LQI = exp(-0.5 * (abs(th_err_lqi) / sigma)^6);
    u_req_unclamped = weight_LQI * u_lqi + (1 - weight_LQI) * u_swing + u_ff;

    % Saturation clamp
    u = max(min(u_req_unclamped, U_MAX), -U_MAX);

    error_i = weight_LQI * (y_meas(1) - ref);
    
    % Anti-Windup
    if (u_req_unclamped >= U_MAX && error_i > 0) || (u_req_unclamped <= -U_MAX && error_i < 0)
        dx_i = 0; % Clamp integral action
    else
        dx_i = error_i;
    end
    
    % Forward Euler integration
    x_i_next = x_i + p.Ts * dx_i;

    % Constant-Gain Nonlinear Observer (Prediction Step)
    
    sin_th = sin(th_hat);
    cos_th = cos(th_hat);

    M11 = p.M + p.m*sin_th^2;
    M21 = p.m*cos_th;
    M22 = p.l*p.m;

    F1 = -F_coulomb_est * tanh(dx_hat*p.k) + dth_hat*p.beta_m*cos_th/p.l - dx_hat*p.beta_M + p.m*(dth_hat^2*p.l - p.g*cos_th)*sin_th + u;
    F2 = -dth_hat*p.beta_m/p.l + p.g*p.m*sin_th;
                  
    det_M = M11 * M22;
    inv_M11 = 1 / M11;
    inv_M21 = -M21 / det_M;
    inv_M22 = 1 / M22;

    ddx_hat_est = inv_M11 * F1;
    ddth_hat_est = inv_M21 * F1 + inv_M22 * F2;

    % Non-linear state derivative
    f_x_hat = [dx_hat; ddx_hat_est; dth_hat; ddth_hat_est];

    % Discrete Constant-Gain Innovation Step
    y_hat = [x_hat_cart; th_hat];
    residual_raw = y_meas - y_hat;
    
    % Normalize angular residual
    innovation = [residual_raw(1); 
                  atan2(sin(residual_raw(2)), cos(residual_raw(2)))];

    % Forward Euler Integration for the discrete observer state
    x_hat_next = x_hat + p.Ts * f_x_hat + Ld * innovation;

end