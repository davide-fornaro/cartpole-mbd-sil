function [u, x_hat_next, x_i_next, u_ff] = controller_discrete_step(x_hat, x_i, y_meas, p, Kd_aug, Ld, ref)
    
    x_hat_cart = x_hat(1);
    dx_hat     = x_hat(2);
    th_hat     = x_hat(3);
    dth_hat    = x_hat(4);

    % Maps angle to [-pi, pi] for continuous LQI feedback
    th_err_lqi = atan2(sin(th_hat), cos(th_hat));
    
    % LQI Error State Vector (Corrected Reference Tracking)
    x_lqi = [x_hat_cart - ref; dx_hat; th_err_lqi; dth_hat];

    % Estimated friction computation
    N_approx_est  = compute_normal_force_approx(th_hat, dth_hat, p);
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
    u = max(min(u_req_unclamped, p.U_MAX), -p.U_MAX);

    error_i = weight_LQI * (y_meas(1) - ref);
    
    % Anti-Windup
    if (u_req_unclamped >= p.U_MAX && error_i > 0) || (u_req_unclamped <= -p.U_MAX && error_i < 0)
        dx_i = 0; % Clamp integral action
    else
        dx_i = error_i;
    end
    
    % Forward Euler integration
    x_i_next = x_i + p.Ts * dx_i;

    % Constant-Gain Nonlinear Observer (Prediction Step)
    F_ext_hat = 0;
    M_ext_hat = 0;
    q_ddot_est = compute_accelerations(dx_hat, th_hat, dth_hat, u, p, F_ext_hat, M_ext_hat);

    f_x_hat = [dx_hat; q_ddot_est(1); dth_hat; q_ddot_est(2)];

    % Discrete Constant-Gain Innovation Step
    y_hat = [x_hat_cart; th_hat];
    residual_raw = y_meas - y_hat;
    
    % Normalize angular residual
    innovation = [residual_raw(1); 
                  atan2(sin(residual_raw(2)), cos(residual_raw(2)))];

    % Forward Euler Integration for the discrete observer state
    x_hat_next = x_hat + p.Ts * f_x_hat + Ld * innovation;

end