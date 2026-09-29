function [u, dx_hat_dot, dx_i, u_ff] = controller(t, x_hat, x_i, y_meas, p, K_aug, L, ref_func)
    r = ref_func(t);

    x_hat_cart = x_hat(1);
    dx_hat     = x_hat(2);
    th_hat     = x_hat(3);
    dth_hat    = x_hat(4);

    % Continuous Trigonometric Normalization
    th_err_lqi = atan2(sin(th_hat), cos(th_hat));
    
    x_lqi = [x_hat_cart - r; dx_hat; th_err_lqi; dth_hat];

    % Feedforward Friction Compensation
    N_approx_est  = compute_normal_force_approx(th_hat, dth_hat, p);
    F_coulomb_est = p.mu_c * N_approx_est;
    
    u_ff = p.ff_compensation * F_coulomb_est * tanh(p.k * dx_hat);

    % LQI Stabilization
    u_lqi = -K_aug(1:4) * x_lqi - K_aug(5) * x_i;

    % Energy-Based Swing-Up
    E_kinetic   = 0.5 * p.m * (p.l * dth_hat)^2;
    E_potential = p.m * p.g * p.l * cos(th_hat); 
    E_current   = E_kinetic + E_potential;
    E_target    = p.m * p.g * p.l;
    E_error     = E_current - E_target;

    x_penalty = p.swing.k_p * ((x_hat_cart - r) / p.TRACK_LIMIT)^3;

    % Symmetry Breaker (C-infinity Starter Kick)
    is_down_and_still = exp(-0.5 * (dth_hat/0.2)^2) * exp(-0.5 * ((cos(th_hat) + 1)/0.05)^2);
    u_kick = p.swing.kick_amp * is_down_and_still;

    u_swing = p.swing.k_E * E_error * dth_hat * cos(th_hat) ...
              - x_penalty ...
              - p.swing.k_d * dx_hat ...
              + u_kick;

    % C-infinity Convex Blending
    sigma = p.swing.th_thresh / 1.5; 
    weight_LQI = exp(-0.5 * (abs(th_err_lqi) / sigma)^6);
    u_req = weight_LQI * u_lqi + (1 - weight_LQI) * u_swing + u_ff;

    % Anti-windup Clamp
    dx_i = weight_LQI * (y_meas(1) - r);
    if (u_req >= p.U_MAX && dx_i > 0) || (u_req <= -p.U_MAX && dx_i < 0)
        dx_i = 0;
    end

    % Actuator Saturation
    u = max(min(u_req, p.U_MAX), -p.U_MAX);

    % Constant-Gain Nonlinear Observer (Prediction Step)
    F_ext_hat = 0;
    M_ext_hat = 0;
    
    q_ddot_est = compute_accelerations(dx_hat, th_hat, dth_hat, u, p, F_ext_hat, M_ext_hat);
    
    f_x_hat = [dx_hat; q_ddot_est(1); dth_hat; q_ddot_est(2)];
    
    dx_hat_dot = f_x_hat + L * [y_meas(1) - x_hat_cart; atan2(sin(y_meas(2) - th_hat), cos(y_meas(2) - th_hat))];
end