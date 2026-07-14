function dz = system_dynamics(t, z, p, K, L, A, B, C, U_MAX, disturbances_func)
    x_phys = z(1:4); % [x, dx, th, dth]
    x_hat  = z(5:8); % [x_hat, dx_hat, th_hat, dth_hat]

    x   = x_phys(1);
    dx  = x_phys(2);
    th  = x_phys(3);
    dth = x_phys(4);

    y_meas = [x; th];

    N_approx = (-p.beta_m*sin(th)*dth + p.l*(p.M*p.g + p.m*(p.g*cos(th) - p.l*dth^2)*cos(th)))/p.l;
    F_coulomb_est = p.mu_c * N_approx;

    % friction compensation
    u_ff = p.ff_compensation * F_coulomb_est * tanh(p.k * x_hat(2));

    u_lqr = -K * x_hat;
    u_req = u_lqr + u_ff;

    u = max(min(u_req, U_MAX), -U_MAX); % u real

    dist_val = disturbances_func(t);
    F_ext   = dist_val(1);
    M_ext = dist_val(2);

    M_mat = [
        p.M + p.m*sin(th)^2, 0;
        p.m*cos(th), p.l*p.m
    ];
    F_vec = [
        u + F_ext - N_approx*(p.mu_c + (-p.mu_c + p.mu_s)*exp(-dx^2/p.v_s^2))*tanh(p.k*dx) - p.beta_M*dx + p.beta_m*cos(th)*dth/p.l + p.m*(-p.g*cos(th) + p.l*dth^2)*sin(th);
        M_ext/p.l - p.beta_m*dth/p.l + p.g*p.m*sin(th)
    ];

    q_ddot = M_mat \ F_vec;

    ddx  = q_ddot(1);
    ddth = q_ddot(2);

    dx_phys_dot = [dx; ddx; dth; ddth];

    dx_hat_dot = A * x_hat + B * u + L * (y_meas - C * x_hat);

    dz = [dx_phys_dot; dx_hat_dot];
end