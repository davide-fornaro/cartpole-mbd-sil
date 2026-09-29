function dx_phys_dot = cartpole_plant(t, x_phys, u, p, disturbances_func)
    x   = x_phys(1);
    dx  = x_phys(2);
    th  = x_phys(3);
    dth = x_phys(4);

    dist_val = disturbances_func(t);
    F_ext   = dist_val(1);
    M_ext = dist_val(2);

    q_ddot = compute_accelerations(dx, th, dth, u, p, F_ext, M_ext);
    dx_phys_dot = [dx; q_ddot(1); dth; q_ddot(2)];
end