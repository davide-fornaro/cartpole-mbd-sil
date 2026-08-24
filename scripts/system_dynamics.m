function dz = system_dynamics(t, z, p, K_aug, L, U_MAX, disturbances_func, ref_func)
% z = [x_phys; dx_phys; th_phys; dth_phys; x_hat; dx_hat; th_hat; dth_hat; x_i]

x_phys = z(1:4); % [x, dx, th, dth]
x_hat  = z(5:8); % [x_hat, dx_hat, th_hat, dth_hat]
x_i    = z(9);

x  = x_phys(1);
th = x_phys(3);
y_meas = [x; th];

[u, dx_hat_dot, dx_i_dot] = controller(t, x_hat, x_i, y_meas, p, K_aug, L, U_MAX, ref_func);

dx_phys_dot = cartpole_plant(t, x_phys, u, p, disturbances_func);

dz = [dx_phys_dot; dx_hat_dot; dx_i_dot];
end