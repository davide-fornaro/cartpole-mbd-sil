function dz = system_dynamics(t, z, p, K, L, A, B, C, U_MAX, disturbances_func, ref_func)
% z = [x_phys; dx_phys; th_phys; dth_phys; x_hat; dx_hat; th_hat; dth_hat; x_i]
x_phys = z(1:4); % [x, dx, th, dth]
x_hat  = z(5:8); % [x_hat, dx_hat, th_hat, dth_hat]
x_i    = z(9);

x   = x_phys(1);
dx  = x_phys(2);
th  = x_phys(3);
dth = x_phys(4);

y_meas = [x; th];

r = ref_func(t);

th_hat  = x_hat(3);
dth_hat = x_hat(4);
N_approx_real = max(0, (-dth.*p.beta_m.*sin(th) + p.l.*(p.M.*p.g + p.m.*(-dth.^2.*p.l + p.g.*cos(th)).*cos(th)))./p.l);
N_approx_est  = max(0, (-p.beta_m*sin(th_hat)*dth_hat + p.l*(p.M*p.g + p.m*(p.g*cos(th_hat) - p.l*dth_hat^2)*cos(th_hat)))/p.l);
F_coulomb_est = p.mu_c * N_approx_est;

% friction compensation
u_ff = p.ff_compensation * F_coulomb_est * tanh(p.k * x_hat(2));

u_lqi = -K(1:4) * x_hat - K(5) * x_i;
u_req = u_lqi + u_ff;

u = max(min(u_req, U_MAX), -U_MAX); % u real

dist_val = disturbances_func(t);
F_ext   = dist_val(1);
M_ext = dist_val(2);

M_mat = [p.M + p.m.*sin(th).^2 0; p.m.*cos(th) p.l.*p.m];
F_vec = [F_ext - N_approx_real.*(p.mu_c + (-p.mu_c + p.mu_s).*exp(-dx.^2./p.v_s.^2)).*tanh(dx.*p.k) + dth.*p.beta_m.*cos(th)./p.l - dx.*p.beta_M + p.m.*(dth.^2.*p.l - p.g.*cos(th)).*sin(th) + u; M_ext./p.l - dth.*p.beta_m./p.l + p.g.*p.m.*sin(th)];

q_ddot = M_mat \ F_vec;

ddx  = q_ddot(1);
ddth = q_ddot(2);

dx_phys_dot = [dx; ddx; dth; ddth];

dx_hat_dot = A * x_hat + B * (u - u_ff) + L * (y_meas - C * x_hat);

dx_i = y_meas(1) - r;

if (u_req >= U_MAX && dx_i > 0) || (u_req <= -U_MAX && dx_i < 0)
    dx_i = 0;
end

dz = [dx_phys_dot; dx_hat_dot; dx_i];
end