function dx_phys_dot = cartpole_plant(t, x_phys, u, p, disturbances_func)

x   = x_phys(1);
dx  = x_phys(2);
th  = x_phys(3);
dth = x_phys(4);

N_approx_real = max(0, (-dth.*p.beta_m.*sin(th) + p.l.*(p.M.*p.g + p.m.*(-dth.^2.*p.l + p.g.*cos(th)).*cos(th)))./p.l);

dist_val = disturbances_func(t);
F_ext   = dist_val(1);
M_ext = dist_val(2);

M_mat = [p.M + p.m.*sin(th).^2 0; p.m.*cos(th) p.l.*p.m];
F_vec = [F_ext - N_approx_real.*(p.mu_c + (-p.mu_c + p.mu_s).*exp(-dx.^2./p.v_s.^2)).*tanh(dx.*p.k) + dth.*p.beta_m.*cos(th)./p.l - dx.*p.beta_M + p.m.*(dth.^2.*p.l - p.g.*cos(th)).*sin(th) + u; M_ext./p.l - dth.*p.beta_m./p.l + p.g.*p.m.*sin(th)];

q_ddot = M_mat \ F_vec;

ddx  = q_ddot(1);
ddth = q_ddot(2);

dx_phys_dot = [dx; ddx; dth; ddth];
end