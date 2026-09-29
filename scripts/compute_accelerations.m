function q_ddot = compute_accelerations(dx, th, dth, u, p, F_ext, M_ext)
    sin_th = sin(th);
    cos_th = cos(th);

    N_approx = compute_normal_force_approx(th, dth, p);
    
    M_mat = [p.M + p.m.*(sin_th.^2), 0; 
             p.m.*cos_th,            p.l.*p.m];
    
    F_vec = [F_ext - N_approx.*(p.mu_c + (-p.mu_c + p.mu_s).*exp(-dx.^2./p.v_s.^2)).*tanh(dx.*p.k) + dth.*p.beta_m.*cos_th./p.l - dx.*p.beta_M + p.m.*(dth.^2.*p.l - p.g.*cos_th).*sin_th + u; 
             M_ext./p.l - dth.*p.beta_m./p.l + p.g.*p.m.*sin_th];
    
    q_ddot = M_mat \ F_vec;
end