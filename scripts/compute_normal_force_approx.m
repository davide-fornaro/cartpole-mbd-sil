function N_approx = compute_normal_force_approx(th, dth, p)
    sin_th = sin(th);
    cos_th = cos(th);

    N_approx = max(0, (-dth.*p.beta_m.*sin_th + p.l.*(p.M.*p.g + p.m.*(-dth.^2.*p.l + p.g.*cos_th).*cos_th))./p.l);
end