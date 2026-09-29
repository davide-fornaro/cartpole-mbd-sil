function plot_cartpole_results(t_out, x_real, x_hat, theta_real, theta_hat, u_history, u_ff_history, p, is_discrete)
    if is_discrete
        fig_name = 'Simulation Results (Discrete)';
    else
        fig_name = 'Simulation Results (Continuous)';
    end
    
    figure('Name', fig_name, 'Position', [100 100 900 800]);

    subplot(3,1,1);
    plot(t_out, x_real, 'b', 'LineWidth', 1.5); hold on;
    plot(t_out, x_hat, 'r--', 'LineWidth', 1.2);
    plot([t_out(1) t_out(end)], [p.TRACK_LIMIT p.TRACK_LIMIT], 'k--', 'LineWidth', 1.2);
    plot([t_out(1) t_out(end)], [-p.TRACK_LIMIT -p.TRACK_LIMIT], 'k--', 'LineWidth', 1.2);
    grid on; ylabel('Position [m]');
    title('Cart Spatial Tracking');
    legend('Cart Position (Real)', 'Cart Position (Estimated)', 'Hardware Limits', 'Location', 'Best');

    subplot(3,1,2);
    plot(t_out, theta_real*180/pi, 'k', 'LineWidth', 1.5); hold on;
    plot(t_out, theta_hat*180/pi, 'g--', 'LineWidth', 1.2);
    grid on; ylabel('Angle [deg]');
    title('Angular Stabilization');
    legend('\theta True', '\theta Estimated', 'Location', 'Best');

    subplot(3,1,3);
    if is_discrete
        stairs(t_out, u_history, 'm', 'LineWidth', 1.5); hold on;
        stairs(t_out, u_ff_history, 'c--', 'LineWidth', 1.2);
    else
        plot(t_out, u_history, 'm', 'LineWidth', 1.5); hold on;
        plot(t_out, u_ff_history, 'c--', 'LineWidth', 1.2);
    end
    plot([t_out(1) t_out(end)], [p.U_MAX p.U_MAX], 'r--');
    plot([t_out(1) t_out(end)], [-p.U_MAX -p.U_MAX], 'r--');
    grid on; ylabel('Force [N]'); xlabel('Time [s]');
    title('Actuator Effort');
    legend('Total Effort (u)', 'Feedforward (u_{ff})', 'U_{MAX}', 'Location', 'Best');
    
    if max(abs(x_real)) >= p.TRACK_LIMIT - 1e-3
        warning('CRITICAL: The cart hit the track limits. System failed.');
    else
        disp('SUCCESS: System stabilized without hardware collisions.');
    end
end