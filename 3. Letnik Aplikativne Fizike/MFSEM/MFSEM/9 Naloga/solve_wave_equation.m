% solve_wave_equation.m
% Solves the 1D wave equation using an explicit finite difference method (CTCS).

clear;
clc;
close all;

% --- Parameters ---
a = 1.0;         % Length of the string (m)
c_wave = 100.0;  % Wave speed (m/s) (c in problem doc)
y0_peak = 0.05;  % Peak amplitude of initial triangle (m)
x_peak_rel = 0.4; % Relative position of the peak (0.4 means at x=0.4a)

% Spatial grid
dx_val = 0.01;   % Spatial step (m) (h in problem doc)
nx = round(a / dx_val) + 1;
x = linspace(0, a, nx);

% Temporal grid
% Courant condition: C_stab = c_wave * dt / dx < 1
C_stab = 1.1; % Courant stability factor
dt_val = C_stab * dx_val / c_wave;

% Simulation time: omega1*t/pi from 0 to 2
% omega1 = pi*c_wave/a (fundamental angular frequency)
% So, (pi*c_wave/a)*t/pi = c_wave*t/a from 0 to 2
% This means t ranges from 0 to 2*a/c_wave
t_max = 2 * a / c_wave;

nt = round(t_max / dt_val) + 1;
t_vec = linspace(0, t_max, nt);
dt_val = t_vec(2) - t_vec(1); % Recalculate dt

fprintf('Spatial step dx = %.4f m, Number of spatial points nx = %d\n', dx_val, nx);
fprintf('Temporal step dt = %.4f s, Number of temporal steps nt = %d\n', dt_val, nt);
fprintf('Courant number C = c*dt/dx = %.4f\n', c_wave * dt_val / dx_val);
fprintf('Total simulation time t_max = %.3f s (corresponds to c*t/a = %.1f)\n', t_max, c_wave*t_max/a);

if c_wave * dt_val / dx_val >= 1
    warning('Courant condition C < 1 might be violated. C = %.4f', c_wave * dt_val / dx_val);
end

% --- Initial Conditions ---
u = zeros(nx, nt); % Displacement matrix u(x,t)

% Initial shape u(x,0) = f(x): triangular with peak y0_peak at x_peak_rel*a
x_peak_abs = x_peak_rel * a;
for m = 1:nx
    if x(m) <= x_peak_abs
        if x_peak_abs == 0 % Avoid division by zero if peak is at x=0
            u(m,1) = y0_peak; % Assume it's a step down if peak is at 0
        else
            u(m,1) = y0_peak * x(m) / x_peak_abs;
        end
    else
        if (a - x_peak_abs) == 0 % Avoid division by zero if peak is at x=a
             u(m,1) = y0_peak; % Assume it's a step down if peak is at a
        else
            u(m,1) = y0_peak * (a - x(m)) / (a - x_peak_abs);
        end
    end
end
% Ensure ends are exactly zero if rounding caused issues
u(1,1) = 0;
u(nx,1) = 0;


% --- Boundary Conditions ---
% Fixed ends: u(0,t) = 0, u(a,t) = 0
u(1, :) = 0;
u(nx, :) = 0;

% --- Initial Velocity Condition (for u_m^1, i.e., u(:,2)) ---
% Initial velocity is zero. Formula: u_m^1 = u_m^0 + 0.5 * (c*dt/dx)^2 * (u_{m+1}^0 - 2u_m^0 + u_{m-1}^0)
beta_sq = (c_wave * dt_val / dx_val)^2; % This is r in python script, or (kc/h)^2 in problem doc

if nt > 1 % Only compute if there's more than one time step
    for m = 2:nx-1 % Internal points
        u(m,2) = u(m,1) + 0.5 * beta_sq * (u(m+1,1) - 2*u(m,1) + u(m-1,1));
    end
    % Boundary conditions u(1,2) and u(nx,2) are already zero.
end

% --- Simulation (CTCS for n >= 2) ---
% Formula: u_m^{n+1} = 2u_m^n - u_m^{n-1} + beta_sq * (u_{m+1}^n - 2u_m^n + u_{m-1}^n)
for n = 2:nt-1 % Time loop (n is current time, computing n+1)
    for m = 2:nx-1 % Spatial loop (internal points)
        u(m, n+1) = 2*u(m,n) - u(m,n-1) + beta_sq * (u(m+1,n) - 2*u(m,n) + u(m-1,n));
    end
    % Boundary conditions u(1,n+1) and u(nx,n+1) are already zero.
end

% --- Plotting ---
figure;
hold on;
% Plot based on c*t/a values (0, 0.5, 1, 1.5, 2)
target_cta_ratios = [0, 0.5, 1.0, 1.5, 2.0];
plotted_indices = false(1, length(target_cta_ratios)); % To avoid duplicate plots if dt is small

% Always plot initial condition
plot(x, u(:,1), 'DisplayName', sprintf('c*t/a = %.2f (t=%.3fs)', c_wave*t_vec(1)/a, t_vec(1)));

colors = lines(length(target_cta_ratios)); % Get distinct colors

for i_plot_time = 1:length(target_cta_ratios)
    target_t = target_cta_ratios(i_plot_time) * a / c_wave;
    [~, best_idx] = min(abs(t_vec - target_t)); % Find closest time step
    
    % Check if this time or a very close one has already been plotted (especially for t=0)
    already_plotted = false;
    if best_idx == 1 && i_plot_time > 1 % t=0 already plotted
        already_plotted = true;
    end

    if ~already_plotted
        plot(x, u(:,best_idx), 'Color', colors(i_plot_time,:), 'DisplayName', sprintf('c*t/a = %.2f (t=%.3fs)', c_wave*t_vec(best_idx)/a, t_vec(best_idx)));
    end
end

hold off;
xlabel('Position x (m)');
ylabel('Displacement u (m)');
title(sprintf('String Shape (dx=%.3f, dt=%.4f, C=%.2f)', dx_val, dt_val, C_stab));
legend show;
grid on;
ylim([-1.1*y0_peak, 1.1*y0_peak]); % Consistent y-axis scaling

% Animation (optional)
% prompt_anim = input('Show animation? (y/n): ', 's');
% if lower(prompt_anim) == 'y'
%     figure;
%     h_anim = plot(x, u(:,1));
%     ax_anim = gca;
%     ax_anim.XLim = [0 a];
%     ax_anim.YLim = [-1.1*y0_peak, 1.1*y0_peak];
%     title_anim = title(sprintf('String Vibration - c*t/a = %.2f', c_wave*t_vec(1)/a));
%     xlabel('Position x (m)');
%     ylabel('Displacement u (m)');
%     grid on;
% 
%     for n_anim = 1:nt
%         set(h_anim, 'YData', u(:,n_anim));
%         cta_ratio = c_wave * t_vec(n_anim) / a;
%         set(title_anim, 'String', sprintf('String Vibration - c*t/a = %.2f (t=%.3fs)', cta_ratio, t_vec(n_anim)));
%         drawnow;
%         % pause(0.01); % Adjust for animation speed
%     end
% end
% 
% fprintf('Simulation finished.\n');