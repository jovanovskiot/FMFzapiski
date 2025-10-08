% solve_heat_equation.m
% Solves the 1D heat equation using FTCS and Crank-Nicolson methods.

clear;
clc;
close all;

% --- Parameters ---
a = 1.0;         % Length of the slab (m)
T0 = 50.0;       % Initial temperature in the heated region (deg C)
D_const = 0.1;   % Thermal diffusivity (m^2/s) (K in Python scripts)

% Spatial grid
dx_val = 0.02;   % Spatial step (m) (h in problem doc)
nx = round(a / dx_val) + 1;
x = linspace(0, a, nx);

% Temporal grid
t_max = 2.0;     % Maximum simulation time (s)
% For FTCS, dt must satisfy D*dt/dx^2 < 0.5
% For Crank-Nicolson, it's unconditionally stable, but accuracy depends on dt.

% Choose method: 'FTCS' or 'CN' (Crank-Nicolson)
method = 'CN'; % Change to 'FTCS' to run the other method

if strcmp(method, 'FTCS')
    % Ensure stability for FTCS: r_ftcs = D*dt/dx^2 < 0.5
    % Let r_ftcs = 0.45
    r_ftcs_stab = 0.45;
    dt_val = r_ftcs_stab * dx_val^2 / D_const;
    fprintf('Using FTCS method.\n');
elseif strcmp(method, 'CN')
    % For Crank-Nicolson, dt can be larger. Let's pick a value.
    % A common choice is dt ~ dx for C-N if D is around 1.
    % Or dt relative to dx^2/D for comparison with FTCS critical dt.
    dt_val = 0.5 * dx_val^2 / D_const; % Example: similar to FTCS critical step
    % dt_val = 0.01; % Or a fixed dt for C-N
    fprintf('Using Crank-Nicolson method.\n');
else
    error('Invalid method specified. Choose "FTCS" or "CN".');
end

nt = round(t_max / dt_val) + 1;
t_vec = linspace(0, t_max, nt);
dt_val = t_vec(2) - t_vec(1); % Recalculate dt based on nt for precision

fprintf('Spatial step dx = %.4f m, Number of spatial points nx = %d\n', dx_val, nx);
fprintf('Temporal step dt = %.4f s, Number of temporal steps nt = %d\n', dt_val, nt);
fprintf('Total simulation time t_max = %.2f s\n', t_max);

% --- Initial Conditions ---
u = zeros(nx, nt); % Temperature matrix u(x,t)
% Ambient temperature is 0, heated region [0.2a, 0.4a] is T0
idx_heated = (x >= 0.2*a) & (x <= 0.4*a);
u(idx_heated, 1) = T0;

% --- Boundary Conditions ---
% Kept at ambient temperature (0)
u(1, :) = 0;
u(nx, :) = 0;

% --- Source Term ---
% q = 0, so Q_term = 0 as per problem task 1 description
% kQ in problem doc is dt_val * Q_source_term. Here Q_source_term is 0.

% --- Simulation ---
if strcmp(method, 'FTCS')
    r_ftcs = D_const * dt_val / dx_val^2;
    fprintf('FTCS stability factor r = D*dt/dx^2 = %.4f\n', r_ftcs);
    if r_ftcs >= 0.5
        warning('FTCS stability condition r < 0.5 might be violated. r = %.4f', r_ftcs);
    end

    for n = 1:nt-1 % Time loop (from t_n to t_{n+1})
        for m = 2:nx-1 % Spatial loop (internal points)
            u(m, n+1) = u(m,n) + r_ftcs * (u(m+1,n) - 2*u(m,n) + u(m-1,n));
        end
        % Boundary conditions are already set and maintained
    end

elseif strcmp(method, 'CN')
    r_cn = D_const * dt_val / dx_val^2; % This is r from problem doc (page 1 bottom, and page 2 top)
    alpha = r_cn / 2.0; % This is D*dt / (2*dx^2)

    % LHS Matrix (for internal points u_2, ..., u_{nx-1})
    % Size of internal_points_matrix is (nx-2) x (nx-2)
    num_internal_pts = nx - 2;
    
    if num_internal_pts <= 0 % Handle cases with very few points (e.g., nx=1 or nx=2)
        % No internal points to solve for, solution is just boundary conditions
        % This loop will effectively do nothing if num_internal_pts is 0 or negative
    else
        main_diag_LHS = (1 + 2*alpha) * ones(num_internal_pts, 1);
        
        if num_internal_pts > 1
            sub_diag_LHS  = -alpha * ones(num_internal_pts - 1, 1);
            sup_diag_LHS  = -alpha * ones(num_internal_pts - 1, 1);
        else % Only one internal point, no sub/super diagonals for the LHS_M matrix
            sub_diag_LHS = []; % Empty
            sup_diag_LHS = []; % Empty
        end
        
        % Construct the B matrix for spdiags
        % B_mat_for_spdiags must have num_internal_pts rows.
        % Column 1: sub-diagonal (d=-1), pad with 0 at the top.
        % Column 2: main-diagonal (d=0).
        % Column 3: super-diagonal (d=1), pad with 0 at the bottom.
        
        % Initialize B_mat_for_spdiags with main diagonal
        B_mat_for_spdiags = main_diag_LHS; 
        diag_indices = 0;

        if num_internal_pts > 1
            B_mat_for_spdiags = [[0; sub_diag_LHS], main_diag_LHS, [sup_diag_LHS; 0]];
            diag_indices = [-1, 0, 1];
        end
            
        LHS_M = spdiags(B_mat_for_spdiags, diag_indices, num_internal_pts, num_internal_pts);

        for n = 1:nt-1 % Time loop
            % RHS vector
            RHS_vec = zeros(num_internal_pts, 1);
            
            % Contribution from u^n terms
            % This loop needs to handle the cases for num_internal_pts = 1 carefully
            if num_internal_pts == 1
                m = 2; % The only internal point index
                RHS_vec(1) = alpha * u(m-1,n) + (1 - 2*alpha) * u(m,n) + alpha * u(m+1,n);
            else
                for m_int = 1:num_internal_pts % m_int is index in internal points vector (1 to nx-2)
                    m = m_int + 1; % m is index in full u vector (2 to nx-1)
                    RHS_vec(m_int) = alpha * u(m-1,n) + (1 - 2*alpha) * u(m,n) + alpha * u(m+1,n);
                end
            end
            
            % Add contributions from known boundary conditions u(1,n+1) and u(nx,n+1)
            % These are zero in this problem setup.
            % If u(1,n+1) (g0(t_{n+1})) was non-zero and num_internal_pts > 0:
            % RHS_vec(1) = RHS_vec(1) + alpha * u(1, n+1); 
            % If u(nx,n+1) (g1(t_{n+1})) was non-zero and num_internal_pts > 0:
            % RHS_vec(num_internal_pts) = RHS_vec(num_internal_pts) + alpha * u(nx, n+1);
            
            % Solve for u at internal points at time n+1
            if num_internal_pts > 0
                u_internal_new = LHS_M \ RHS_vec;
                u(2:nx-1, n+1) = u_internal_new;
            end
            
            % Boundary conditions u(1,n+1) and u(nx,n+1) are already set.
        end
    end
end

% --- Plotting ---
figure;
hold on;
plot_indices = round(linspace(1, nt, min(nt, 6))); % Plot up to 6 profiles
colors = lines(length(plot_indices)); % Get distinct colors

for i = 1:length(plot_indices)
    idx = plot_indices(i);
    plot(x, u(:,idx), 'Color', colors(i,:), 'DisplayName', sprintf('t = %.2f s', t_vec(idx)));
end
hold off;
xlabel('Position x (m)');
ylabel('Temperature T (deg C)');
title(sprintf('Temperature Profile (%s, dx=%.3f, dt=%.4f)', method, dx_val, dt_val));
legend show;
grid on;

%% Animation (optional, can be slow for large grids)
% prompt_anim = input('Show animation? (y/n): ', 's');
% if lower(prompt_anim) == 'y'
%     figure;
%     h_anim = plot(x, u(:,1));
%     ax_anim = gca;
%     ax_anim.XLim = [0 a];
%     ax_anim.YLim = [min(u(:))-0.1*T0, max(u(:))+0.1*T0]; % Adjust YLims based on data
%     if ax_anim.YLim(1) == ax_anim.YLim(2) % Handle case of all zero temperature
%         ax_anim.YLim = [-1 T0+1];
%     end
%     title_anim = title(sprintf('Temperature Evolution (%s) - t = %.2f s', method, t_vec(1)));
%     xlabel('Position x (m)');
%     ylabel('Temperature T (deg C)');
%     grid on;
% 
%     for n_anim = 1:nt
%         set(h_anim, 'YData', u(:,n_anim));
%         set(title_anim, 'String', sprintf('Temperature Evolution (%s) - t = %.2f s', method, t_vec(n_anim)));
%         drawnow;
%         % pause(0.01); % Adjust for animation speed
%     end
% end
% 
% fprintf('Simulation finished.\n');