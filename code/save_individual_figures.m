function save_individual_figures()
% SAVE_INDIVIDUAL_FIGURES Generates and saves Figs 5A-5E with TIME GRADIENTS.
% Based on Hastings & Powell (1991) parameters.

    % --- General Settings ---
    % Increased time slightly to ensure the gradient is clearly visible
    tspan = [0, 200000]; 
    u0 = [0.76; 0.16; 9.9];
    % High precision settings for chaotic systems
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 

    %% ==========================================
    %  GROUP 1: b1 = 3.0 (Figs 5A, 5B, 5E)
    %  ==========================================
    fprintf('Simulating b1 = 3.0 ...\n');
    p_b3 = [5.0, 0.1, 3.0, 2.0, 0.4, 0.01]; 
    z_plane_b3 = 9.0;
    
    % Event detection for Poincare section
    opts_b3 = odeset(options, 'Events', @(t,y) event_z_plane(t,y, z_plane_b3));
    
    % Run Simulation
    [T, Y, TE, YE, ~] = ode45(@(t,y) food_chain_ode(t,y,p_b3), tspan, u0, opts_b3);
    
    % Filter points (Handle of the teacup)
    % Ranges: x in [0.9, 1.0], y in [0.0, 0.1]
    mask_b3 = (YE(:,1) >= 0.9) & (YE(:,1) <= 1.0) & ...
              (YE(:,2) >= 0.0) & (YE(:,2) <= 0.1);
          
    xn_b3 = YE(mask_b3, 1);
    yn_b3 = YE(mask_b3, 2);

    % --- SAVE FIGURE 5A (Poincare Section with Gradient) ---
    f1 = figure('Visible', 'off'); 
    
    % Create color vector based on order of appearance
    colors_A = 1:length(xn_b3);
    
    % Scatter plot with gradient
    scatter(xn_b3, yn_b3, 10, colors_A, 'filled');
    colormap(turbo); 
    cb = colorbar; cb.Label.String = 'Iteration (Time)';
    
    xlim([0.95, 0.983]); ylim([0.015, 0.04]);
    xlabel('x(n)'); ylabel('y(n)');
    title('Figure 5A: Poincaré Section (b_1=3.0)');
    grid off; box on;
    
    fprintf('Saving Fig5A_gradient.png...\n');
    exportgraphics(f1, 'Fig5A_gradient.png', 'Resolution', 300);
    close(f1);

    % --- SAVE FIGURE 5B (Poincare Map with Gradient) ---
    f2 = figure('Visible', 'off');
    
    % Prepare data for map (x_n vs x_n+1)
    x_current = xn_b3(1:end-1);
    x_next = xn_b3(2:end);
    colors_B = 1:length(x_current); % Gradient matches the number of pairs
    
    % Scatter plot
    scatter(x_current, x_next, 10, colors_B, 'filled');
    colormap(turbo);
    cb = colorbar; cb.Label.String = 'Iteration (Time)';
    
    hold on; 
    plot([0.95 0.98], [0.95 0.98], 'k-', 'LineWidth', 0.5); % Diagonal
    hold off;
    
    xlim([0.95, 0.98]); ylim([0.95, 0.98]);
    xlabel('x(n)'); ylabel('x(n+1)');
    title('Figure 5B: Poincaré Map (b_1=3.0)');
    grid off; box on;
    
    fprintf('Saving Fig5B_gradient.png...\n');
    exportgraphics(f2, 'Fig5B_gradient.png', 'Resolution', 300);
    close(f2);

    % --- SAVE FIGURE 5E (3D Attractor with Colored Intersections) ---
    f3 = figure('Visible', 'off');
    
    % 1. Plot continuous trajectory (Black line)
    start_idx = round(length(T) * 0.2); % Skip first 20% transients
    plot3(Y(start_idx:end,1), Y(start_idx:end,2), Y(start_idx:end,3), 'k-', 'LineWidth', 0.1);
    hold on;
    
    % 2. Draw the transparent plane z=9
    [Xp, Yp] = meshgrid([0.5 1.0], [0 0.5]);
    surf(Xp, Yp, ones(size(Xp))*z_plane_b3, 'FaceColor', 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    % 3. Plot intersection points with Gradient
    % We use scatter3 here to apply the same color logic to the 3D dots
    z_points = ones(size(xn_b3)) * z_plane_b3;
    scatter3(xn_b3, yn_b3, z_points, 15, colors_A, 'filled');
    colormap(turbo);
    
    view(-35, 20); 
    xlabel('x'); ylabel('y'); zlabel('z');
    zlim([7.5 10.5]);
    title('Figure 5E: 3D Attractor');
    grid on;
    
    fprintf('Saving Fig5E_gradient.png...\n');
    exportgraphics(f3, 'Fig5E_gradient.png', 'Resolution', 300);
    close(f3);

    %% ==========================================
    %  GROUP 2: b1 = 6.0 (Figs 5C, 5D)
    %  ==========================================
    fprintf('Simulating b1 = 6.0 ...\n');
    p_b6 = [5.0, 0.1, 6.0, 2.0, 0.4, 0.01];
    z_plane_b6 = 3.0;
    
    opts_b6 = odeset(options, 'Events', @(t,y) event_z_plane(t,y, z_plane_b6));
    [~, ~, ~, YE6, ~] = ode45(@(t,y) food_chain_ode(t,y,p_b6), tspan, u0, opts_b6);
    
    % Filter points
    % Ranges: x in [0.93, 1.003], y in [0.0, 0.085]
    mask_b6 = (YE6(:,1) >= 0.93) & (YE6(:,1) <= 1.003) & ...
              (YE6(:,2) >= 0.0) & (YE6(:,2) <= 0.085);
          
    xn_b6 = YE6(mask_b6, 1);
    yn_b6 = YE6(mask_b6, 2);

    % --- SAVE FIGURE 5C (Poincare Section with Gradient) ---
    f4 = figure('Visible', 'off');
    
    colors_C = 1:length(xn_b6);
    scatter(xn_b6, yn_b6, 10, colors_C, 'filled');
    colormap(turbo);
    cb = colorbar; cb.Label.String = 'Iteration (Time)';

    xlim([0.93, 1.003]); ylim([-0.003, 0.09]);
    xlabel('x(n)'); ylabel('y(n)');
    title('Figure 5C: Poincaré Section (b_1=6.0)');
    grid off; box on;
    
    fprintf('Saving Fig5C_gradient.png...\n');
    exportgraphics(f4, 'Fig5C_gradient.png', 'Resolution', 300);
    close(f4);

    % --- SAVE FIGURE 5D (Poincare Map with Gradient) ---
    f5 = figure('Visible', 'off');
    
    x_curr_6 = xn_b6(1:end-1);
    x_next_6 = xn_b6(2:end);
    colors_D = 1:length(x_curr_6);
    
    scatter(x_curr_6, x_next_6, 10, colors_D, 'filled');
    colormap(turbo);
    cb = colorbar; cb.Label.String = 'Iteration (Time)';
    
    hold on; 
    plot([0.93 1.003], [0.93 1.003], 'k-', 'LineWidth', 0.5); 
    hold off;
    
    xlim([0.93, 1.003]); ylim([0.93, 1.003]);
    xlabel('x(n)'); ylabel('x(n+1)');
    title('Figure 5D: Poincaré Map (b_1=6.0)');
    grid off; box on;
    
    fprintf('Saving Fig5D_gradient.png...\n');
    exportgraphics(f5, 'Fig5D_gradient.png', 'Resolution', 300);
    close(f5);
    
    fprintf('All figures saved successfully.\n');
end

%% --- Helper Functions ---

function dydt = food_chain_ode(~, u, p)
    x = u(1); y = u(2); z = u(3);
    a1=p(1); a2=p(2); b1=p(3); b2=p(4); d1=p(5); d2=p(6);
    
    % Functional Responses
    f1 = (a1 * x) / (1 + b1 * x);
    f2 = (a2 * y) / (1 + b2 * y);
    
    % ODE System
    dydt = [x*(1 - x) - f1*y; 
            f1*y - f2*z - d1*y; 
            f2*z - d2*z];
end

function [value, isterminal, direction] = event_z_plane(~, u, z_val)
    value = u(3) - z_val;
    isterminal = 0;   
    direction = -1; % Detect crossings from above
end