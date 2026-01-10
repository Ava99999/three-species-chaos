function visualize_poincare_map(b1_in, b2_in, z_plane_in, u0_in)

    %% --- 1. Configuration & Defaults ---
    
    % Standard fixed parameters: [a1, a2, b1, b2, d1, d2]
    p_base = [5, 0.1, 6.2, 2.0, 0.4, 0.01];
    
    % Handle Input Arguments
    if nargin < 1 || isempty(b1_in), b1 = 3.0; else, b1 = b1_in; end
    if nargin < 2 || isempty(b2_in), b2 = 2.0; else, b2 = b2_in; end
    if nargin < 3 || isempty(z_plane_in)
        % Auto-guess plane if not provided (High b1 usually needs lower z)
        if b1 > 5.0, z_plane = 3.0; else, z_plane = 9.0; end
    else
        z_plane = z_plane_in;
    end
    if nargin < 4 || isempty(u0_in), u0 = [0.76; 0.16; 9.9]; else, u0 = u0_in; end

    % Construct parameter vector
    p = p_base;
    p(3) = b1; 
    p(4) = b2;

    fprintf('Settings: b1=%.2f, b2=%.2f, Plane z=%.2f\n', b1, b2, z_plane);

    %% --- 2. Simulation ---
    tspan = [0, 80000]; % Long duration to gather enough points
    
    % Event Detection: Cross z_plane from above
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, ...
                     'Events', @(t,y) event_z_plane(t,y, z_plane));

    fprintf('Simulating trajectory... ');
    [T, ~, ~, YE, ~] = ode45(@(t,y) food_chain_ode(t,y,p), tspan, u0, options);
    fprintf('Done.\n');

    %% --- 3. Data Processing ---
    
    if isempty(YE)
        warning('No intersection points found! The trajectory never crossed z=%.2f.', z_plane);
        return;
    end

    % Remove transients: Discard the first 20% of intersection points
    % This ensures we are plotting the attractor, not the approach.
    cutoff_idx = round(size(YE, 1) * 0.2);
    if cutoff_idx == 0, cutoff_idx = 1; end
    
    valid_YE = YE(cutoff_idx:end, :);
    
    if isempty(valid_YE)
        warning('All intersections occurred during the transient phase.');
        return;
    end

    xn = valid_YE(:, 1);
    yn = valid_YE(:, 2);

    %% --- 4. Visualization ---
    figure('Color', 'w', 'Name', sprintf('Poincaré Analysis (b1=%.1f, b2=%.1f)', b1, b2));

    % --- Left Plot: Poincaré Section (y vs x on the plane) ---
    subplot(1, 2, 1);
    plot(xn, yn, 'k.', 'MarkerSize', 4);
    
    xlabel('Prey x(n)'); 
    ylabel('Consumer y(n)');
    title(sprintf('Poincaré Section\n(z = %.1f)', z_plane));
    grid on; box on; axis tight;
    
    % Add a tiny margin to axes so points aren't on the edge
    xl = xlim; yl = ylim;
    xlim([xl(1)-0.005, xl(2)+0.005]);
    ylim([yl(1)-0.005, yl(2)+0.005]);

    % --- Right Plot: Poincaré Map (x_n+1 vs x_n) ---
    subplot(1, 2, 2);
    
    % We plot x(n) vs x(n+1)
    x_current = xn(1:end-1);
    x_next    = xn(2:end);
    
    plot(x_current, x_next, 'k.', 'MarkerSize', 4);
    hold on;
    
    % Add 1:1 Diagonal Line for reference
    min_val = min(xn); max_val = max(xn);
    plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 1);
    hold off;
    
    xlabel('x(n)'); 
    ylabel('x(n+1)');
    title('Poincaré Map (1D Return)');
    grid on; box on; axis tight;
    
    % Add margins
    xl = xlim; yl = ylim;
    xlim([xl(1)-0.002, xl(2)+0.002]);
    ylim([yl(1)-0.002, yl(2)+0.002]);

end

%% --- Helper Functions ---

function dydt = food_chain_ode(~, u, p)
    % Unpack state and parameters
    x = u(1); y = u(2); z = u(3);
    a1=p(1); a2=p(2); b1=p(3); b2=p(4); d1=p(5); d2=p(6);
    
    % Functional responses
    f1 = (a1 * x) / (1 + b1 * x);
    f2 = (a2 * y) / (1 + b2 * y);
    
    % Differential equations
    dx = x*(1 - x) - f1*y;
    dy = f1*y - f2*z - d1*y;
    dz = f2*z - d2*z;
    
    dydt = [dx; dy; dz];
end

function [value, isterminal, direction] = event_z_plane(~, u, z_val)
    % Event function to detect when trajectory crosses z = z_val
    value = u(3) - z_val;
    isterminal = 0;   % Continue simulation
    direction = -1;   % Only detect crossings from above
end