function visualize_food_chain(b1_override,b2_override, u0_override, z_plane_override)
% VISUALIZE_FOOD_CHAIN Interactive 3D trajectory plotter (No Gradient).
%
% Usage:
%   visualize_food_chain()              -> Standard Chaos (b1=3.0)
%   visualize_food_chain(2.5)           -> Limit Cycle (b1=2.5)
%   visualize_food_chain(3.0, [0.5;0.5;8]) -> Custom Start

    %% --- 1. Configuration & Defaults ---
    p_default = [5.0, 0.1, 3.0, 2.0, 0.4, 0.01]; % [a1, a2, b1, b2, d1, d2]
    u0_default = [0.76; 0.16; 9.9];
    z_default = 9.0;

    % Handle User Inputs
    if nargin < 1 || isempty(b1_override)
        p = p_default;
        fprintf('Using default b1 = 3.0 (Chaos)\n');
    else
        p = p_default;
        p(3) = b1_override; 
        fprintf('Using custom b1 = %.2f\n', b1_override);
    end

    if nargin < 2 || isempty(b2_override)
        fprintf('Using default b2 = 2.0 (Chaos)\n');
    else 
        p(4) = b2_override; 
        fprintf('Using custom b2 = %.2f\n', b2_override);
    end

    if nargin < 3 || isempty(u0_override)
        u0 = u0_default;
    else
        u0 = u0_override(:);
    end

    if nargin < 4 || isempty(z_plane_override)
        % Auto-adjust plane if b1 is very high (e.g. b1=6 case)
        if p(3) > 5.0
            z_plane = 3.0;
        else
            z_plane = z_default;
        end
    else
        z_plane = z_plane_override;
    end

    %% --- 2. Simulation ---
    tspan = [0, 10000]; 
    
    % Solver Options (High precision + Event Detection)
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, ...
                     'Events', @(t,y) event_z_plane(t,y, z_plane));

    fprintf('Simulating... ');
    [T, Y, TE, YE, ~] = ode45(@(t,y) food_chain_ode(t,y,p), tspan, u0, options);
    fprintf('Done.\n');

    %% --- 3. Visualization ---
    figure('Color', 'w', 'Name', sprintf('Trajectory b1=%.2f', p(3)));
    
    % Remove initial transients (first 20%)
    start_idx = round(length(T) * 0.2); 
    
    % A. Plot Continuous Trajectory (Black Line)
    plot3(Y(start_idx:end,1), Y(start_idx:end,2), Y(start_idx:end,3), ...
          'Color', [0,0,0, 0.6], 'LineWidth', 0.5); 
    hold on;

    % B. Plot Poincaré Plane (Transparent Blue Surface)
    x_min = min(Y(:,1)); x_max = max(Y(:,1));
    y_min = min(Y(:,2)); y_max = max(Y(:,2));
    
    [Xp, Yp] = meshgrid([x_min x_max], [y_min y_max]);
    surf(Xp, Yp, ones(size(Xp))*z_plane, ...
         'FaceColor', 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    % C. Plot Intersection Points (Red Dots - No Gradient)
    valid_pts = TE > T(start_idx);
    
    if any(valid_pts)
        pts_x = YE(valid_pts, 1);
        pts_y = YE(valid_pts, 2);
        pts_z = YE(valid_pts, 3);
        
        % Simple red markers
        plot3(pts_x, pts_y, pts_z, 'r.', 'MarkerSize', 12);
        
        % Optional: Add text label for the plane
        text(x_max, y_max, z_plane + 0.1, sprintf(' z=%.1f', z_plane), ...
             'Color', 'b', 'FontSize', 8);
    else
        warning('No intersections found with plane z=%.2f', z_plane);
    end

    % D. Styling
    grid on; box on;
    view(-35, 20); 
    xlabel('Prey (x)'); ylabel('Consumer (y)'); zlabel('Predator (z)');
    title(sprintf('Food Chain Dynamics (b_1 = %.2f)', p(3)));
    
    zlim([min(Y(:,3))*0.9, max(Y(:,3))*1.1]);
    
    hold off;
end

%% --- Helper Functions ---
function dydt = food_chain_ode(~, u, p)
    x = u(1); y = u(2); z = u(3);
    a1=p(1); a2=p(2); b1=p(3); b2=p(4); d1=p(5); d2=p(6);
    
    f1 = (a1 * x) / (1 + b1 * x);
    f2 = (a2 * y) / (1 + b2 * y);
    
    dydt = [x*(1 - x) - f1*y; 
            f1*y - f2*z - d1*y; 
            f2*z - d2*z];
end

function [value, isterminal, direction] = event_z_plane(~, u, z_val)
    value = u(3) - z_val;
    isterminal = 0;   
    direction = -1; % Detect crossings from above
end