function replicate_fig5()
% REPLICATE_FIG5 Reproduces Figure 5 from Hastings & Powell (1991)
% Translating Julia script logic to MATLAB using Event Detection for Poincaré maps.

    % --- General Settings ---
    tspan = [0, 500000];
    u0 = [0.76; 0.16; 9.9]; % Initial conditions [x, y, z]
    % Use high precision to capture chaotic dynamics accurately
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Create Figure Layout
    figure('Color', 'w', 'Name', 'Hastings & Powell Fig 5', 'Position', [100, 100, 1000, 800]);

    %% ==========================================
    %  PART 1: b1 = 3.0 (Figures 5A and 5B)
    %  ==========================================
    % Parameters: [a1, a2, b1, b2, d1, d2]
    p_b3 = [5.0, 0.1, 3.0, 2.0, 0.4, 0.01]; 
    z_plane_b3 = 9.0;
    
    % Configure Event to detect crossing z=9 from above
    opts_b3 = odeset(options, 'Events', @(t,y) event_z_plane(t,y, z_plane_b3));
    
    % Run Simulation
    fprintf('Simulating b1 = 3.0 ...\n');
    [~, ~, ~, YE, ~] = ode45(@(t,y) food_chain_ode(t,y,p_b3), tspan, u0, opts_b3);
    
    % Filter points to isolate the "handle" of the attractor
    % Julia conditions: x in [0.9, 1.0], y in [0.0, 0.1]
    mask_b3 = (YE(:,1) >= 0.9) & (YE(:,1) <= 1.0) & ...
              (YE(:,2) >= 0.0) & (YE(:,2) <= 0.1);
          
    xn_b3 = YE(mask_b3, 1);
    yn_b3 = YE(mask_b3, 2);
    
    % --- Plot 5A: Poincaré Section ---
    subplot(2, 2, 1);
    plot(xn_b3, yn_b3, 'k.', 'MarkerSize', 4);
    xlim([0.95, 0.983]); ylim([0.015, 0.04]);
    xlabel('x(n)'); ylabel('y(n)');
    title('A: Poincaré Section (b_1 = 3.0)');
    grid off; box on;

    % --- Plot 5B: Poincaré Map (x_n vs x_{n+1}) ---
    subplot(2, 2, 2);
    plot(xn_b3(1:end-1), xn_b3(2:end), 'k.', 'MarkerSize', 4);
    hold on;
    plot([0.95 0.98], [0.95 0.98], 'k-'); % 1:1 diagonal line
    hold off;
    xlim([0.95, 0.98]); ylim([0.95, 0.98]);
    xlabel('x(n)'); ylabel('x(n+1)');
    title('B: Poincaré Map (b_1 = 3.0)');
    grid off; box on;

    %% ==========================================
    %  PART 2: b1 = 6.0 (Figures 5C and 5D)
    %  ==========================================
    p_b6 = [5.0, 0.1, 6.0, 2.0, 0.4, 0.01];
    z_plane_b6 = 3.0;
    
    % Configure Event to detect crossing z=3 from above
    opts_b6 = odeset(options, 'Events', @(t,y) event_z_plane(t,y, z_plane_b6));
    
    % Run Simulation
    fprintf('Simulating b1 = 6.0 ...\n');
    [~, ~, ~, YE6, ~] = ode45(@(t,y) food_chain_ode(t,y,p_b6), tspan, u0, opts_b6);
    
    % Filter points
    % Julia conditions: x in [0.93, 1.0], y in [0.0, 0.085]
    mask_b6 = (YE6(:,1) >= 0.93) & (YE6(:,1) <= 1.003) & ...
              (YE6(:,2) >= 0.0) & (YE6(:,2) <= 0.085);
          
    xn_b6 = YE6(mask_b6, 1);
    yn_b6 = YE6(mask_b6, 2);

    % --- Plot 5C: Poincaré Section ---
    subplot(2, 2, 3);
    plot(xn_b6, yn_b6, 'k.', 'MarkerSize', 4);
    xlim([0.93, 1.003]); ylim([-0.003, 0.09]);
    xlabel('x(n)'); ylabel('y(n)');
    title('C: Poincaré Section (b_1 = 6.0)');
    grid off; box on;

    % --- Plot 5D: Poincaré Map ---
    subplot(2, 2, 4);
    plot(xn_b6(1:end-1), xn_b6(2:end), 'k.', 'MarkerSize', 4);
    hold on;
    plot([0.93 1.003], [0.93 1.003], 'k-'); % 1:1 diagonal line
    hold off;
    xlim([0.93, 1.003]); ylim([0.93, 1.003]);
    xlabel('x(n)'); ylabel('x(n+1)');
    title('D: Poincaré Map (b_1 = 6.0)');
    grid off; box on;
    
    fprintf('Done.\n');
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
    isterminal = 0;   % Do not stop the integration
    direction = -1;   % Only detect crossings from positive to negative (from above)
end