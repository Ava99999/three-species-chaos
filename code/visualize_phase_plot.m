function visualize_phase_plot(b1_override, b2_override, u0_override)
% VISUALIZE_PHASE_PLOT Robust 3D plotter with diagnostics.
%
% Usage:
%   visualize_phase_plot()

    %% --- 1. Configuration ---
    p_default = [5.0, 0.1, 3.0, 2.0, 0.4, 0.01]; 
    u0_default = [0.76; 0.16; 9.9];
    
    % Handle Inputs
    if nargin < 1 || isempty(b1_override), p = p_default;
    else, p = p_default; p(3) = b1_override; end

    if nargin < 2 || isempty(b2_override), p(4) = 2.0;
    else, p(4) = b2_override; end

    if nargin < 3 || isempty(u0_override), u0 = u0_default;
    else, u0 = u0_override(:); end

    %% --- 2. Simulation ---
    % Reduced tspan slightly to ensure responsiveness. 
    % 100,000 can sometimes freeze slower computers.
    tspan = [0, 20000]; 
    
    % Options: Enforce NonNegative to prevent singularities (x < 0)
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'NonNegative', [1 2 3]);

    fprintf('Simulating (t=0 to %d)...\n', tspan(2));
    [T, Y] = ode45(@(t,y) food_chain_ode(t,y,p), tspan, u0, options);
    
    % --- DIAGNOSTICS: Check Data Validity ---
    if isempty(Y)
        error('Error: Simulation returned no data.');
    end
    
    fprintf('Data Points Generated: %d\n', length(T));
    fprintf('X Range: [%.4f, %.4f]\n', min(Y(:,1)), max(Y(:,1)));
    fprintf('Y Range: [%.4f, %.4f]\n', min(Y(:,2)), max(Y(:,2)));
    fprintf('Z Range: [%.4f, %.4f]\n', min(Y(:,3)), max(Y(:,3)));

    if any(isnan(Y(:))) || any(isinf(Y(:)))
        error('Error: Simulation contains NaN or Inf values. System unstable.');
    end

    %% --- 3. Visualization ---
    fprintf('Plotting...\n');
    
    % Create Figure with specific renderer to avoid "Grey Window" freeze
    f = figure('Color', 'w', 'Name', sprintf('Phase Plot b1=%.2f', p(3)));
    set(f, 'Renderer', 'painters'); % safer for simple line plots
    clf; % Clear any existing garbage
    
    % Remove transients
    start_idx = round(length(T) * 0.2); 
    if start_idx == 0, start_idx = 1; end
    
    % Plot
    plot3(Y(start_idx:end,1), Y(start_idx:end,2), Y(start_idx:end,3), ...
          'Color', [0.2, 0.2, 0.2], 'LineWidth', 0.5); 
    
    % Force draw
    grid on; box on;
    xlabel('Prey (x)'); ylabel('Consumer (y)'); zlabel('Predator (z)');
    title(sprintf('Phase Plot (b_1=%.2f)', p(3)));
    view(-35, 20); 
    axis tight;
    
    drawnow; % Force MATLAB to render immediately
    fprintf('Done.\n');
end

%% --- Helper Functions ---
function dydt = food_chain_ode(~, u, p)
    x = u(1); y = u(2); z = u(3);
    f1 = (p(1)*x)/(1+p(3)*x); 
    f2 = (p(2)*y)/(1+p(4)*y);
    dydt = [x*(1-x)-f1*y; f1*y-f2*z-p(5)*y; f2*z-p(6)*z];
end