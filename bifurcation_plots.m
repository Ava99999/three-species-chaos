clear all
close all

chaos_parameters

%% defining three species system
f1 = @(u) param.a1*u/(1 + param.b1*u); % functional response f1
f2 = @(u) param.a2*u/(1 + param.b2*u); % functional response f2

%% parameters ode
tspan = param.t; %[t0,T];
ic = param.ic1; % initial condition from reprod. paper
options = odeset('RelTol',1e-8,'abstol',1e-8*[1,1,1]);

%% bifurcation
Z = zeros(param.NB, length(param.b1_r)); % record timestep values
Z_locmax = zeros(param.NB, length(param.b1_r)); % record local maxima
Z_max = zeros(length(param.b1_r),1); % record global maxima
Z_min = zeros(length(param.b1_r),1); % record global minima
errors = zeros(length(param.b1_r),1); % in case no cycling behavior is found
margin = 0.66;
% TODO i don't think that will occur here because i forced the length to be
% 10000

tic
parfor j = 1:length(param.b1_r)
    % solve the ode 
    f1_loc = @(u) param.a1*u/(1 + param.b1_r(j)*u); % functional response f1 local
    [t,X] = ode45(@(t,y) superpredators(t,y,f1_loc,f2,param.d1,param.d2),tspan, ic, options);
    
    % sanity check if using dynamic number of steps
    if length(X(:,3)) > param.NB
        Z(:,j) = X(end-param.NB+1:end,3); % select last NB steps
        % global max and min
        Z_max(j) = max(Z(:,j));
        Z_min(j) = min(Z(:,j));
    else
        errors(j) = j; % save index of unstable solution
    end
end
toc

tic
% select the local maxima that are large enough compared to the cycle
for j = 1:length(param.b1_r)
   for i = 2:(param.NB-1) % todo does the first index matter?
       % identify maxima
       if Z(i,j) >= Z(i-1,j) && Z(i,j) >= Z(i+1,j)
            Z_range = Z_max(j) - Z_min(j);
            % only select those values that are large enough
            if Z(i,j) > (Z_min(j) + margin*Z_range)
                Z_locmax(i,j) = Z(i,j);
                
                if param.b1_r(j)==3.0
                    fprintf("value of b1 ")
                    param.b1_r(j)
                    fprintf("maximum") 
                    Z(i,j)
                    fprintf("value before")
                    Z(i-1,j)
                    fprintf("value after")
                    Z(i+1,j)
                end
            end
       end
   end
end
toc

%% plot results
% copies of the parameter range, long row vector
B1 = repelem(param.b1_r,1,param.NB); 
Z_temp = reshape(Z_locmax,1,[]);
plot_matrix = [B1; Z_temp];
% remove duplicate maxima/minima
plot_points = unique(plot_matrix',"rows", 'stable'); 
indices = find(plot_points(:,2)==0);
plot_points(indices,:) = [];

% plot of all points
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 4;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.25 2.6]);
ylim([11.4 12.8]);
xlabel("b_1");
ylabel("z_{max}")
title('c');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gcf,'Renderer','zbuffer')
sz = 2;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.2 6.5]);
ylim([3 13]);
xlabel("b_1");
ylabel("z_{max}")

% figure 4a
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 4;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.25 2.6]);
ylim([11.4 12.8]);
xlabel("b_1");
ylabel("z_{max}")
title('c');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gcf,'Renderer','zbuffer')
sz = 2;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.2 3.2]);
ylim([9.5 13]);
xlabel("b_1");
ylabel("z_{max}")
title('a');
ax = gca;
ax.TitleHorizontalAlignment = 'left';

% figure 4b
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 4;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.25 2.6]);
ylim([11.4 12.8]);
xlabel("b_1");
ylabel("z_{max}")
title('c');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gcf,'Renderer','zbuffer')
sz = 1;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([3 6.5]);
ylim([3 10]);
xlabel("b_1");
ylabel("z_{max}")
title('b');
ax = gca;
ax.TitleHorizontalAlignment = 'left';

% figure 4c
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 4;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.25 2.6]);
ylim([11.4 12.8]);
xlabel("b_1");
ylabel("z_{max}")
title('c');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
