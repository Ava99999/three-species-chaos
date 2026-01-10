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
ic_x = rand(param.ic_r,1);
ic_y = rand(param.ic_r,1);
ic_z = 7 + 4.*rand(param.ic_r,1);
ic_range = [ic_x ic_y ic_z];
b_par = param.b2_r; % choice for b1 or b2


%% bifurcation
Z = zeros(param.NB, length(b_par),size(ic_range,1)); % record timestep values
Z_locmax = zeros(param.NB, length(b_par),size(ic_range,1)); % record local maxima
Z_max = zeros(length(b_par),size(ic_range,1)); % record global maxima
Z_min = zeros(length(b_par),size(ic_range,1)); % record global minima
errors = zeros(length(b_par),size(ic_range,1)); % in case no cycling behavior is found
margin = 0.66;

tic
parfor j = 1:length(b_par)

    Z_j = zeros(param.NB, length(ic_x));
    Zmax_j = zeros(1, length(ic_x));
    Zmin_j = zeros(1, length(ic_x));

    for k = 1:length(ic_x)
        f1_loc = @(u) param.a1*u/(1 + param.b1_r(j)*u);
        f2_loc = @(u) param.a2*u/(1 + param.b2_r(j)*u);

        [t,X] = ode45(@(t,y) superpredators( ...
            t,y,f1,f2_loc,param.d1,param.d2), ...
            tspan, ic_range(k,:), options);

        if size(X,1) > param.NB
            tail = X(end-param.NB+1:end,3);
            Z_j(:,k) = tail;
            Zmax_j(k) = max(tail);
            Zmin_j(k) = min(tail);
        end
    end

    Z(:,j,:)   = Z_j;
    Z_max(j,:) = Zmax_j;
    Z_min(j,:) = Zmin_j;
end
toc

tic
% select the local maxima that are large enough compared to the cycle
for k = 1:length(ic_x)
    for j = 1:length(b_par)
       for i = 2:(param.NB-1) % todo does the first index matter?
           % identify maxima
           if Z(i,j,k) >= Z(i-1,j,k) && Z(i,j,k) >= Z(i+1,j,k)
                Z_range = Z_max(j,k) - Z_min(j,k);
                % only select those values that are large enough
                if Z(i,j,k) > (Z_min(j,k) + margin*Z_range)
                    Z_locmax(i,j,k) = Z(i,j,k);
                    
                    % temp check
                    %{
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
                    %}
                end
           end
       end
    end
end
toc

%% plot results
B1_layer = repelem(b_par, 1, param.NB); 
N_ic = size(ic_range, 1);
B1_flat = repmat(B1_layer, 1, N_ic);
Z_flat = reshape(Z_locmax, 1, []);
plot_matrix = [B1_flat', Z_flat'];
plot_matrix = plot_matrix(plot_matrix(:,2) ~= 0, :);
plot_points = unique(plot_matrix, 'rows', 'stable'); 


% temporary plot
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 0.5;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([1.5 2.15]);
%ylim([0 12]);
xlabel("b_2");
ylabel("z_{max}")
%title('a');
%ax = gca;
%ax.TitleHorizontalAlignment = 'left';


% Plot of all points
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 1;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.2 6.5]);
ylim([3 13]);
xlabel("b_1");
ylabel("z_{max}");

% Figure 4a
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 0.5;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.2 3.2]);
ylim([9.5 13]);
xlabel("b_1");
ylabel("z_{max}");
title('a');
ax = gca;
ax.TitleHorizontalAlignment = 'left';

% Figure 4b
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 1;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([3.2 6.5]);
ylim([3 10]);
xlabel("b_1");
ylabel("z_{max}");
title('b');
ax = gca;
ax.TitleHorizontalAlignment = 'left';

% Figure 4c
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 1;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.25 2.6]);
ylim([11.4 12.8]);
xlabel("b_1");
ylabel("z_{max}");
title('c');
ax = gca;
ax.TitleHorizontalAlignment = 'left';

%{
%% plot results
% copies of the parameter range, long row vector
% 2d version (1 ic)
B1 = repelem(param.b1_r,1,param.NB); 
Z_temp = reshape(Z_locmax,1,[]);
plot_matrix = [B1; Z_temp];
% remove duplicate maxima/minima
plot_points = unique(plot_matrix',"rows", 'stable'); 
indices = find(plot_points(:,2)==0);
plot_points(indices,:) = []; 

Z_flat = Z_locmax(:);
nb1 = length(param.b1_r);
nk  = size(ic_range,1);
B1_flat = repelem(param.b1_r, param.NB * nk).';
mask = Z_flat ~= 0;
plot_matrix = [B1_flat(mask), Z_flat(mask)];
plot_points = unique(round(plot_matrix,6), 'rows', 'stable');




% plot of all points
figure('Units','inches','Position',[1 1 6 4]);
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
sz = 1;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([3.2 6.5]);
ylim([3 10]);
xlabel("b_1");
ylabel("z_{max}")
title('b');
ax = gca;
ax.TitleHorizontalAlignment = 'left';

% figure 4c
figure('Units','inches','Position',[1 1 6 4]);
set(gcf,'Renderer','zbuffer')
sz = 1;
scatter(plot_points(:,1), plot_points(:,2), sz,'black','filled');
xlim([2.25 2.6]);
ylim([11.4 12.8]);
xlabel("b_1");
ylabel("z_{max}")
title('c');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
%}