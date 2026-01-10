clear all
close all

chaos_parameters

%% defining three species system
f1 = @(u) param.a1*u/(1 + param.b1*u); % functional response f1
f2 = @(u) param.a2*u/(1 + param.b2*u); % functional response f2

%% solve system in time
t0 = 0;
T = 10000;
tspan = param.t; %[t0,T];
ic = param.ic3; % initial condition from reprod. paper
ic_pert = ic + [0.01, 0, 0];
options = odeset('RelTol',1e-8,'abstol',1e-8*[1,1,1]);
[t,X] = ode45(@(t,y) superpredators(t,y,f1,f2,param.d1,param.d2),tspan, ic, options);

% illustrate dependency ic's
[t2,X2] = ode45(@(t,y) superpredators(t,y,f1,f2,param.d1,param.d2),tspan, ic_pert, options);

%% time series plots
ts = 4999;
te = 6499;
%ts = 1999;
%te = 9499;

figure;
plot(t(ts:te),X(ts:te,1));
xlabel('t');
ylabel('x');
axis tight;
ylim([0.1 1]);

figure;
plot(t(ts:te),X(ts:te,2));
xlabel('t');
ylabel('y');
axis tight;
ylim([0 0.45]);

figure;
plot(t(ts:te),X(ts:te,3));
xlabel('t');
ylabel('z');
axis tight;
ylim([7.2 10.7]);

%% 3d plot
figure;
p = plot3(X(ts:te,1), X(ts:te,2), X(ts:te,3));
xlabel('x');
ylabel('y');
zlabel('z');
p.Color = "black";
grid on;
set(gcf, 'Position', [100 100 600 600])  % [left bottom width height]
title(sprintf('b_1 = %.4g', param.b1));
axis tight;

view([45 30]);
set(gca, 'XDir', 'normal');
set(gca, 'YDir', 'reverse');


%% perturbation plot time series
figure;
plot(t(1:501),X(1:501,1));
hold on;
plot(t2(1:501),X2(1:501,1));

xlabel('t');
ylabel('x');
axis tight;
ylim([0.1 1]);
hold off

% visualize b_1*x and b_2*y
mean_b1 = mean(param.b1*X(ts:te,1));
mean_b2 = mean(param.b2*X(ts:te,2));
max_b1 = max(param.b1*X(ts:te,1));
min_b1 = min(param.b1*X(ts:te,1));
max_b2 = max(param.b2*X(ts:te,2));
min_b2 = min(param.b2*X(ts:te,2));

figure;
plot(t(ts:te), param.b1 * X(ts:te,1), 'r', 'DisplayName', 'b_1*x');
hold on;
plot(t(ts:te), mean_b1 * ones(size(t(ts:te))), 'r--', 'DisplayName', 'mean b_1*x')
plot(t(ts:te), param.b2 * X(ts:te,2), 'b', 'DisplayName', 'b_2*y');
plot(t(ts:te), mean_b2 * ones(size(t(ts:te))), 'b--', 'DisplayName', 'mean b_2*x')
xlabel('t');
%ylabel('Population');
legend show;
axis tight;
hold off;

% visualize f1 and f2
figure;
xplot = linspace(0,1);
plot(xplot, f1_plot(xplot,param.a1,param.b1), 'r', 'DisplayName', 'f_1');
hold on
plot(xplot, param.b1*xplot, 'g', 'DisplayName', 'b_1 * x')
plot(xplot, f2_plot(xplot,param.a2,param.b2), 'b', 'DisplayName', 'f_2');
xlabel('x');
%ylabel('Population');
legend show;
axis tight;
hold off;

figure;
plot(param.b1*xplot, f1_plot(xplot, param.a1,param.b1), 'LineWidth', 1);
ylabel('f_1(x)')
xlabel('b_1*x');

figure;
yplot = linspace(0,0.5);
plot(param.b2*yplot, f2_plot(yplot, param.a2,param.b2), 'LineWidth', 1);
ylabel('f_2(y)')
xlabel('b_2*y');

%{
% PARAMETERS
t_end = min(20000, te);
step  = 1;
fps   = 30;
outfile = 'trajectory_3d.mp4';

idx = ts:step:t_end;

% FIXED AXIS LIMITS (computed once)
xlims = [min(X(idx,1)) max(X(idx,1))];
ylims = [min(X(idx,2)) max(X(idx,2))];
zlims = [min(X(idx,3)) max(X(idx,3))];

% VIDEO SETUP
v = VideoWriter(outfile, 'MPEG-4');
v.FrameRate = fps;
open(v);

% FIGURE SETUP
fig = figure('Position', [100 100 600 600]);
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');

xlabel(ax, 'x');
ylabel(ax, 'y');
zlabel(ax, 'z');

xlim(ax, xlims);
ylim(ax, ylims);
zlim(ax, zlims);

view(ax, [45 30]);
set(ax, 'XDir', 'normal', 'YDir', 'reverse');

% PREALLOCATED LINE (FAST)
p = plot3(ax, NaN, NaN, NaN, 'k', 'LineWidth', 1);

% OPTIONAL: moving point
h = plot3(ax, NaN, NaN, NaN, 'ro', 'MarkerFaceColor','r');

% ANIMATION LOOP
for k = 1:numel(idx)
    i = idx(1:k);

    p.XData = X(i,1);
    p.YData = X(i,2);
    p.ZData = X(i,3);

    h.XData = X(idx(k),1);
    h.YData = X(idx(k),2);
    h.ZData = X(idx(k),3);

    drawnow limitrate
    writeVideo(v, getframe(fig));
end

close(v);
%}

function y = f1_plot(u,a1,b1)
    y = a1.*u./(1 + b1.*u); 
end 

function y = f2_plot(u,a2,b2)
    y = a2.*u./(1 + b2.*u); 
end 